/*
 * Author: Benoit Sklenard benoit.sklenard@cea.fr 
 * 
 * 3D Poisson solver.
 * Neumann and Dirichlet boundary conditions are supported and read from Kernel::Mesh class.
 *
 * Matrix operations rely on boost::numeric::ublas library. Atlas library should be more appropriate.
 *
 * Conjugate gradient is used to solve linear system Ax = b
 * Newton-Raphson method is used to get self-consistent Thomas-Fermi/Poisson solution
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <iomanip>

#include "kernel/Domain.h" 
#include "kernel/Mesh.h"
#include "kernel/Constants.h"

#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include "io/FileParameters.h"
#include "io/Parameters.h"

#include "okmc/MobileParticleParam.h"

#include "lkmc/LatticeDiamondParam.h"

#include "FermiDirac.h"
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/timer.hpp>

#include "Poisson.h"

using namespace Kernel;
using namespace boost::numeric;
using namespace Electrostatics;


Poisson::Poisson(Tcl_Interp *pTcl, Domain *pDomain) {
	LOWMSG("Loading Poisson...");

	_pDomain            = pDomain;
	_pMesh              = pDomain->_pMesh;
	_pTcl               = pTcl;

	_T                  = 0;

	_nx                 = (_pMesh->getPeriodicX() ? _pMesh->getnx() - 1 : _pMesh->getnx());
	_ny                 = (_pMesh->getPeriodicY() ? _pMesh->getny() - 1 : _pMesh->getny());
	_nz                 = (_pMesh->getPeriodicZ() ? _pMesh->getnz() - 1 : _pMesh->getnz());

	_CG_tol             = Domains::global()->getFileParameters()->getFloat("MC/Electrostatic/cg.tolerance");       // 1e-6;
	_CG_maxIter         = Domains::global()->getFileParameters()->getInt("MC/Electrostatic/cg.max.iteration");     // 2000;
	_CG_verbose         = Domains::global()->getFileParameters()->getInt("MC/Electrostatic/cg.verbose");

	_NR_tol             = Domains::global()->getFileParameters()->getFloat("MC/Electrostatic/newton.tolerance");   // 1e-4;
	_NR_maxIter         = Domains::global()->getFileParameters()->getInt("MC/Electrostatic/newton.max.iteration"); // 50;
	_NR_verbose         = Domains::global()->getFileParameters()->getInt("MC/Electrostatic/newton.verbose");

	_FermiDirac         = Domains::global()->getFileParameters()->getBool("MC/Electrostatic/FermiDirac");
	_partialIonization  = Domains::global()->getFileParameters()->getBool("MC/Electrostatic/partial.ionization");
	_selfConsistentLKMC = Domains::global()->getFileParameters()->getBool("MC/Electrostatic/self.consistent.lkmc");

	_Qe                 = ublas::zero_vector<double>(_nx * _ny * _nz);
	_Qh                 = ublas::zero_vector<double>(_nx * _ny * _nz);
	_Qokmc              = ublas::zero_vector<double>(_nx * _ny * _nz);
	_psi                = ublas::zero_vector<double>(_nx * _ny * _nz);

	_Qlkmc_M            = ublas::zero_vector<double>(_nx * _ny * _nz);
	_Qlkmc_P            = ublas::zero_vector<double>(_nx * _ny * _nz);

	_Eg                 = ublas::zero_vector<double>(_nx * _ny * _nz);
	_Eg0                = ublas::zero_vector<double>(_nx * _ny * _nz);
	_Nc                 = ublas::zero_vector<double>(_nx * _ny * _nz);
	_Nv                 = ublas::zero_vector<double>(_nx * _ny * _nz);
	_DEc                = ublas::zero_vector<double>(_nx * _ny * _nz);
	_DEv                = ublas::zero_vector<double>(_nx * _ny * _nz);

	const IO::ParameterManager      *pPM  = Domains::global()->PM();
	const OKMC::MobileParticleParam *pMPP = _pDomain->_pMPPar;
	const LKMC::LatticeParam             **pLP  = _pDomain->_pLaPar;

	for(M_TYPE mt = 0; mt < pPM->getNMaterials(); ++mt)
	{
		// OKMC
		for(P_TYPE pt = 0; pt < pPM->getNParticles(); ++pt)
		{
			if(!pPM->isParticleDefined(pt, mt))
				continue;
			if (pMPP->_mapToGrid[mt][pt])
				_Nokmc[mt][pt] = ublas::zero_vector<double>(_nx * _ny * _nz);
		}
		// LKMC
		if (pLP[mt] != NULL && pLP[mt]->_mapToGrid)
			_Nlkmc[mt] = ublas::zero_vector<double>(_nx * _ny * _nz);
	}
}

Poisson::~Poisson() {
}

void Poisson::computeElectronicParameters(double TKelvin) {

	MeshNode           ***pNode = _pMesh->getNodes();
	IO::FileParameters   *pPar  = Domains::global()->getFileParameters();
	IO::ParameterManager *pPM   = Domains::global()->PM();

	for (Kernel::M_TYPE mt = 0; mt < Domains::global()->PM()->getNMaterials(); ++mt) {
		std::string path = pPM->getMaterialName(mt) + "/ElectronicStructure";

		if (!pPar->specified(path)) {
			WARNINGMSG(path << " does not exist!");
			_matParam[mt]._Nc  = 0;
			_matParam[mt]._Nv  = 0;
			_matParam[mt]._Eg  = 0;
			_matParam[mt]._Eg0 = 0;

			continue ;
		}
		std::string eDOSPath           = path + "/eDOSMass";
		std::string hDOSPath           = path + "/hDOSMass";
		std::string EgPath             = path + "/Bandgap";
		std::string EcDilatationalPath = path + "/Ec.dilatational";
		std::string EvDilatationalPath = path + "/Ev.dilatational";
		std::string EcDeviatoricPath   = path + "/Ec.deviatoric";
		std::string EvDeviatoricPath   = path + "/Ev.deviatoric";

		std::string T  = boost::lexical_cast<std::string>(_T);
		std::string T0 = boost::lexical_cast<std::string>(300);

		pPar->loadProcedure(_pTcl, eDOSPath, 1);
		pPar->loadProcedure(_pTcl, hDOSPath, 1);
		pPar->loadProcedure(_pTcl, EgPath, 1);

		double mdc = pPar->getFloatProc(_pTcl, eDOSPath, T);
		double mdv = pPar->getFloatProc(_pTcl, hDOSPath, T);

		_matParam[mt]._Nc  = 2. * pow(M0 * mdc * KB * _T / (2 * M_PI * pow(PLANCK_BAR, 2)), 3./2.);  // m3
		_matParam[mt]._Nv  = 2. * pow(M0 * mdv * KB * _T / (2 * M_PI * pow(PLANCK_BAR, 2)), 3./2.);  // m3
		_matParam[mt]._Eg  = ELECTRONVOLT_TO_HARTREE(pPar->getFloatProc(_pTcl, EgPath, T));     // Ha
		_matParam[mt]._Eg0 = ELECTRONVOLT_TO_HARTREE(pPar->getFloatProc(_pTcl, EgPath, T0));    // Ha

		if (pPar->specified(EcDilatationalPath) && pPar->specified(EvDilatationalPath)) {
			_matParam[mt]._Dc  = ublas::zero_vector<double>(3);
			_matParam[mt]._Dcx = ublas::zero_vector<double>(3);
			_matParam[mt]._Dcy = ublas::zero_vector<double>(3);
			_matParam[mt]._Dcz = ublas::zero_vector<double>(3);
			_matParam[mt]._Dv  = ublas::zero_vector<double>(2);
			_matParam[mt]._Dvb = ublas::zero_vector<double>(2);
			_matParam[mt]._Dvd = ublas::zero_vector<double>(2);

			std::map<std::string, float> Dc  = pPar->getFloatMap(EcDilatationalPath);
			std::map<std::string, float> Dcx = pPar->getFloatMap(EcDeviatoricPath + "(1)");
			std::map<std::string, float> Dcy = pPar->getFloatMap(EcDeviatoricPath + "(2)");
			std::map<std::string, float> Dcz = pPar->getFloatMap(EcDeviatoricPath + "(3)");

			std::map<std::string, float> Dv  = pPar->getFloatMap(EvDilatationalPath);
			std::map<std::string, float> Dvb = pPar->getFloatMap(EvDeviatoricPath + "(1)");
			std::map<std::string, float> Dvd = pPar->getFloatMap(EvDeviatoricPath + "(2)");

			for (unsigned i = 0; i < 3; ++i) {
				_matParam[mt]._Dc(i)  = Dc[boost::lexical_cast<std::string>(i + 1)];
				_matParam[mt]._Dcx(i) = Dcx[boost::lexical_cast<std::string>(i + 1)];
				_matParam[mt]._Dcy(i) = Dcy[boost::lexical_cast<std::string>(i + 1)];
				_matParam[mt]._Dcz(i) = Dcz[boost::lexical_cast<std::string>(i + 1)];
			}

			for (unsigned i = 0; i < 2; ++i) {
				_matParam[mt]._Dv(i)  = Dv[boost::lexical_cast<std::string>(i + 1)];
				_matParam[mt]._Dvb(i) = Dvb[boost::lexical_cast<std::string>(i + 1)];
				_matParam[mt]._Dvd(i) = Dvd[boost::lexical_cast<std::string>(i + 1)];
			}
			LOWMSG("got stress data for " << pPM->getMaterialName(mt));
		}
	}

	// initialize electronic parameters for the structure
	for (unsigned ix = 0; ix < _nx; ++ix) {
		for (unsigned iy = 0; iy < _ny; ++iy) {
			for (unsigned iz = 0; iz < _nz; ++iz) {
				if (pNode[ix][iy][iz]._active) {
					const unsigned idx = pNode[ix][iy][iz]._index;

					_Nc(idx)  = pNode[ix][iy][iz]._volume * 1e-27 * getNc(&pNode[ix][iy][iz]);
					_Nv(idx)  = pNode[ix][iy][iz]._volume * 1e-27 * getNv(&pNode[ix][iy][iz]);

					_Eg(idx)  = getEg(&pNode[ix][iy][iz]);
					_Eg0(idx) = getEg0(&pNode[ix][iy][iz]);

					if (_Eg(pNode[ix][iy][iz]._index) <= 0. || _Eg0(pNode[ix][iy][iz]._index) <= 0.)
						ERRORMSG("Poisson::computeElectronicParameters: band gap cannot be zero!");
				}
			}
		}
	}
}

void Poisson::setFirstGuess() {
	MeshNode ***pNode                     = _pMesh->getNodes();
	const IO::ParameterManager *pPM       = Domains::global()->PM();
	const OKMC::MobileParticleParam *pMPP = _pDomain->_pMPPar;

	for (size_t ix = 0; ix < _nx; ++ix) {
		for (size_t iy = 0; iy < _ny; ++iy) {
			for (size_t iz = 0; iz < _nz; ++iz) {
				if (pNode[ix][iy][iz]._active) {
					const unsigned i = pNode[ix][iy][iz]._index;
					double n = 0;

					for(M_TYPE mt = 0; mt < pPM->getNMaterials(); ++mt)
						for(P_TYPE pt = pPM->getNFamilies(); pt < pPM->getNParticles(); ++pt)
						{
							P_TYPE fam = pPM->getFamily(pt);
							P_POS  pos = pPM->getPPos(pt);
							if(fam == V_TYPE || fam == pPM->getMaterial(mt)._pt[0] || fam == pPM->getMaterial(mt)._pt[1] ||
								pos == POS_I || pos == POS_V)
								continue;
							if(!_Nokmc[mt][pt].empty() && _Nokmc[mt][pt](i) > 0)
								n += pMPP->_state2charge[mt][pt][0]; //impurities and dopants only have state 0!
						}

					if (n > 0.)
						_psi(i) = 0.5 * getEg(&pNode[ix][iy][iz]);
					else if (n < 0.)
						_psi(i) = - 0.5 * getEg(&pNode[ix][iy][iz]);
					else
						_psi(pNode[ix][iy][iz]._index) = 0;

					for(M_TYPE mt = 0; mt < pPM->getNMaterials(); ++mt) {
						if (!_Nlkmc[mt].empty() && _Nlkmc[mt](i) > 0) {
							_psi(i) = 0;
						}
					}
				}
			}
		}
	}
}

void Poisson::setT(double T)
{
	if (_T != T)
	{
		_T = T;
		MEDMSG("Building vector for OKMC particles...");
		buildOKMCVector();

		if (_selfConsistentLKMC)
			buildLKMCVector();

		MEDMSG("Building Laplacian matrix...");
		buildLaplacianMatrix();
		computeElectronicParameters(T);
		setFirstGuess();
	}
}

bool Poisson::compute() {
	MEDMSG("Building vector for OKMC particles...");
	buildOKMCVector();

	if (_selfConsistentLKMC)
		buildLKMCVector();

	MEDMSG("Solving non-linear Poisson equation...");

	if (! NewtonRaphson())
		ERRORMSG("Poisson::NewtonRaphson did not converged!");

	map2Nodes();
	map2Elements();

	return true;
}

double Poisson::getDEc(MeshNode *pNode) {
	const double kT = KB / Q * _T;

	std::set<MeshElement *> s;

	double DEc = 0;
	_pMesh->getElementsFromNode(pNode, s);

	for (std::set<MeshElement *>::const_iterator it = s.begin(); it != s.end(); ++it) {
		const M_TYPE mt = (*it)->getMaterial();

		// ublas::vector<double> DEci = ublas::zero_vector<double>(3);

		double tmp = 0;
		for (unsigned i = 0; i < 3; ++i) {
			// const ublas::vector<double> E = (*it)->strain();
			ublas::vector<double> E = ublas::zero_vector<double>(3);

			const double Dc         = (_matParam[mt]._Dc.size()  > 0 ? _matParam[mt]._Dc(i)  : 0);
			const double Dcx        = (_matParam[mt]._Dcx.size() > 0 ? _matParam[mt]._Dcx(i) : 0);
			const double Dcy        = (_matParam[mt]._Dcy.size() > 0 ? _matParam[mt]._Dcy(i) : 0);
			const double Dcz        = (_matParam[mt]._Dcz.size() > 0 ? _matParam[mt]._Dcz(i) : 0);

			tmp += exp(- (Dc * (E(0) + E(1) + E(2)) + Dcx * E(0) + Dcy * E(1) + Dcz * E(2))/ kT);
		}

		DEc += ELECTRONVOLT_TO_HARTREE(-kT * log(1./3. * tmp));

		if (Domains::global()->PM()->getMaterialName(mt) == "Silicon")
			LOWMSG("DEc:" << DEc / static_cast<double>(s.size()));
	}

	return (DEc / static_cast<double>(s.size()));
}

double Poisson::getNc(MeshNode *pNode) {
	std::set<MeshElement *> s;
	double Nc = 0;
	_pMesh->getElementsFromNode(pNode, s);

	for (std::set<MeshElement *>::const_iterator it = s.begin(); it != s.end(); ++it)
		Nc += _matParam[(*it)->getMaterial()]._Nc;

	return (Nc / static_cast<double>(s.size()));
}

double Poisson::getNv(MeshNode *pNode) {
	std::set<MeshElement *> s;
	double Nv = 0;
	_pMesh->getElementsFromNode(pNode, s);

	for (std::set<MeshElement *>::const_iterator it = s.begin(); it != s.end(); ++it)
		Nv += _matParam[(*it)->getMaterial()]._Nv;

	return (Nv / static_cast<double>(s.size()));
}

double Poisson::getEg(MeshNode *pNode) {
	std::set<MeshElement *> s;
	double Eg = 0;
	_pMesh->getElementsFromNode(pNode, s);

	for (std::set<MeshElement *>::const_iterator it = s.begin(); it != s.end(); ++it)
		Eg += _matParam[(*it)->getMaterial()]._Eg;

	return (Eg / static_cast<double>(s.size()));
}

double Poisson::getEg0(MeshNode *pNode) {
	std::set<MeshElement *> s;
	double Eg0 = 0;
	_pMesh->getElementsFromNode(pNode, s);

	for (std::set<MeshElement *>::const_iterator it = s.begin(); it != s.end(); ++it)
		Eg0 += _matParam[(*it)->getMaterial()]._Eg0;

	return (Eg0 / static_cast<double>(s.size()));
}

// OBSOLETE METHOD
// Do NOT call!
/* double Poisson::getLA(MeshNode *pNode) {
	double w = 0.;

	for (std::map<LKMC::LatticeAtom *, double>::iterator it = pNode->_mLA.begin(); it !=pNode->_mLA.end(); ++it)
		w += it->second;

	return w;
} */

void Poisson::map2Nodes() {
	MeshNode ***pNode = _pMesh->getNodes();

	for (unsigned ix = 0; ix < _nx; ++ix) {
		for (unsigned iy = 0; iy < _ny; ++iy) {
			for (unsigned iz = 0; iz < _nz; ++iz) {
				if (pNode[ix][iy][iz]._active) {
					pNode[ix][iy][iz]._V        = HARTREE_TO_ELECTRONVOLT(_psi(pNode[ix][iy][iz]._index));
					pNode[ix][iy][iz]._Eg       = HARTREE_TO_ELECTRONVOLT(_Eg(pNode[ix][iy][iz]._index));
					pNode[ix][iy][iz]._Eg0      = HARTREE_TO_ELECTRONVOLT(_Eg0(pNode[ix][iy][iz]._index));
					pNode[ix][iy][iz]._Ec       = HARTREE_TO_ELECTRONVOLT(0.5 * getEg(&pNode[ix][iy][iz]) - _psi(pNode[ix][iy][iz]._index));
					pNode[ix][iy][iz]._Ev       = HARTREE_TO_ELECTRONVOLT(- 0.5 * getEg(&pNode[ix][iy][iz]) - _psi(pNode[ix][iy][iz]._index));
					pNode[ix][iy][iz]._eDensity = _Qe(pNode[ix][iy][iz]._index) / (pNode[ix][iy][iz]._volume * 1e-21);
					pNode[ix][iy][iz]._hDensity = _Qh(pNode[ix][iy][iz]._index) / (pNode[ix][iy][iz]._volume * 1e-21);

					// pNode[ix][iy][iz]._N_okmc   = _Qokmc(pNode[ix][iy][iz]._index) / (pNode[ix][iy][iz]._volume * 1e-27); // FIXME

					if (ix > 0 && ix < (_nx - 1) && iy > 0 && iy < (_ny - 1) && iz > 0 && iz < (_nz - 1)) {
						if (pNode[ix-1][iy][iz]._index < 0)
							continue ;
						if (pNode[ix+1][iy][iz]._index < 0)
							continue ;
						if (pNode[ix][iy-1][iz]._index < 0)
							continue ;
						if (pNode[ix][iy+1][iz]._index < 0)
							continue ;
						if (pNode[ix][iy][iz-1]._index < 0)
							continue ;
						if (pNode[ix][iy][iz+1]._index < 0)
							continue ;

						const double dxp = pNode[ix+1][iy][iz]._coord._x - pNode[ix][iy][iz]._coord._x;
						const double dxm = pNode[ix][iy][iz]._coord._x - pNode[ix-1][iy][iz]._coord._x;
						const double dyp = pNode[ix][iy+1][iz]._coord._y - pNode[ix][iy][iz]._coord._y;
						const double dym = pNode[ix][iy][iz]._coord._y - pNode[ix][iy-1][iz]._coord._y;
						const double dzp = pNode[ix][iy][iz+1]._coord._z - pNode[ix][iy][iz]._coord._z;
						const double dzm = pNode[ix][iy][iz]._coord._z - pNode[ix][iy][iz-1]._coord._z;

						const double Vxp = HARTREE_TO_ELECTRONVOLT(_psi(pNode[ix+1][iy][iz]._index));
						const double Vxm = HARTREE_TO_ELECTRONVOLT(_psi(pNode[ix-1][iy][iz]._index));
						const double Vyp = HARTREE_TO_ELECTRONVOLT(_psi(pNode[ix][iy+1][iz]._index));
						const double Vym = HARTREE_TO_ELECTRONVOLT(_psi(pNode[ix][iy-1][iz]._index));
						const double Vzp = HARTREE_TO_ELECTRONVOLT(_psi(pNode[ix][iy][iz+1]._index));
						const double Vzm = HARTREE_TO_ELECTRONVOLT(_psi(pNode[ix][iy][iz-1]._index));

						const double V   = HARTREE_TO_ELECTRONVOLT(_psi(pNode[ix][iy][iz]._index));

						pNode[ix][iy][iz]._E    = ublas::zero_vector<double>(3); // electric field in eV/nm
						pNode[ix][iy][iz]._E(0) = - 0.5 * ((Vxp - V) / dxp + (V  - Vxm) / dxm);
						pNode[ix][iy][iz]._E(1) = - 0.5 * ((Vyp - V) / dyp + (V  - Vym) / dym);
						pNode[ix][iy][iz]._E(2) = - 0.5 * ((Vzp - V) / dzp + (V  - Vzm) / dzm);
					}
				}
			}
		}
	}
}

void Poisson::map2Elements() {
	for(Kernel::Mesh::iterator it = _pMesh->begin(); it != _pMesh->end(); ++it)
	{
		std::set<MeshNode *> s;
		const M_TYPE mt = it->getMaterial();

		_pMesh->getNodesFromElement(&(*it), s);

		if (s.empty())
			ERRORMSG("Poisson::map2Elements: cannot find nodes of element " << it->getIndex());

		double V  = 0;
		for (std::set<MeshNode *>::const_iterator i = s.begin(); i != s.end(); ++i)
			V  += (*i)->_V;

		it->electrostaticPotential() = V / static_cast<double>(s.size());
		it->bandGap() = HARTREE_TO_ELECTRONVOLT(_matParam[mt]._Eg);
	}
}

bool Poisson::NewtonRaphson() {
	MeshNode ***pNode                     = _pMesh->getNodes();
	const double kT                       = ELECTRONVOLT_TO_HARTREE(KB / Q * _T);
	const IO::ParameterManager *pPM       = Domains::global()->PM();
	const OKMC::MobileParticleParam *pMPP = _pDomain->_pMPPar;
	bool         ret                      = false;

	ublas::vector<double> deltaPsi  = ublas::zero_vector<double>(_nx * _ny * _nz);
	ublas::vector<double> err       = ublas::zero_vector<double>(_nx * _ny * _nz);
	ublas::vector<double> totCharge = ublas::zero_vector<double>(_nx * _ny * _nz);

	// ublas::vector<double> Qtmp      = ublas::zero_vector<double>(_nx * _ny * _nz);

	boost::timer t;
	t.restart();

	for (size_t it = 0; it < _NR_maxIter; ++it) {
		ublas::compressed_matrix<double> J = _laplacian;

		for (size_t ix = 0; ix < _nx; ++ix) {
			for (size_t iy = 0; iy < _ny; ++iy) {
				for (size_t iz = 0; iz < _nz; ++iz) {
					if (pNode[ix][iy][iz]._active) {
						const size_t i   = pNode[ix][iy][iz]._index;

						double dQe        = 0.;
						double dQh        = 0.;
						double dQokmc     = 0;
						double dQlkmc_M   = 0.;
						double dQlkmc_P   = 0.;

						const double Ef = 0;
						const double CBM = 0.5 * _Eg(i) - _psi(i);
						const double VBM = - 0.5 * _Eg(i) - _psi(i);

						const double EfEc = Ef - CBM;
						const double EvEf = VBM - Ef;
						// const double EfEc = _psi(i) - 0.5 * _Eg(i);
						// const double EvEf = - _psi(i) - 0.5 * _Eg(i);

						_Qokmc(i)   = 0;
						_Qlkmc_M(i) = 0.;
						_Qlkmc_P(i) = 0.;

						if (_FermiDirac) {
							dQe    = _Nc(i) / kT * FermiDirac::mhalf(EfEc / kT);
							dQh    = - _Nv(i) / kT * FermiDirac::mhalf(EvEf / kT);
							_Qe(i) = _Nc(i) * FermiDirac::phalf(EfEc / kT);
							_Qh(i) = _Nv(i) * FermiDirac::phalf(EvEf / kT);
						}
						else {
							dQe    = _Nc(i) / kT * exp(EfEc / kT);
							dQh    = - _Nv(i) / kT * exp(EvEf / kT);
							_Qe(i) = _Nc(i) * exp(EfEc / kT);
							_Qh(i) = _Nv(i) * exp(EvEf / kT);
						}

						for(M_TYPE mt = 0; mt < pPM->getNMaterials(); ++mt)
						{// OKMC
							for(P_TYPE pt = pPM->getNFamilies(); pt < pPM->getNParticles(); ++pt)
							{
								P_TYPE fam = pPM->getFamily(pt);
								P_POS  pos = pPM->getPPos(pt);
								if(fam == V_TYPE || fam == pPM->getMaterial(mt)._pt[0] || fam == pPM->getMaterial(mt)._pt[1] ||
										pos == POS_I || pos == POS_V)
									continue;
								if (! _Nokmc[mt][pt].empty() && _Nokmc[mt][pt](i) > 0)
								{
									if (pMPP->_state2charge[mt][pt][0] == 1) { //DONOR
										if (_partialIonization) {
											// const double E = _psi(i) - 0.5 * _Eg(i) + _Eg(i) / _Eg0(i) * ELECTRONVOLT_TO_HARTREE(pMPP->_stateEnergy[mt][pt]);
											const double E = Ef - (CBM - _Eg(i) / _Eg0(i) * ELECTRONVOLT_TO_HARTREE(pMPP->_stateEnergy[mt][pt]));
											_Qokmc(i) += _Nokmc[mt][pt](i) / (1. + pMPP->_stateDegeneracy[mt][pt] * exp(E / kT));
											dQokmc    += - _Nokmc[mt][pt](i) * pMPP->_stateDegeneracy[mt][pt] / kT * exp(E / kT) / pow(1. + pMPP->_stateDegeneracy[mt][pt] * exp(E / kT), 2);
										}
										else {
											_Qokmc(i) += _Nokmc[mt][pt](i);
										}
									}
									else if (pMPP->_state2charge[mt][pt][0] == -1) { //ACCEPTOR
										if (_partialIonization) {
											// const double E = - _psi(i) - 0.5 * _Eg(i) + _Eg(i) / _Eg0(i) * ELECTRONVOLT_TO_HARTREE(pMPP->_stateEnergy[mt][pt]);
											const double E = (VBM + _Eg(i) / _Eg0(i) * ELECTRONVOLT_TO_HARTREE(pMPP->_stateEnergy[mt][pt])) - Ef;
											_Qokmc(i) -= _Nokmc[mt][pt](i) / (1. + pMPP->_stateDegeneracy[mt][pt] * exp(E / kT));
											dQokmc    -= _Nokmc[mt][pt](i) * pMPP->_stateDegeneracy[mt][pt] / kT * exp(E / kT) / pow(1. + pMPP->_stateDegeneracy[mt][pt] * exp(E / kT), 2);
										}
										else {
											_Qokmc(i) -= _Nokmc[mt][pt](i);
										}
									}
								}
							}

							// LKMC
							/*
							if (_selfConsistentGFLS && _Nlkmc[mt].empty() == false && _Nlkmc[mt](i) > 0.) {
								if (_pDomain->_pLat[mt]->getType() != LKMC::LatticeParam::DIAMOND_1) {
									ERRORMSG("Poisson::NewtonRaphson does not support " << pPM->getMaterialName(mt) << " lattice.");
								}
								LKMC::LatticeDiamondParam *pLDP = static_cast<LKMC::LatticeDiamondParam *>(_pDomain->_pLaPar[mt]);

								const double E_lkmc_M = _psi(i) - 0.5 * _Eg(i) + _Eg(i) / _Eg0(i) * ELECTRONVOLT_TO_HARTREE(pLDP->_E_M0);   // Ef - E(0,-)
								const double E_lkmc_P = - _psi(i) - 0.5 * _Eg(i) + _Eg(i) / _Eg0(i) * ELECTRONVOLT_TO_HARTREE(pLDP->_E_P0); // E(+,0) - Ef

								const double g_M = pLDP->_g_M;
								const double g_0 = pLDP->_g_0;
								const double g_P = pLDP->_g_P;

								dQlkmc_M    +=   _Nc(i) / _Nlkmc[mt](i) * 1. / kT * g_M / g_0 * exp(E_lkmc_M / kT);
								dQlkmc_P    += - _Nv(i) / _Nlkmc[mt](i) * 1. / kT * g_P / g_0 * exp(E_lkmc_P / kT);
								_Qlkmc_M(i) += _Nc(i) / _Nlkmc[mt](i) * g_M / g_0 * exp(E_lkmc_M / kT);
								_Qlkmc_P(i) += _Nv(i) / _Nlkmc[mt](i) * g_P / g_0 * exp(E_lkmc_P / kT);
							}

							if (_selfConsistentGFLS2 && _Nlkmc[mt].empty() == false && _Nlkmc[mt](i) > 0.) {
								if (_pDomain->_pLat[mt]->getType() != LKMC::LatticeParam::DIAMOND_1) {
									ERRORMSG("Poisson::NewtonRaphson does not support " << pPM->getMaterialName(mt) << " lattice.");
								}
								LKMC::LatticeDiamondParam *pLDP = static_cast<LKMC::LatticeDiamondParam *>(_pDomain->_pLaPar[mt]);

								// const double Ef  = 0;
								// const double VBM = - V - 0.5 * Eg;

							    // Note: defect level scales with bandgap
								const double E_lkmc_M    = Ef - (VBM + _Eg(i) / _Eg0(i) * ELECTRONVOLT_TO_HARTREE(pLDP->_E_M0));      // Ef  - E(0,-)
								const double E_lkmc_P    = (VBM + _Eg(i) / _Eg0(i) * ELECTRONVOLT_TO_HARTREE(pLDP->_E_P0)) - Ef; // Ef  - E(0,-)

								const double g_M = pLDP->_g_M;
								const double g_0 = pLDP->_g_0;
								const double g_P = pLDP->_g_P;

								const double den_M = 1. + (g_0 / g_M) * exp( - E_lkmc_M / kT) + (g_P / g_M) * exp((E_lkmc_P - E_lkmc_M) / kT);
								const double den_P = 1. + (g_0 / g_P) * exp( - E_lkmc_P / kT) + (g_M / g_P) * exp((E_lkmc_M - E_lkmc_P) / kT);

								dQlkmc_M    +=   _Nlkmc[mt](i) * (1. / kT * (g_0 / g_M) * exp( - E_lkmc_M / kT) + 2. / kT * (g_P / g_M) * exp((E_lkmc_P - E_lkmc_M)/ kT)) / (den_M * den_M);
								dQlkmc_P    += - _Nlkmc[mt](i) * (1. / kT * (g_0 / g_P) * exp( - E_lkmc_P / kT) + 2. / kT * (g_M / g_P) * exp((E_lkmc_M - E_lkmc_P)/ kT)) / (den_P * den_P);
								_Qlkmc_M(i) += _Nlkmc[mt](i) / den_M;
								_Qlkmc_P(i) += _Nlkmc[mt](i) / den_P;
							}*/

							if (_selfConsistentLKMC && _Nlkmc[mt].empty() == false && _Nlkmc[mt](i) > 0.) {
								if (_pDomain->_pLat[mt]->getType() != LKMC::LatticeParam::DIAMOND) {
									ERRORMSG("Poisson::NewtonRaphson does not support " << pPM->getMaterialName(mt) << " lattice.");
								}
								const LKMC::LatticeDiamondParam *pLDP = static_cast<const LKMC::LatticeDiamondParam *>(_pDomain->_pLaPar[mt]);

								const double E_lkmc_M = _psi(i) - 0.5 * _Eg(i) + _Eg(i) / _Eg0(i) * ELECTRONVOLT_TO_HARTREE(pLDP->_E_M0);   // Ef - E(0,-)
								const double E_lkmc_P = - _psi(i) - 0.5 * _Eg(i) + _Eg(i) / _Eg0(i) * ELECTRONVOLT_TO_HARTREE(pLDP->_E_P0); // E(+,0) - Ef

								const double g_M = pLDP->_g_M;
								const double g_0 = pLDP->_g_0;
								const double g_P = pLDP->_g_P;

								const double den_M = 1. + (g_0 / g_M) * exp( - E_lkmc_M / kT) + (g_P / g_M) * exp((E_lkmc_P - E_lkmc_M) / kT);
								const double den_P = 1. + (g_0 / g_P) * exp( - E_lkmc_P / kT) + (g_M / g_P) * exp((E_lkmc_M - E_lkmc_P) / kT);

								dQlkmc_M    +=   _Nlkmc[mt](i) * (1. / kT * (g_0 / g_M) * exp( - E_lkmc_M / kT) + 2. / kT * (g_P / g_M) * exp((E_lkmc_P - E_lkmc_M)/ kT)) / (den_M * den_M);
								dQlkmc_P    += - _Nlkmc[mt](i) * (1. / kT * (g_0 / g_P) * exp( - E_lkmc_P / kT) + 2. / kT * (g_M / g_P) * exp((E_lkmc_M - E_lkmc_P)/ kT)) / (den_P * den_P);
								_Qlkmc_M(i) += _Nlkmc[mt](i) / den_M;
								_Qlkmc_P(i) += _Nlkmc[mt](i) / den_P;
							}
						}

						J(i,i) -= 4. * M_PI * (dQh - dQe + dQokmc + dQlkmc_P - dQlkmc_M);
					}
				}
			}
		}

		ublas::noalias(totCharge) = _Qh - _Qe + _Qokmc + _Qeff + _Qlkmc_P - _Qlkmc_M;
		ublas::noalias(err)       = - ublas::prod(_laplacian, _psi) + 4 * M_PI * totCharge;

		if (ublas::norm_2(err) < _NR_tol) {
			ret = true;
			if (_NR_verbose > 0)
				LOWMSG("Poisson::NewtonRaphson converged in " << it << " iterations (" << t.elapsed() << " s)");

			break ;
		}
		if (_NR_verbose > 1)
			LOWMSG("Poisson::NewtonRaphson iteration " << it << " error= " << ublas::norm_2(err) << " (" << t.elapsed() << " s)");

		deltaPsi = ublas::zero_vector<double>(_nx * _ny * _nz);

		if (! ConjugateGradient(J, deltaPsi, err))
			ERRORMSG("Poisson::CG did not converged");

		_psi += deltaPsi;
	}

	return ret;
}

bool Poisson::ConjugateGradient(const ublas::compressed_matrix<double> &A, ublas::vector<double> &x, const ublas::vector<double> &b) {
	boost::timer t;
	t.restart();

	ublas::vector<double> a;

	bool                  ret = false;
	ublas::vector<double> p   = ublas::zero_vector<double>(x.size());
	ublas::vector<double> r   = ublas::zero_vector<double>(x.size());


	ublas::noalias(r) = b - ublas::prod(A, x);
	ublas::noalias(p) = r;

	for (size_t it = 0; it < _CG_maxIter; ++it) {
		if (ublas::norm_2(r) < _CG_tol) {

			if (_CG_verbose > 0)
				LOWMSG("Poisson::CG solved in " << it << " iterations (" << t.elapsed() << " s)");

			ret = true;
			break ;
		}
		if (_CG_verbose > 1) {
			LOWMSG("Poisson::CG iteration " << it << " error= " << ublas::norm_2(r) << " (" << t.elapsed() << " s)");
		}

		a                 = ublas::zero_vector<double>(x.size());
		ublas::noalias(a) = ublas::prod(A, p);

		double  lambda = ublas::inner_prod(r, p) / ublas::inner_prod(a, p);

		x += lambda * p;
		r -= lambda * a;

		p = r - (ublas::inner_prod(r, a) / ublas::inner_prod(a, p)) * p;
	}

	return (ret);
}

// This method is for test purpose only
void Poisson::testPN() {
	MeshNode ***pNode = _pMesh->getNodes();

	for (size_t ix = 0; ix < _nx; ++ix)
		for (size_t iy = 0; iy < _ny; ++iy)
			for (size_t iz = 0; iz < _nz; ++iz)
				if (pNode[ix][iy][iz]._active) {
					if (ix < _nx / 2)
						_Qokmc(pNode[ix][iy][iz]._index) = .01;
					else
						_Qokmc(pNode[ix][iy][iz]._index) = -.005;

					if (ix == 0 || ix == _nx - 1 || iy == 0 || iy == _ny - 1 || iz == 0 || iz == _nz - 1)
						_Qokmc(pNode[ix][iy][iz]._index) *= .5;
				}
}

void Poisson::buildOKMCVector() {
	MeshNode ***pNode = _pMesh->getNodes();
	const IO::ParameterManager *pPM       = Domains::global()->PM();

	for(M_TYPE mt = 0; mt < pPM->getNMaterials(); ++mt) {
		for(P_TYPE pt = 0; pt < pPM->getNParticles(); ++pt) {
			if (!_Nokmc[mt][pt].empty()) {
				_Nokmc[mt][pt] = ublas::zero_vector<double>(_nx * _ny * _nz);
			}
		}
	}

	for (size_t ix = 0; ix < _nx; ++ix) {
		for (size_t iy = 0; iy < _ny; ++iy) {
			for (size_t iz = 0; iz < _nz; ++iz) {
				if (pNode[ix][iy][iz]._active && !pNode[ix][iy][iz]._mPart.empty()) {
					for (std::map<OKMC::Particle *, double>::const_iterator it = pNode[ix][iy][iz]._mPart.begin(); it != pNode[ix][iy][iz]._mPart.end(); ++it) {
						const OKMC::Particle * pPart = it->first;
						const M_TYPE mt        = pPart->getElement()->getMaterial();
						const P_TYPE pt        = pPart->getPType();
						if(_Nokmc[mt][pt].empty()) //particle not mapped.
							continue;
						_Nokmc[mt][pt](pNode[ix][iy][iz]._index) += it->second;
					}
				}
			}
		}
	}
}

void Poisson::buildLKMCVector() {
	MeshNode ***pNode = _pMesh->getNodes();
	const IO::ParameterManager *pPM       = Domains::global()->PM();

	for(M_TYPE mt = 0; mt < pPM->getNMaterials(); ++mt)
		if (!_Nlkmc[mt].empty())
		{
			_Nlkmc[mt] = ublas::zero_vector<double>(_nx * _ny * _nz);
			for (size_t ix = 0; ix < _nx; ++ix)
				for (size_t iy = 0; iy < _ny; ++iy)
					for (size_t iz = 0; iz < _nz; ++iz)
						if (pNode[ix][iy][iz]._active && !pNode[ix][iy][iz]._mLA.empty())
							for (std::map<LKMC::LatticeAtom *, double>::const_iterator it = pNode[ix][iy][iz]._mLA.begin(); it != pNode[ix][iy][iz]._mLA.end(); ++it)
								if (/*!it->first->getDefective() && */!it->first->getPerformed())
									_Nlkmc[mt](pNode[ix][iy][iz]._index) += it->second;
		}
}

double Poisson::getPermittivity(MeshNode *n1, MeshNode *n2, MeshNode *n3) {
	double eps = 0.;

	std::set<MeshElement *> s1;
	std::set<MeshElement *> s2;
	std::set<MeshElement *> s3;

	_pMesh->getElementsFromNode(n1, s1);
	_pMesh->getElementsFromNode(n2, s2);
	_pMesh->getElementsFromNode(n3, s3);

	for(std::set<MeshElement *>::const_iterator it = s1.begin(); it != s1.end(); ++it)
		if (s2.find(*it) != s2.end() && s3.find(*it) != s3.end()) {
			// eps = Domains::global()->PM()->getPermittivity((*it)->getMaterial());
			eps = Domains::global()->PM()->getMaterial((*it)->getMaterial())._permittivity;
			break ;
		}

	return eps;
}

void Poisson::buildLaplacianMatrix() {
	double dxm;
	double dxp;
	double dym;
	double dyp;
	double dzm;
	double dzp;

	MeshNode ***pNode = _pMesh->getNodes();

	_laplacian = ublas::compressed_matrix<double> (_nx * _ny * _nz, _nx * _ny * _nz);
	_Qeff      = ublas::zero_vector<double> (_nx * _ny * _nz);

	for (unsigned ix = 0; ix < _pMesh->getnx(); ++ix) {
		for (unsigned iy = 0; iy < _pMesh->getny(); ++iy) {
			for (unsigned iz = 0; iz < _pMesh->getnz(); ++iz) {

				if (pNode[ix][iy][iz]._active) {
					double a[6] = {};

					unsigned ixm = ix - 1;
					unsigned ixp = ix + 1;
					unsigned iym = iy - 1;
					unsigned iyp = iy + 1;
					unsigned izm = iz - 1;
					unsigned izp = iz + 1;

					// xm
					if (ix == 0) {
						ixm = (_pMesh->getPeriodicX() ? _nx-1 : ix);
						dxm = (_pMesh->getPeriodicX() ? pNode[_nx][iy][iz]._coord._x - pNode[_nx-1][iy][iz]._coord._x : 0);
					} else
						dxm = pNode[ix][iy][iz]._coord._x   - pNode[ixm][iy][iz]._coord._x;
					// xp
					if (ix == _nx - 1) {
						ixp = (_pMesh->getPeriodicX() ? 0 : _nx-1);
						dxp = (_pMesh->getPeriodicX() ? pNode[_nx][iy][iz]._coord._x - pNode[_nx-1][iy][iz]._coord._x : 0);
					}
					else
						dxp = pNode[ixp][iy][iz]._coord._x - pNode[ix][iy][iz]._coord._x;
					// ym
					if (iy == 0) {
						iym = (_pMesh->getPeriodicY() ? _ny-1 : iy);
						dym = (_pMesh->getPeriodicY() ? pNode[ix][_ny][iz]._coord._y - pNode[ix][_ny-1][iz]._coord._y : 0);
					} else
						dym = pNode[ix][iy][iz]._coord._y   - pNode[ix][iym][iz]._coord._y;
					// yp
					if (iy == _ny - 1) {
						iyp = (_pMesh->getPeriodicY() ? 0 : _ny-1);
						dyp = (_pMesh->getPeriodicY() ? pNode[ix][_ny][iz]._coord._y - pNode[ix][_ny-1][iz]._coord._y : 0);
					}
					else
						dyp = pNode[ix][iyp][iz]._coord._y - pNode[ix][iy][iz]._coord._y;
					// zm
					if (iz == 0) {
						izm = (_pMesh->getPeriodicZ() ? _nz-1 : iz);
						dzm = (_pMesh->getPeriodicZ() ? pNode[ix][iy][_nz]._coord._z - pNode[ix][iy][_nz-1]._coord._z : 0);
					} else
						dzm = pNode[ix][iy][iz]._coord._z   - pNode[ix][iy][izm]._coord._z;
					// zp
					if (iz == _nz - 1) {
						izp = (_pMesh->getPeriodicZ() ? 0 : _nz-1);
						dzp = (_pMesh->getPeriodicZ() ? pNode[ix][iy][_nz]._coord._z - pNode[ix][iy][_nz-1]._coord._z : 0);
					}
					else
						dzp = pNode[ix][iy][izp]._coord._z - pNode[ix][iy][iz]._coord._z;


					dxm /= BOHR_RADIUS * 1e9;
					dxp /= BOHR_RADIUS * 1e9;
					dym /= BOHR_RADIUS * 1e9;
					dyp /= BOHR_RADIUS * 1e9;
					dzm /= BOHR_RADIUS * 1e9;
					dzp /= BOHR_RADIUS * 1e9;

					if (ix != ixp) {
						const double eps1 = getPermittivity(&pNode[ixp][iy][iz], &pNode[ix][iyp][iz], &pNode[ix][iy][izp]);
						const double eps2 = getPermittivity(&pNode[ixp][iy][iz], &pNode[ix][iym][iz], &pNode[ix][iy][izp]);
						const double eps3 = getPermittivity(&pNode[ixp][iy][iz], &pNode[ix][iyp][iz], &pNode[ix][iy][izm]);
						const double eps4 = getPermittivity(&pNode[ixp][iy][iz], &pNode[ix][iym][iz], &pNode[ix][iy][izm]);

						a[0] = 1. / (4. * dxp) * (eps1 * dyp * dzp + eps2 * dym * dzp + eps3 * dyp * dzm + eps4 * dym * dzm);
					}
					if (ix != ixm) {
						const double eps1 = getPermittivity(&pNode[ixm][iy][iz], &pNode[ix][iyp][iz], &pNode[ix][iy][izp]);
						const double eps2 = getPermittivity(&pNode[ixm][iy][iz], &pNode[ix][iym][iz], &pNode[ix][iy][izp]);
						const double eps3 = getPermittivity(&pNode[ixm][iy][iz], &pNode[ix][iyp][iz], &pNode[ix][iy][izm]);
						const double eps4 = getPermittivity(&pNode[ixm][iy][iz], &pNode[ix][iym][iz], &pNode[ix][iy][izm]);

						a[1] = 1. / (4. * dxm) * (eps1 * dyp * dzp + eps2 * dym * dzp + eps3 * dyp * dzm + eps4 * dym * dzm);
					}
					if (iy != iyp) {
						const double eps1 = getPermittivity(&pNode[ix][iyp][iz], &pNode[ixp][iy][iz], &pNode[ix][iy][izp]);
						const double eps2 = getPermittivity(&pNode[ix][iyp][iz], &pNode[ixm][iy][iz], &pNode[ix][iy][izp]);
						const double eps3 = getPermittivity(&pNode[ix][iyp][iz], &pNode[ixp][iy][iz], &pNode[ix][iy][izm]);
						const double eps4 = getPermittivity(&pNode[ix][iyp][iz], &pNode[ixm][iy][iz], &pNode[ix][iy][izm]);

						a[2] = 1. / (4. * dyp) * (eps1 * dxp * dzp + eps2 * dxm * dzp + eps3 * dxp * dzm + eps4 * dxm * dzm);
					}
					if (iy != iym) {
						const double eps1 = getPermittivity(&pNode[ix][iym][iz], &pNode[ixp][iy][iz], &pNode[ix][iy][izp]);
						const double eps2 = getPermittivity(&pNode[ix][iym][iz], &pNode[ixm][iy][iz], &pNode[ix][iy][izp]);
						const double eps3 = getPermittivity(&pNode[ix][iym][iz], &pNode[ixp][iy][iz], &pNode[ix][iy][izm]);
						const double eps4 = getPermittivity(&pNode[ix][iym][iz], &pNode[ixm][iy][iz], &pNode[ix][iy][izm]);

						a[3] = 1. / (4. * dym) * (eps1 * dxp * dzp + eps2 * dxm * dzp + eps3 * dxp * dzm + eps4 * dxm * dzm);
					}
					if (iz != izp) {
						const double eps1 = getPermittivity(&pNode[ix][iy][izp], &pNode[ixp][iy][iz], &pNode[ix][iyp][iz]);
						const double eps2 = getPermittivity(&pNode[ix][iy][izp], &pNode[ixp][iy][iz], &pNode[ix][iym][iz]);
						const double eps3 = getPermittivity(&pNode[ix][iy][izp], &pNode[ixm][iy][iz], &pNode[ix][iyp][iz]);
						const double eps4 = getPermittivity(&pNode[ix][iy][izp], &pNode[ixm][iy][iz], &pNode[ix][iym][iz]);

						a[4] = 1. / (4. * dzp) * (eps1 * dyp * dxp + eps2 * dym * dxp + eps3 * dyp * dxm + eps4 * dym * dxm);
					}
					if (iz != izm) {
						const double eps1 = getPermittivity(&pNode[ix][iy][izm], &pNode[ixp][iy][iz], &pNode[ix][iyp][iz]);
						const double eps2 = getPermittivity(&pNode[ix][iy][izm], &pNode[ixp][iy][iz], &pNode[ix][iym][iz]);
						const double eps3 = getPermittivity(&pNode[ix][iy][izm], &pNode[ixm][iy][iz], &pNode[ix][iyp][iz]);
						const double eps4 = getPermittivity(&pNode[ix][iy][izm], &pNode[ixm][iy][iz], &pNode[ix][iym][iz]);

						a[5] = 1. / (4. * dzm) * (eps1 * dyp * dxp + eps2 * dym * dxp + eps3 * dyp * dxm + eps4 * dym * dxm);
					}

					if (!pNode[ixp][iy][iz]._active && _pMesh->isDirichletX())
						_Qeff(pNode[ix][iy][iz]._index) += a[0] * ELECTRONVOLT_TO_HARTREE(pNode[ixp][iy][iz]._V) / (4. * M_PI);
					if (!pNode[ixm][iy][iz]._active && _pMesh->isDirichletX())
						_Qeff(pNode[ix][iy][iz]._index) += a[1] * ELECTRONVOLT_TO_HARTREE(pNode[ixm][iy][iz]._V) / (4. * M_PI);

					if (!pNode[ix][iyp][iz]._active && _pMesh->isDirichletY())
						_Qeff(pNode[ix][iy][iz]._index) += a[2] * ELECTRONVOLT_TO_HARTREE(pNode[ix][iyp][iz]._V) / (4. * M_PI);
					if (!pNode[ix][iym][iz]._active && _pMesh->isDirichletY())
						_Qeff(pNode[ix][iy][iz]._index) += a[3] * ELECTRONVOLT_TO_HARTREE(pNode[ix][iym][iz]._V) / (4. * M_PI);

					if (!pNode[ix][iy][izm]._active && _pMesh->isDirichletZ())
						_Qeff(pNode[ix][iy][iz]._index) += a[4] * ELECTRONVOLT_TO_HARTREE(pNode[ix][iy][izp]._V) / (4. * M_PI);
					if (!pNode[ix][iy][izp]._active && _pMesh->isDirichletZ())
						_Qeff(pNode[ix][iy][iz]._index) += a[5] * ELECTRONVOLT_TO_HARTREE(pNode[ix][iy][izm]._V) / (4. * M_PI);

					_laplacian(pNode[ix][iy][iz]._index, pNode[ix][iy][iz]._index)      = a[0] + a[1] + a[2] + a[3] + a[4] + a[5];
					if (ix != ixp)
						_laplacian(pNode[ix][iy][iz]._index, pNode[ixp][iy][iz]._index) = - a[0];
					if (ix != ixm)
						_laplacian(pNode[ix][iy][iz]._index, pNode[ixm][iy][iz]._index) = - a[1];
					if (iy != iyp)
						_laplacian(pNode[ix][iy][iz]._index, pNode[ix][iyp][iz]._index) = - a[2];
					if (iy != iym)
						_laplacian(pNode[ix][iy][iz]._index, pNode[ix][iym][iz]._index) = - a[3];
					if (iz != izp)
						_laplacian(pNode[ix][iy][iz]._index, pNode[ix][iy][izp]._index) = - a[4];
					if (iz != izm)
						_laplacian(pNode[ix][iy][iz]._index, pNode[ix][iy][izm]._index) = - a[5];
				}
			}
		}
	}
}
