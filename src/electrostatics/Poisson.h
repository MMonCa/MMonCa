/*
 * Author: Benoit Sklenard benoit.sklenard@cea.fr 
 * 
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

#ifndef POISSON_H
#define POISSON_H

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace IO { class ParameterManager; class FileParameters; }
namespace Kernel { class Mesh; class Domain; }

struct Tcl_Interp;

namespace Electrostatics {

class MeshNode;

class Poisson {

	typedef struct {
		double                                _Nc;
		double                                _Nv;
		double                                _Eg;
		double                                _Eg0;
		boost::numeric::ublas::vector<double> _Dc;
		boost::numeric::ublas::vector<double> _Dv;

		boost::numeric::ublas::vector<double> _Dcx;
		boost::numeric::ublas::vector<double> _Dcy;
		boost::numeric::ublas::vector<double> _Dcz;

		boost::numeric::ublas::vector<double> _Dvb;
		boost::numeric::ublas::vector<double> _Dvd;
	} MATERIAL_PARAM;


        Kernel::Domain     *_pDomain;
        Kernel::Mesh       *_pMesh;
	Tcl_Interp *_pTcl;

	unsigned    _nx;
	unsigned    _ny;
	unsigned    _nz;

	double      _T;

	double      _CG_tol;
	unsigned    _CG_maxIter;
	unsigned    _CG_verbose;

	double      _NR_tol;
	unsigned    _NR_maxIter;
	unsigned    _NR_verbose;

	bool        _FermiDirac;
	bool        _partialIonization;
	bool        _selfConsistentLKMC;
	// bool        _selfConsistentGFLS;
	// bool        _selfConsistentGFLS2;

	boost::numeric::ublas::compressed_matrix<double> _laplacian;

	boost::numeric::ublas::vector<double>            _psi;
	boost::numeric::ublas::vector<double>            _Qe;      // electrons charge
	boost::numeric::ublas::vector<double>            _Qh;      // holes charge
	boost::numeric::ublas::vector<double>            _Qokmc;   // fixed impurities charge

	boost::numeric::ublas::vector<double>            _Qeff;    // effective charge (resulting from Dirichlet BC)

	boost::numeric::ublas::vector<double>            _Qlkmc_M;
	boost::numeric::ublas::vector<double>            _Qlkmc_P;

	boost::numeric::ublas::vector<double>            _DEc;
	boost::numeric::ublas::vector<double>            _DEv;

	boost::numeric::ublas::vector<double>            _Eg;
	boost::numeric::ublas::vector<double>            _Eg0;
	boost::numeric::ublas::vector<double>            _Nc;
	boost::numeric::ublas::vector<double>            _Nv;

	boost::numeric::ublas::vector<double>            _Nokmc[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES];
	boost::numeric::ublas::vector<double>            _Nlkmc[Kernel::MAX_MATERIALS];

	MATERIAL_PARAM                                   _matParam[Kernel::MAX_MATERIALS];

	// void printPotential3D(std::string &);
	// void printPotential2D(std::string &);
	// void printPotential1D(std::string &);
	// void map2Elements();


	void map2Nodes();
	void map2Elements();

	void buildLaplacianMatrix();
	void buildOKMCVector();
	void buildLKMCVector();

	// void updateOKMCVector();
	// void buildCharge();

	bool NewtonRaphson();

	bool ConjugateGradient(const boost::numeric::ublas::compressed_matrix<double> &, boost::numeric::ublas::vector<double> &, const boost::numeric::ublas::vector<double> &);

	void test(unsigned, unsigned, unsigned);
	void testPrint(unsigned, unsigned, unsigned);
	void testPN();

	double getPermittivity(MeshNode *, MeshNode *, MeshNode *);
	void   computeElectronicParameters();
	void   setFirstGuess();

	// double getCharge(OKMC::Particle *); // OBSOLETE
	double getLA(MeshNode *); // OBSOLETE
	void updateLA();

	double getNc(MeshNode *);
	double getNv(MeshNode *);
	double getEg(MeshNode *);
	double getEg0(MeshNode *);
	double getDEc(MeshNode *);
	double getDEv(MeshNode *);

public:
	Poisson(Tcl_Interp *, Kernel::Domain *);
	~Poisson();

	bool compute(const double);
};

}
#endif /* ! POISSON_H */
