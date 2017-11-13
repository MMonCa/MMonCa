/*
 * Author: ignacio.martin@imdea.org
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

#include "SaveCmd.h"
#include "kernel/Domain.h"
#include "kernel/Mesh.h"
#include "kernel/Coordinates.h"
#include "domains/SimData.h"
#include "okmc/Particle.h"
#include "lkmc/LatticeSite.h"
#include "lkmc/LatticeParam.h"
#include "domains/MCClient.h"
#include "io/ParameterManager.h"
#include "io/VtkWriter.h"
#include <string>
#include <fstream>
#include <sstream>
#include <map>


#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

using std::string;
using std::endl;
using std::vector;
using std::pair;
using Kernel::Coordinates;

namespace IO {

namespace ublas = boost::numeric::ublas;

SaveCmd::SaveCmd(Tcl_Interp *p, int argc, const char *argv[]) : Command(p, argc, argv)
{
}

int SaveCmd::operator()()
{	
	bool bDefects = specified("lkmc.defect");
	vector<string> defects;
	if(specified("defects"))
		defects = getStrings("defects", 1);

	if(specified("xyz"))
	{
		string filename = getString("xyz");
		filename+=".xyz";
		save(filename, XYZ_TYPE, bDefects, defects);
	}
	if(specified("lammps")) {
		string filename = getString("lammps");
		filename+=".dump";
		save(filename, LAMMPS_TYPE, bDefects,defects);
	}
	if(specified("atomeye")) {
		string filename = getString("atomeye");
		filename+=".cfg";
		save(filename, ATOMEYE_TYPE, bDefects, defects);
	}
	if(specified("csv"))
	{
		string filename = getString("csv");
		filename+=".csv";
		save(filename, CSV_TYPE, bDefects, defects);
	}
    if(specified("vtk"))
    {
        string filename = getString("vtk");
        filename+=".vti";
        save(filename, VTK_TYPE, bDefects, defects);
    }
	return TCL_OK;
}

// Be careful, experimental feature!
void SaveCmd::writeAtomeye(std::ofstream &ofs, std::vector<Domains::ParticleData> &part) const {
	Coordinates m, M;
	Domains::global()->client()->getCell(m, M);

	// Note about scale: there is no reason to use the parameter issued by the user
	double scale = 10.; // nm -> A

	ublas::matrix<double> H    = ublas::zero_matrix<double>(3,3);
	ublas::matrix<double> Hinv = ublas::zero_matrix<double>(3,3);

	H(0,0) = (M._y - m._y) * scale;
	H(0,1) = 0.;
	H(0,2) = 0.;

	H(1,0) = 0.;
	H(1,1) = (M._x - m._x) * scale;
	H(1,2) = 0.;

	H(2,0) = 0.;
	H(2,1) = 0.;
	H(2,2) = (M._z - m._z) * scale;

	// Compute invert matrix
	ublas::permutation_matrix<std::size_t> pm(H.size1());

	int res = lu_factorize(H,pm);
	if( res != 0 ) {
		WARNINGMSG("SaveCmd::writeAtomeye: aborting save command because of a problem during LU factorization!");
		return ;
	}
	Hinv.assign(ublas::identity_matrix<double>(H.size1()));
	lu_substitute(H, pm, Hinv);

	ofs << "Number of particles = " << part.size() << std::endl;
	ofs << "A = 1. Angstrom (basic length-scale)" << std::endl;

	for (unsigned i = 0; i < 3; ++i) {
		for (unsigned j = 0; j < 3; ++j) {
		  ofs << "H0(" << i + 1 << "," << j + 1 << ") = " <<  H(i,j) << " A" << std::endl;
		}
	}

	ofs << "Transform(1,1) = 1" << std::endl;
	ofs << "Transform(1,2) = 0" << std::endl;
	ofs << "Transform(1,3) = 0" << std::endl;
	ofs << "Transform(2,1) = 0" << std::endl;
	ofs << "Transform(2,2) = 1" << std::endl;
	ofs << "Transform(2,3) = 0" << std::endl;
	ofs << "Transform(3,1) = 0" << std::endl;
	ofs << "Transform(3,2) = 0" << std::endl;
	ofs << "Transform(3,3) = 1" << std::endl;

	for(vector<Domains::ParticleData>::iterator it = part.begin(); it != part.end(); ++it) {
		ublas::vector<double> x = ublas::zero_vector<double>(3);
		ublas::vector<double> s;

		x(0) = it->get<1>()._y * scale;
		x(1) = it->get<1>()._x * scale;
		x(2) = it->get<1>()._z * scale;

		s = ublas::prod(Hinv, x);

		// Only Si export supported for the moment
		// need to store atomic mass (in a.m.u) in config to support other species
		ofs << "28.0855 Si " << s(0) << " " << s(1) << " " << s(2) <<  " 0 0 0 " << std::endl;
	}
}

void SaveCmd::save(const string &filename, FORMAT_TYPES ft, bool bDefects, const vector<string> &defects) const
{
	Coordinates m, M;
	Domains::global()->client()->getCell(m, M);

	vector<Domains::ParticleData> data;
    if(specified("lattice"))
		for(Domains::MeshElementIterator it= Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
		{
			LKMC::LatticeSite *part = it->getFirstLS();
			if(part)
			{
				while(part)
				{
					std::string name = Domains::global()->PM()->getParticleName(it->getMaterial(), part->getPType());
					data.push_back(Domains::ParticleData(it->getMaterial(), part->getCoordinates(), part->getPType(), Kernel::Event::LATTICEATOM, name));
					part = part->getNext();
				}
			}
			else // OKMC Case (build the whole lattice)
			{
				Kernel::M_TYPE mt = it->getMaterial();
				if(Domains::global()->PM()->getMaterialName(mt) == "Gas")
					continue;
				Coordinates m, M;
				it->getCorners(m, M);
				vector<LKMC::Lattice::LatticeInformation> atoms;
				it->getDomain()->_pLat[mt]->fill(m, M, atoms, false);
				float x = it->getAlloyFraction();
				for(vector<LKMC::Lattice::LatticeInformation>::iterator it = atoms.begin(); it!=atoms.end();++it)
				{
					if(Domains::global()->client()->rand() < x)
					{
						std::string name = Domains::global()->PM()->getAlloyName(mt);
						data.push_back(Domains::ParticleData(mt, it->_coords, 1, Kernel::Event::LATTICEATOM, name));
					}
					else
					{
						std::string name = Domains::global()->PM()->getParticleName(mt, it->_type);
						data.push_back(Domains::ParticleData(mt, it->_coords, 0, Kernel::Event::LATTICEATOM, name));
					}
				}
			}
		}
	else
	{
		Domains::global()->data()->getLKMCInterface(data, bDefects);
    	Domains::global()->data()->getOKMCParticles(data, defects);
	}

	std::stringstream ff;
	ff << filename;
	if(!specified("no.print"))
		LOWMSG2("Writing " << data.size() << " atoms in "<< ff.str());
	std::ios_base::openmode mode = std::ios_base::out;
	if(specified("append"))
		mode |= std::ios_base::app;
	std::ofstream ofs(ff.str().c_str(), mode);

	float scale = specified("scale") ? getFloat("scale") : 1.;

	if (ft == LAMMPS_TYPE) {
		m *= scale;
		M *= scale;
		ofs << "ITEM: TIMESTEP\n";
		ofs << Domains::global()->getTime() << std::endl;
		ofs << "ITEM: NUMBER OF ATOMS\n";
		ofs << data.size() << std::endl;
		ofs << "ITEM: BOX BOUNDS\n";
		ofs << m._x << "\t" << M._x << std::endl;
		ofs << m._y << "\t" << M._y << std::endl;
		ofs << m._z << "\t" << M._z << std::endl;
		ofs << "ITEM: ATOMS id type x y z defect iddef\n";
    }
    else if (ft == ATOMEYE_TYPE) {
		writeAtomeye(ofs, data);
	}
    else if (ft == XYZ_TYPE) {
		if(data.size())
		{
			ofs << data.size() << endl;
			ofs << "type x y z def ID defectID" << endl;
		}
		else
			WARNINGMSG("Not writing 0 atoms in " << ff.str());
	}
    else if (ft == VTK_TYPE) {
        VtkWriter::VtkFormat format = VtkWriter::BINARY;
        if (specified("ascii")) format = VtkWriter::ASCII;
        
        Kernel::Mesh *pme = Domains::global()->getDomain(0)->_pMesh;

        const std::vector<float> &xl = pme->getLines(0);
        const std::vector<float> &yl = pme->getLines(1);
        const std::vector<float> &zl = pme->getLines(2);

        double dx = xl[1]-xl[0];
        double dy = yl[1]-yl[0];
        double dz = zl[1]-zl[0];

        int nx = xl.size()-1;
        int ny = yl.size()-1;
        int nz = zl.size()-1;

        std::vector<int> mat;
        std::vector<double> pot;
        std::vector<double> strainxx, strainyy, strainzz, strainxz, strainxy, strainyz;
        std::vector<double> stressxx, stressyy, stresszz, stressxz, stressxy, stressyz;

        for (int k=0; k<nz; ++k) {
            for (int j=0; j<ny; ++j) {
                for (int i=0; i<nx; ++i) {
                    int idx=pme->getIndexFromIndices(i,j,k);
                    mat.push_back((int)pme->getElement(idx)->getMaterial());
                    pot.push_back(pme->getElement(idx)->electrostaticPotential());
                    strainxx.push_back(pme->getElement(idx)->strain_xx());
                    strainyy.push_back(pme->getElement(idx)->strain_yy());
                    strainzz.push_back(pme->getElement(idx)->strain_zz());
                    strainxy.push_back(pme->getElement(idx)->strain_xy());
                    strainxz.push_back(pme->getElement(idx)->strain_xz());
                    strainyz.push_back(pme->getElement(idx)->strain_yz());
                    stressxx.push_back(pme->getElement(idx)->stress_xx());
                    stressyy.push_back(pme->getElement(idx)->stress_yy());
                    stresszz.push_back(pme->getElement(idx)->stress_zz());
                    stressxy.push_back(pme->getElement(idx)->stress_xy());
                    stressxz.push_back(pme->getElement(idx)->stress_xz());
                    stressyz.push_back(pme->getElement(idx)->stress_yz());
                }
            }
        }

        VtkWriterImageData vti;
        
        vti.setOrigin(m._x,m._y,m._z);
        vti.setSpacing(dx,dy,dz);
        vti.setSize(nx, ny, nz);

        vti.addCellData("Materials", mat, format);
        vti.addCellData("ElectrostaticPotential", pot, format);
        vti.addCellData("StrainXX", strainxx, format);
        vti.addCellData("StrainYY", strainyy, format);
        vti.addCellData("StrainZZ", strainzz, format);
        vti.addCellData("StrainXY", strainxy, format);
        vti.addCellData("StrainXZ", strainxz, format);
        vti.addCellData("StrainYZ", strainyz, format);
        vti.addCellData("StressXX", stressxx, format);
        vti.addCellData("StressYY", stressyy, format);
        vti.addCellData("StressZZ", stresszz, format);
        vti.addCellData("StressXY", stressxy, format);
        vti.addCellData("StressXZ", stressxz, format);
        vti.addCellData("StressYZ", stressyz, format);
        vti.write(ofs);
    }
    else {
		ofs << " x, y, z, particle, defect, name\n";
	}

	const unsigned nParticles = Domains::global()->PM()->getNParticles();
	for(vector<Domains::ParticleData>::iterator it=data.begin(); it!=data.end(); ++it) {
		switch (ft)
		{
		case XYZ_TYPE:
			ofs << Domains::global()->PM()->getParticleName(it->get<0>(), it->get<2>() % nParticles) << " " <<
			it->get<1>()._x*scale << " " << it->get<1>()._y*scale << " " << it->get<1>()._z*scale << " " <<
			it->get<3>() << " " << it->get<4>() << " " << it->get<5>() << endl;
			break;
		case LAMMPS_TYPE:
			ofs << it-data.begin() +1 << "\t" << it->get<2>() << "\t" <<
			it->get<1>()._x*scale << "\t" << it->get<1>()._y*scale << "\t" << it->get<1>()._z*scale << "\t" <<
			it->get<3>() << "\t" << it->get<5>() << std::endl;
			break;
		case ATOMEYE_TYPE:
			break;
        case VTK_TYPE:
            break;
		case CSV_TYPE:
			ofs << it->get<1>()._x*scale << ", " << it->get<1>()._y*scale << ", " << it->get<1>()._z*scale << ", " <<
			it->get<2>() << ", " <<
			it->get<3>() << ", " <<
			it->get<4>() << ", " <<
			it->get<5>() << std::endl;
			break;
		}
	}
	ofs.close();
	if(!specified("no.print"))
		LOWMSG(" done!.");
}

}
