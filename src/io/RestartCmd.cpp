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
 *
 * RestartCmd.cpp
 *
 *  Created on: Sep 25, 2014
 */

#include "RestartCmd.h"

using std::string;

#include "version.h"
#include "domains/MCClient.h"
#include "io/FileParameters.h"
#include "kernel/Mesh.h"
#include "kernel/MeshElement.h"
#include "kernel/RateManager.h"

#include "okmc/MobileParticle.h"
#include "okmc/Cluster.h"
#include "okmc/Interface.h"

#include "lkmc/LatticeAtomBCC1.h"
#include "lkmc/LatticeAtomDEpi.h"
#include "lkmc/LatticeAtomDSPER.h"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>
#include <sstream>

namespace IO {

RestartCmd::RestartCmd(Tcl_Interp *p, int argc, const char *argv[]) : Command(p, argc, argv)
{


}

int RestartCmd::operator()()
{
	std::string versionString = "MMonCa";
	versionString += VERSION;
	if(specified("save"))
	{
		string filename = getString("save");
		filename += ".mmonca";
		LOWMSG2("Writing restart information in " << filename);
		std::ios_base::openmode mode = std::ofstream::binary;
		std::ofstream file(filename.c_str(), mode);
		if(file.fail())
			ERRORMSG("Cannot open " << filename << " for writing!");
		boost::iostreams::filtering_streambuf< boost::iostreams::input> in;
		in.push( boost::iostreams::gzip_compressor());

		std::stringstream os;
		os << versionString << std::endl;
		restart(*static_cast<std::ostream *>(&os));
		in.push(os);

		boost::iostreams::copy(in, file);
		LOWMSG(". Done");
	}

	if(specified("load"))
	{
		string filename = getString("load");
		filename += ".mmonca";
		LOWMSG("Reading restart information from " << filename << "\n{\n");

		std::ifstream file(filename.c_str(), std::ofstream::binary);
		if(file.fail())
			ERRORMSG("Cannot open " << filename << " for reading!");
		boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
		in.push(boost::iostreams::gzip_decompressor());
		in.push(file);
		std::stringstream is;
		boost::iostreams::copy(in, is);

		std::string vers;
		is >> vers;
		if(vers != versionString)
			WARNINGMSG("MMonCa restart file version " << vers <<
					" might not be compatible with binary version " << versionString);
		restart(*static_cast<std::istream *>(&is));
		LOWMSG("} . Restarting done");
	}
	return TCL_OK;
}

void RestartCmd::restart(std::istream &is) const
{
	{
		string check;
		is >> check;
		if(check != "Parameters")
			ERRORMSG("Restarting MMonCa file is corrupted. " << "Parameters" << " not found");
		Domains::global()->setFileParameters(new FileParameters(is));
	}
	// restart client and set temperature
	{
		string check;
		is >> check;
		if(check != "MCClient")
			ERRORMSG("Restarting MMonCa file is corrupted. " << "MCClient" << " not found");
		Domains::MCClient *pCli = new Domains::MCClient(is);
		Kernel::Coordinates m, M;
		pCli->getCell(m, M);
		Domains::global()->setClient(pCli, m, M);
		float temp;
		is >> temp;
		Domains::global()->setTempK(temp);
	}
	//update MeshElements
	{
		string check;
		is >> check;
		if(check != "MeshElements")
			ERRORMSG("Restarting MMonCa file is corrupted. " << "MeshElement" << " not found");
		for(Domains::MeshElementIterator it = Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
			it.modify()->restart(is);
	}
	//Lattice atoms
	//the LKMCModel has detected the initialization using a istream, and has not created any LatticeAtom
	{
		string check;
		is >> check;
		if(check != "LatticeAtoms")
			ERRORMSG("Restarting MMonCa file is corrupted. " << "LatticeAtoms" << " not found");
		unsigned nLA;
		is >> nLA;
		for(unsigned i=0; i<nLA; ++i)
		{
			string atom_type;
			is >> atom_type;
			LKMC::LatticeAtom * pLA = 0;
			if(atom_type == "BCC1")
				pLA = new LKMC::LatticeAtomBCC1(is);
			else if(atom_type == "DSPE")
				pLA = new LKMC::LatticeAtomDSPER(is);
			else if(atom_type == "DEPI")
				pLA = new LKMC::LatticeAtomDEpi(is);
			else
				LOWMSG("Restarting MMonCa file is corrupted. Defect " << atom_type << " not expected.");
			pLA->getElement()->incLA(pLA->getPerformed());
			pLA->getDomain()->_pRM->insert(pLA, pLA->getElement());
		}
		//rebuild rates for atoms
		for(Domains::MeshElementIterator it = Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
		{
			LKMC::LatticeSite *pLS = it->getFirstLS();
			while(pLS)
			{
				LKMC::LatticeAtom *pLA = dynamic_cast<LKMC::LatticeAtom *>(pLS);
				if(pLA)
					pLA->getDomain()->_pRM->update(pLA, pLA->getElement());
				pLS = pLS->getNext();
			}
		}
	}
	//update events
	{
		string check;
		is >> check;
		if(check != "Defects")
			ERRORMSG("Restarting MMonCa file is corrupted. " << "Event" << " not found");
		unsigned nDefects;
		is >> nDefects;
		for(unsigned i=0; i<nDefects; ++i)
		{
			string defect_type;
			is >> defect_type;
			if(defect_type == "MP")
				new OKMC::MobileParticle(is);
			else if(defect_type == "CL")
			{
				OKMC::Cluster * pCl = new OKMC::Cluster(is);
				pCl->getDomain()->_pRM->insert(pCl, pCl->getElement());
			}
			else if(defect_type == "IF")
				OKMC::Interface::restart(is);
			else
				LOWMSG("Restarting MMonCa file is corrupted. Defect " << defect_type << " not expected.");
		}
	}
	{
		//time and rate managers
		double time;
		is >> time;
		Domains::global()->setTime(time);
		string check;
		is >> check;
		if(check != "RateManager")
			ERRORMSG("Restarting MMonCa file is corrupted. " << "RateManager" << " not found");

		for(unsigned d=0; d < Domains::global()->getDomains(); ++d)
			Domains::global()->getDomain(d)->_pRM->restart(is);
	}
}

//saves the information in a format that allows restarting.
void RestartCmd::restart(std::ostream &os) const
{
	os << "Parameters ";
	Domains::global()->getFileParameters()->restart(os);
	os << "\nMCClient ";
	//save Client information
	Domains::global()->client()->restart(os);
	os << Domains::global()->getTempK();
	//save MeshElement information
	os << "\nMeshElements ";
	for(Domains::MeshElementIterator it = Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
		it->restart(os);
	os << "\nLatticeAtoms ";
	{
		unsigned nLA = 0;
		for(Domains::MeshElementIterator it = Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
		{
			LKMC::LatticeSite *pLS = it->getFirstLS();
			while(pLS)
			{
				if(dynamic_cast<LKMC::LatticeAtom *>(pLS))
					nLA++;
				pLS = pLS->getNext();
			}
		}
		os << nLA << "\n";
		for(Domains::MeshElementIterator it = Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
		{
			LKMC::LatticeSite *pLS = it->getFirstLS();
			while(pLS)
			{
				LKMC::LatticeAtom *pLA = dynamic_cast<LKMC::LatticeAtom *>(pLS);
				if(pLA)
				{
					os << pLA->getClass() << " ";
					pLA->restart(os);
					os << std::endl;
				}
				pLS = pLS->getNext();
			}
		}
	}
	os << "\nDefects ";
	unsigned nDef = 0;
	for(Domains::DefectIterator it = Domains::global()->beginDI(); it != Domains::global()->endDI(); ++it)
		nDef++;
	os << nDef << "\n";
	for(Domains::DefectIterator it = Domains::global()->beginDI(); it != Domains::global()->endDI(); ++it)
	{
		switch(it->getEType())
		{
		case Kernel::Event::MOBILEPARTICLE:
			os << "MP ";
			break;
		case Kernel::Event::CLUSTER:
			os << "CL ";
			break;
		case Kernel::Event::INTERFACE:
			os << "IF ";
			break;
		default:
			ERRORMSG("Event " << int(it->getEType()) << " not supported in restart.");
		}
		it->restart(os);
		os << std::endl;
	}
	{
		//RateManagers
		os << Domains::global()->getTime() << std::endl;
		os << "RateManager ";
		for(unsigned d=0; d < Domains::global()->getDomains(); ++d)
			Domains::global()->getDomain(d)->_pRM->restart(os);
		os << std::endl;
	}
}


} /* namespace IO */
