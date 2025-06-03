/*
 * SimData.cpp
 *
 *  Created on: Feb 15, 2011
 *      Author: ignacio.martin@imdea.org
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

#include "SimData.h"
#include "MCClient.h"
#include "kernel/Mesh.h"
#include "kernel/Domain.h"
#include "kernel/RateManager.h"
#include "kernel/Constants.h"
#include "lkmc/LatticeParam.h"
#include "lkmc/LatticeAtom.h"
#include "lkmc/LatticeAtomDiamond.h"
#include "lkmc/Lattice.h"
#include "okmc/Particle.h"
#include "okmc/Defect.h"
#include "okmc/EDType.h"
#include "okmc/Interface.h"
#include "okmc/MobileParticle.h"
#include "okmc/Cluster.h"
#include "io/ParameterManager.h"
#include "kernel/StateManager.h"
#include "io/Diagnostic.h"
#include "io/Polynomial.h"
#include "io/FileParameters.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iomanip>

using std::vector;
using std::pair;
using std::map;
using std::string;
using std::setw;
using Kernel::Coordinates;
using Kernel::Mesh;

using namespace boost::numeric;

namespace Domains {

SimData::SimData(IO::FileParameters *pFP) : _lastTime(0), _dim(1)
{
	_odelta[0] = pFP->getFloat("MC/Mesh/delta.x");
	_odelta[1] = pFP->getFloat("MC/Mesh/delta.y");
	_odelta[2] = pFP->getFloat("MC/Mesh/delta.z");
}

void SimData::getLKMCInterface(vector<ParticleData>&c, bool bDefective) const
{
	const bool interface = true;
	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		LKMC::LatticeSite *first = it->getFirstLS();
		while(first)
		{
			LKMC::LatticeAtom *pLA = dynamic_cast<LKMC::LatticeAtom *>(first);
			if(pLA)
			{
				std::string name = Domains::global()->PM()->getElementName(pLA->getPType());
				if(interface)
				{
					// a/c interface = crystalline atoms having an amorphous NN
					if(pLA->isDefective())
						c.push_back(ParticleData(it->getMaterial(), pLA->getCoordinates(), pLA->getPType(), pLA->getEType(), name));
					else
						for(unsigned j=0; j < pLA->first(); ++j)
							if(pLA->getNeighbor(j) && pLA->getPerformed() && pLA->getNeighbor(j)->getPerformed() != pLA->getPerformed())
							{
								c.push_back(ParticleData(it->getMaterial(), pLA->getCoordinates(), pLA->getPType(), pLA->getEType(), name));
								break ;
							}
				}
				else
					c.push_back(ParticleData(it->getMaterial(), pLA->getCoordinates(), pLA->getPType(), pLA->getEType(), name));
			}
			first = first->getNext();
		}
	}
}

void SimData::getACInterfaceMean(Kernel::Coordinates &min, Kernel::Coordinates &max, Kernel::Coordinates &meanO) const
{
	std::set<LKMC::LatticeAtom *> ac;
	LKMC::LatticeAtom *pLA;
	Kernel::CoordinatesT<double> mean;

	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
	{
		pLA = static_cast<LKMC::LatticeAtom * >(it->getFirstLS());

		while (pLA) {
			for(unsigned j=0; j < pLA->first(); ++j)
				if(pLA->getNeighbor(j) /*&& amorph.find(pLA->getNeighbor(j)) == amorph.end()*/
						&& !pLA->getNeighbor(j)->getPerformed() && /*!pLA->getNeighbor(j)->getDefective() &&*/ pLA->getPerformed()) {
					double x = pLA->getCoordinates()._x;
					double y = pLA->getCoordinates()._y;
					double z = pLA->getCoordinates()._z;

					if (x >= min._x && x <= max._x && y >= min._y && y <= max._y && z >= min._z && z <= max._z) {
						ac.insert(pLA);
						break ;
					}
				}

			pLA = static_cast<LKMC::LatticeAtom * >(pLA->getNext());
		}
	}

	mean._x = 0.;
	mean._y = 0.;
	mean._z = 0.;

	for (std::set<LKMC::LatticeAtom *>::iterator it = ac.begin(); it != ac.end(); ++it) {
		mean._x += static_cast<double>((*it)->getCoordinates()._x);
		mean._y += static_cast<double>((*it)->getCoordinates()._y);
		mean._z += static_cast<double>((*it)->getCoordinates()._z);
	}
	meanO._x = ac.size() ? mean._x / ac.size() : 0.;
	meanO._y = ac.size() ? mean._y / ac.size() : 0.;
	meanO._z = ac.size() ? mean._z / ac.size() : 0.;
}

float SimData::getACInterfaceCoverage(Kernel::Coordinates &min, Kernel::Coordinates &max) const
{
	std::set<LKMC::LatticeAtom *> ac;
	LKMC::LatticeAtom *pLA;
	Kernel::CoordinatesT<double> mean;

	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
	{
		pLA = static_cast<LKMC::LatticeAtom * >(it->getFirstLS());

		while (pLA) {
			for(unsigned j=0; j < pLA->first(); ++j)
				if(pLA->getNeighbor(j) && pLA->getNeighbor(j)->getPerformed() && !pLA->getPerformed()) {
					double x = pLA->getCoordinates()._x;
					double y = pLA->getCoordinates()._y;
					double z = pLA->getCoordinates()._z;

					if (x >= min._x && x <= max._x && y >= min._y && y <= max._y && z >= min._z && z <= max._z) {
						ac.insert(pLA);
						break ;
					}
				}
			pLA = static_cast<LKMC::LatticeAtom * >(pLA->getNext());
		}
	}

	unsigned precursor = 0;

	for (std::set<LKMC::LatticeAtom *>::iterator it = ac.begin(); it != ac.end(); ++it)
		if((*it)->getState() == LKMC::LatticeAtom::LS_PRECURSOR)
			precursor++;
	return float(precursor)/float(ac.size());
}

void SimData::getACInterfaceMinMax(Kernel::Coordinates &min, Kernel::Coordinates &max, Kernel::Coordinates &acm, Kernel::Coordinates &acM) const
{
	std::set<LKMC::LatticeAtom *> ac;
	LKMC::LatticeAtom *pLA;
	acm = max;
	acM = min;

	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
	{
		pLA = static_cast<LKMC::LatticeAtom * >(it->getFirstLS());

		while (pLA)
		{
			for(unsigned j=0; j < pLA->first(); ++j)
				if(pLA->getNeighbor(j) /*&& amorph.find(pLA->getNeighbor(j)) == amorph.end()*/
						&& !pLA->getNeighbor(j)->getPerformed() && /*!pLA->getNeighbor(j)->getDefective() &&*/ pLA->getPerformed())
				{
					double x = pLA->getCoordinates()._x;
					double y = pLA->getCoordinates()._y;
					double z = pLA->getCoordinates()._z;

					if (x >= min._x && x <= max._x && y >= min._y && y <= max._y && z >= min._z && z <= max._z)
					{
						if( x <= acm._x) acm._x = x;
						if( y <= acm._y) acm._y = y;
						if( z <= acm._z) acm._z = z;
						if( x >= acM._x) acM._x = x;
						if( y >= acM._y) acM._y = y;
						if( z >= acM._z) acM._z = z;
						break;
					}
				}
			pLA = static_cast<LKMC::LatticeAtom * >(pLA->getNext());
		}
	}
}

void SimData::getACInterfaceStdev(Kernel::Coordinates &min, Kernel::Coordinates &max, Kernel::Coordinates &stdevO) const
{
	std::set<LKMC::LatticeAtom *> ac;
	LKMC::LatticeAtom *pLA;
	Kernel::CoordinatesT<double> mean, stdev; //to increase precision a little.

	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
	{
		pLA = static_cast<LKMC::LatticeAtom * >(it->getFirstLS());

		while (pLA) {
		/*	std::set<const LKMC::LatticeAtom *> amorph;
			for(unsigned i=0; i<pLA->first(); ++i) {
				if (pLA->getNeighbor(i) && !pLA->getNeighbor(i)->getPerformed()) {
					// amorphous
					unsigned def = 0;
					for(unsigned j=0; j<pLA->getNeighbor(i)->first(); ++j) {
						if (pLA->getNeighbor(i)->getNeighbor(j) && pLA->getNeighbor(i)->getNeighbor(j)->getDefective())
							++def;
					}

					if (def >= 1) {
						amorph.insert(pLA->getNeighbor(i));
					}
				}
			}*/
			for(unsigned j=0; j < pLA->first(); ++j) {
				if(pLA->getNeighbor(j) /*&& amorph.find(pLA->getNeighbor(j)) == amorph.end() */&& !pLA->getNeighbor(j)->getPerformed()/* && !pLA->getNeighbor(j)->getDefective()*/ && pLA->getPerformed()) {
					double x = pLA->getCoordinates()._x;
					double y = pLA->getCoordinates()._y;
					double z = pLA->getCoordinates()._z;

					if (x >= min._x && x <= max._x && y >= min._y && y <= max._y && z >= min._z && z <= max._z) {
						ac.insert(pLA);
						break ;
					}
				}
			}
			pLA = static_cast<LKMC::LatticeAtom * >(pLA->getNext());
		}
	}

	mean._x = 0.;
	mean._y = 0.;
	mean._z = 0.;

	for (std::set<LKMC::LatticeAtom *>::iterator it = ac.begin(); it != ac.end(); ++it) {
		mean._x += static_cast<double>((*it)->getCoordinates()._x) / static_cast<double>(ac.size());
		mean._y += static_cast<double>((*it)->getCoordinates()._y) / static_cast<double>(ac.size());
		mean._z += static_cast<double>((*it)->getCoordinates()._z) / static_cast<double>(ac.size());
	}

	stdev._x = 0.;
	stdev._y = 0.;
	stdev._z = 0.;
	for (std::set<LKMC::LatticeAtom *>::iterator it = ac.begin(); it != ac.end(); ++it) {
		stdev._x += pow(static_cast<double>((*it)->getCoordinates()._x) - mean._x, 2);
		stdev._y += pow(static_cast<double>((*it)->getCoordinates()._y) - mean._y, 2);
		stdev._z += pow(static_cast<double>((*it)->getCoordinates()._z) - mean._z, 2);
	}

	stdevO._x = ac.size() ? sqrt(stdev._x/ac.size()) : 0.;
	stdevO._y = ac.size() ? sqrt(stdev._y/ac.size()) : 0.;
	stdevO._z = ac.size() ? sqrt(stdev._z/ac.size()) : 0.;
}

void SimData::getOKMCParticles(std::vector<Domains::ParticleData> &data, const vector<string> &defects) const
{
	const unsigned nParticles = Domains::global()->PM()->getNParticles();
	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		OKMC::Particle *pPart = it->getFirstPart();
		while(pPart)
		{
			Kernel::Event::E_TYPE     et = pPart->getDefect()->getEType();
			Kernel::P_TYPE              pt = pPart->getPType();
			const Kernel::Coordinates co = pPart->getCoordinates();
			Kernel::M_TYPE            mt = it->getMaterial();
			string name;
			unsigned defectType = 0;
			switch(et)
			{
			case Kernel::Event::CLUSTER:
			{
				OKMC::Cluster *mc = static_cast<OKMC::Cluster *>(pPart->getDefect());
				name = mc->getDomain()->_pClPar->getParams(it->getMaterial(), mc->getEDType())->_name;
				defectType = 2+mc->getEDType();
			}
			break;
			case Kernel::Event::INTERFACE:
				defectType = 1;
				name = Kernel::Event::getEName(et);
			break;

			case Kernel::Event::MOBILEPARTICLE:
			default:
				defectType = 0;
				name = Kernel::Event::getEName(et);
			break;
			}
			if(defects.empty() || std::find(defects.begin(), defects.end(), name) != defects.end())
				data.push_back(Domains::ParticleData(mt, co, pt+defectType*nParticles, et, name, pPart->getDefect()->getIndex(0)));
			pPart = pPart->getNext();
		}

		unsigned howMany = it->getBAtoms();
		Kernel::Domain *pDomain = it->getDomain();
		Kernel::M_TYPE mt = it->getMaterial();
		Kernel::P_TYPE alloyPt = Domains::global()->PM()->getMaterial(mt)._alloy;
		while(howMany--)
		{
			Kernel::Coordinates co;
			pDomain->_pMesh->getRandomPosition(it->getIndex(), co,
				pDomain->_rng_dom.rand(), pDomain->_rng_dom.rand(), pDomain->_rng_dom.rand());
			if(defects.empty() || std::find(defects.begin(), defects.end(), "Alloy") != defects.end())
				data.push_back(Domains::ParticleData(mt, co, alloyPt, Kernel::Event::MOBILEPARTICLE, "Alloy", 0));
		}
	}

}

//include alloy particles
void SimData::collectDefects(map<string, unsigned> maps[Kernel::MAX_MATERIALS]) const
{
	Domains::DefectIterator itEnd = Domains::global()->endDI();
	for(Domains::DefectIterator it = Domains::global()->beginDI(); it!= Domains::global()->endDI(); ++it)
	{
		Kernel::ID theMap = it->getID();
		assert(theMap._mt < Domains::global()->PM()->getNMaterials());
		string name = Domains::global()->PM()->getIDName(theMap);
		Kernel::M_TYPE mt = it->getElement()->getMaterial();
		std::stringstream ss;
		switch(it->getEType())
		{
		case Kernel::Event::CLUSTER:
			ss << it->getDomain()->_pClPar->getParams(mt, static_cast<const OKMC::Cluster *>(*it)->getEDType())->_name;
			break;
		case Kernel::Event::INTERFACE:
			continue; //to be treated differently
		default:
			ss << Kernel::Event::getEName(it->getEType());
			break;
		}
		ss << "/" << name;
		if(maps[mt].find(ss.str()) == maps[mt].end())
			maps[mt][ss.str()] = 1;
		else
			maps[mt][ss.str()]++;
	}
	//interfaces only
	for(Domains::DefectIterator it = Domains::global()->beginDI(); it!= Domains::global()->endDI(); ++it)
	{
		const OKMC::Interface *pInt = dynamic_cast<const OKMC::Interface *> (*it);
		if(pInt)
		{
			Kernel::ID theMap = pInt->getID();
			string name0 = Domains::global()->PM()->getMaterialName(pInt->getElement(0)->getMaterial());
			string name1 = Domains::global()->PM()->getMaterialName(pInt->getElement(1)->getMaterial());
			std::stringstream ss;
			Kernel::M_TYPE mt;
			if(name0 < name1)
			{
				mt = pInt->getElement(0)->getMaterial();
				ss << name0 << "_" << name1;
			}
			else
			{
				mt = pInt->getElement(1)->getMaterial();
				ss << name1 << "_" << name0;
			}
			for(map<Kernel::P_TYPE, unsigned>::iterator it2 = theMap._pt.begin(); it2 != theMap._pt.end(); ++it2)
			{
				std::stringstream sss;
				sss << ss.str() << "/" << Domains::global()->PM()->getParticleName(theMap._mt, it2->first);
				if(maps[mt].find(sss.str()) == maps[mt].end())
					maps[mt][sss.str()] = it2->second;
				else
					maps[mt][sss.str()] += it2->second;
			}
		}
	}
	//alloy particles
	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		unsigned howMany = it->getBAtoms();
		if(howMany)
		{
			Kernel::M_TYPE mt = it->getMaterial();
			Kernel::P_TYPE pt = Domains::global()->PM()->getMaterial(mt)._alloy;
			std::stringstream ss;
			ss << Kernel::Event::getEName(Kernel::Event::MOBILEPARTICLE) << "/" << Domains::global()->PM()->getParticleName(mt, pt);
			if(maps[mt].find(ss.str()) == maps[mt].end())
				maps[mt][ss.str()] = it->getBAtoms();
			else
				maps[mt][ss.str()] += it->getBAtoms();
		}
	}
}

/*
 * orig_mt: Silicon, Iron, etc...
 * pt: I_TYPE, Carbon, etc...
 * defe: <100>, Cluster, etc...
 * theMap: B2I3, I35, etc...
 * defects. false:particle, true:defects
 */
//include alloy particles
unsigned SimData::count(Kernel::M_TYPE orig_mt, const string &partName, const std::string &defe,
		const string &sID, unsigned minSize, bool defects, const std::string &st) const
{
	unsigned counter=0;
	Kernel::P_TYPE names[Kernel::MAX_MATERIALS];
	for(unsigned mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		names[mt] = Domains::global()->PM()->getParticleNumber(mt, partName);
	Kernel::ID IDs[Kernel::MAX_MATERIALS];
	for(unsigned mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		IDs[mt] = Domains::global()->PM()->getID(mt, sID);
	for(Domains::DefectIterator it = Domains::global()->beginDI(); it!= Domains::global()->endDI(); ++it)
	{
		Kernel::M_TYPE mt = it->getElement()->getMaterial();
		Kernel::P_TYPE partType = names[mt];
		vector<const OKMC::Particle *> parts = it->getParticles();
		if(parts.size() < minSize)
			continue;
		for(vector<const OKMC::Particle *>::iterator itPart = parts.begin(); itPart!=parts.end(); ++itPart)
		{
			const OKMC::Particle *pPart = *itPart;
			Kernel::P_TYPE pt = pPart->getPType();
			Kernel::Event::E_TYPE ev = pPart->getDefect()->getEType();
			Kernel::M_TYPE mt = pPart->getElement()->getMaterial();
			OKMC::Cluster *pMC = dynamic_cast<OKMC::Cluster *>(pPart->getDefect());
			OKMC::Interface *pIF = dynamic_cast<OKMC::Interface *>(pPart->getDefect());
			OKMC::MobileParticle *pMP = dynamic_cast<OKMC::MobileParticle *>(pPart->getDefect());
			if(Domains::global()->PM()->getMaterial(mt)._binary == false && (pMC || pIF))
				pt = Domains::global()->PM()->getParticleForCluster(mt, pt, Kernel::POS_0);
			bool isMt = orig_mt == mt;
			bool isPart = pt == partType;
			bool isGenericDefect  = Kernel::Event::getEName(ev) == defe;
			bool isExtendedDefect = pMC && it->getDomain()->_pClPar->getParams(mt, pMC->getEDType())->_name == defe;
			bool isID = IDs[mt] == pPart->getDefect()->getID();
			bool isGenericState = st.empty();
			bool isState = (pMP && st == Domains::global()->PM()->getStateName(mt, pMP->getPType(), pMP->getState()));

			if( (orig_mt == Kernel::MAX_MATERIALS || isMt)   &&
					(partName.empty() || isPart) &&
					(defe.empty()     || isGenericDefect || isExtendedDefect) &&
					(sID.empty()      || isID) &&
					(isState          || isGenericState))
			{
				counter++;
				if(defects)
					break;
			}
		}
	}
	//alloy particles
	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		Kernel::M_TYPE mt = it->getMaterial();
		Kernel::P_TYPE partType = names[mt];
		unsigned howMany = it->getBAtoms();
		if(howMany)
		{
			Kernel::M_TYPE mt = it->getMaterial();
			Kernel::P_TYPE pt = Domains::global()->PM()->getMaterial(mt)._alloy;
			if(Domains::global()->PM()->getMaterial(mt)._binary == false)
				pt = Domains::global()->PM()->getParticle(mt, pt, Kernel::POS_0);
			Kernel::Event::E_TYPE ev = Kernel::Event::MOBILEPARTICLE;
			bool isMt = orig_mt == mt;
			bool isPart = pt == partType;
			bool isGenericDefect  = Kernel::Event::getEName(ev) == defe;

			if( (orig_mt == Kernel::MAX_MATERIALS || isMt)   &&
					(partName.empty() || isPart) &&
					(defe.empty()     || isGenericDefect))
			{
				counter+= it->getBAtoms();
				if(defects)
					break;
			}
		}
	}
	return counter;
}

//count number of positions in the simulation
unsigned SimData::countPos(Kernel::M_TYPE orig_mt, const std::string &spos,
		const std::string &defe, const std::string &sID) const
{
	unsigned counter=0;
	Kernel::P_POS positions[Kernel::MAX_MATERIALS];
	for(unsigned mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		positions[mt] = Domains::global()->PM()->getPositionNumber(mt, spos);
	Kernel::ID IDs[Kernel::MAX_MATERIALS];
	for(unsigned mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
			IDs[mt] = Domains::global()->PM()->getID(mt, sID);

	for(Domains::DefectIterator it = Domains::global()->beginDI(); it!= Domains::global()->endDI(); ++it)
	{
		Kernel::M_TYPE mt = it->getElement()->getMaterial();
		Kernel::ID theID = it->getID();
		bool isMt = orig_mt == mt;
		Kernel::Event::E_TYPE ev = (*it)->getEType();
		bool isGenericDefect  = Kernel::Event::getEName(ev) == defe;
		const OKMC::Cluster *pMC = dynamic_cast<const OKMC::Cluster *>(*it);
		bool isExtendedDefect = pMC && it->getDomain()->_pClPar->getParams(mt, pMC->getEDType())->_name == defe;
		bool isID = IDs[mt] == (*it)->getID();
		const OKMC::MobileParticle *pMP = dynamic_cast<const OKMC::MobileParticle *>(*it);
		
		if(pMP)
		{
			if( (orig_mt == Kernel::MAX_MATERIALS || isMt) &&
		        	 Domains::global()->PM()->getPPos(pMP->getPType()) == positions[mt] &&
		        	 (defe.empty()     || isGenericDefect || isExtendedDefect) &&
		        	 (sID.empty()      || isID)
		        	 )
		        	 	counter++;
		}
		else
			for(map<Kernel::P_POS, unsigned>::iterator posIt=theID._pos.begin(); posIt!=theID._pos.end(); ++posIt)
				if( (orig_mt == Kernel::MAX_MATERIALS || isMt) &&
		        	 posIt->first == positions[mt] &&
		        	 (defe.empty()     || isGenericDefect || isExtendedDefect) &&
		        	 (sID.empty()      || isID)
		        	 )
		        	 	counter += posIt->second;
	}
	return counter;
}

//computes the average radius of the specified defects
double SimData::radius(Kernel::M_TYPE orig_mt, const string &defect, const string &sID, double minRad) const
{
	unsigned count = 0;
	double radius = 0;

	Kernel::ID IDs[Kernel::MAX_MATERIALS];
	for(unsigned mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		IDs[mt] = Domains::global()->PM()->getID(mt, sID);

	for(Domains::DefectIterator it = Domains::global()->beginDI(); it!= Domains::global()->endDI(); ++it)
	{
		Kernel::M_TYPE mt = it->getElement()->getMaterial();
		const OKMC::Cluster *pMC = dynamic_cast<const OKMC::Cluster *>(*it);
		bool isMt = orig_mt == mt;
		bool isDefect = pMC && it->getDomain()->_pClPar->getParams(mt, pMC->getEDType())->_name == defect;
		bool isID = IDs[mt] == it->getID();
		if( (orig_mt == Kernel::MAX_MATERIALS || isMt)   && pMC &&
			(defect.empty()   || isDefect) &&
			(sID.empty()      || isID))
		{
			double rad = pMC->getRadius();
			if(rad >= minRad)
			{
				count++;
				radius += rad;
			}
		}
	}
	return (count ? radius / count : 0);
}

//static
//given a defect list, it fills it into a vector with the defect types.
void SimData::createDefectList(Kernel::M_TYPE mt, const string &def, vector<unsigned> &defects)
{
	vector<string> texts;
	IO::ParameterManager::getTokens(def, ',', texts);
	IO::ParameterManager * pPM = Domains::global()->PM();
	for(vector<string>::iterator it=texts.begin(); it!=texts.end(); ++it)
	{
		int defect = pPM->getEDType(mt, *it);
		if(defect == -1)
			ERRORMSG("defect name '" << *it << "'not recognized in material" << pPM->getMaterialName(mt));
		defects.push_back(defect);
	}
}

//returns histograms
//for mean values, a tcl script can be used
string SimData::getHistogram(Kernel::M_TYPE orig_mt, const string &orig_defect) const
{
	map<string, unsigned> histogram;

	vector<unsigned> defects;
	createDefectList(orig_mt, orig_defect, defects);
	for(Domains::DefectIterator it = Domains::global()->beginDI(); it!= Domains::global()->endDI(); ++it)
	{
		Kernel::M_TYPE mt = it->getElement()->getMaterial();
		const OKMC::Cluster *pMC = dynamic_cast<const OKMC::Cluster *>(*it);

		if( pMC && orig_mt == mt && std::find(defects.begin(), defects.end(), pMC->getEDType()) != defects.end())
		{
			Kernel::ID id = pMC->getID();
			string name = pMC->name();
			if(histogram.find(name) == histogram.end())
				histogram[name] = 1;
			else
				histogram[name]++;
		}
	}

	std::stringstream ss;
	for(map<string, unsigned>::const_iterator it=histogram.begin(); it!=histogram.end(); ++it)
		ss << it->first << " " << it->second << std::endl;
	return ss.str();
}

// Returns the dose implanted in the material as #FP / #atoms
// Count only the V's
double SimData::getDose(Kernel::M_TYPE mt) 
{
    IO::ParameterManager * pPM = Domains::global()->PM();
    map<string, unsigned> theMap = Domains::global()->client()->getCreatedDefects(mt);
    unsigned nFP(0); 
    unsigned nAtm(0);
    for(map<string, unsigned>::const_iterator it = theMap.begin(); it != theMap.end(); ++it)   
    {
        if(pPM->getDefectType(it->first, mt) == Kernel::Event::MOBILEPARTICLE && 
               it->first == "V")                    
                nFP += it->second;
        else
        {
            Kernel::ID theID = pPM->getID(mt, it->first);
            for(map<Kernel::P_TYPE, unsigned>::const_iterator pIt = theID._pt.begin();
                    pIt != theID._pt.end(); ++pIt)            
                if(pPM->getParticleName(mt, pIt->first) == "V")
                    nFP += it->second * pIt->second;        
        }
    }
    for(Domains::MeshElementIterator it = Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
        nAtm += it->getMaterial() == mt ? it->getAtoms() : 0;
    return double(nFP) / double(nAtm);
}

//obtains macroscopic diffusivity
double SimData::getDiffusivity(const string &mate, const string &name, bool bMacro) const
{
	Kernel::M_TYPE mt = Domains::global()->PM()->getMaterialNumber(mate);
	if(mt == Kernel::MAX_MATERIALS)
		ERRORMSG("Material " << mate << " not recognized.");
	Kernel::P_TYPE pt = Domains::global()->PM()->getParticleNumber(mt, name);
	if(pt == Kernel::UNDEFINED_TYPE)
		ERRORMSG("Particle name not recognized");
	if(pt < 2)
		ERRORMSG("I and V particles not allowed");
	if(Domains::global()->PM()->getMaterial(mt)._binary == false && bMacro)
		pt = Domains::global()->PM()->getFamily(pt);
	if(bMacro && pt >= Domains::global()->PM()->getNFamilies())
		ERRORMSG("Species needs to be simple when specifying macroscopic");


	double r2 = 0;
	unsigned howMany = 0;
	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		if(mt != it->getMaterial())
			continue;
		OKMC::Particle *part = it->getFirstPart();
		while(part)
		{
			if( (bMacro && Domains::global()->PM()->getFamily(part->getPType()) == pt) ||
					pt == part->getPType())
			{
				Kernel::Coordinates r=part->getCoordinates() - part->getOrig();
				r2 += r*r;
				howMany++;
			}
			part = part->getNext();
		}
	}
	double time = Domains::global()->data()->getDeltaTime();
	if(time && howMany)
		return 1./6.*r2*1e-14/howMany/time;
	else
	{
		WARNINGMSG("No time passed or no particles found for diffusion calculation, returning 0");
		return 0;
	}
}

IO::OutDataVectorC<double> SimData::getJumps(const string &name) const
{
	std::vector<std::string> tokens;
	IO::ParameterManager *PM=Domains::global()->PM();
	PM->getTokens(name,'_',tokens);
	Kernel::P_TYPE names[Kernel::MAX_MATERIALS];
	for(unsigned mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		names[mt] = Domains::global()->PM()->getParticleNumber(mt, name);
	IO::OutDataVectorC<double> odv3;
	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		Kernel::M_TYPE mt = (*it)->getMaterial();
		Kernel::P_TYPE pt = names[mt];
		unsigned state = PM->getStateNumber(mt, pt, tokens[1] );
		Coordinates m,M;
		it->getCorners(m, M);
		odv3.push(3, m, M, double(it->getHops(pt, state)) );
	}
	return collapse(odv3);
}

//takes a 3D information field and "collapses" it into 2D or 1D
//Assumes that the total mesh is a uniform tensor mesh even if there are several domains
IO::OutDataVectorC<double> SimData::collapse(const IO::OutDataVectorC<double> &odv_temp) const
{
	IO::OutDataVectorC<double> odv_c;
	IO::OutDataVectorC<double> odv;

	//filter for min, max
	for(IO::OutDataVectorC<double>::const_iterator it = odv_temp.begin(); it!=odv_temp.end(); ++it)
		if(it->_cm.isIntoEqual(_min, _max) && it->_cM.isIntoEqual(_min, _max))
			odv.push_back(*it);

	if(_dim == 3)
		return odv;

	float eps = 1e-2;
	if(_dim == 2) //collapse z
	{
		//stores all x and all y
		vector<float> linesX = Domains::global()->getDomain(0)->_pMesh->getLines(0);
		vector<float> linesY = Domains::global()->getDomain(0)->_pMesh->getLines(1);

		for(vector<float>::iterator itx = linesX.begin(); itx != linesX.end()-1; ++itx)
			for(vector<float>::iterator ity = linesY.begin(); ity != linesY.end()-1; ++ity)
			{
				double howMany = 0;
				double accumulated = 0;
				for(IO::OutDataVectorC<double>::const_iterator it=odv.begin(); it!=odv.end(); ++it)
					if( it->_cm._x < *itx + eps  && it->_cm._x > *itx - eps &&
						it->_cm._y < *ity + eps  && it->_cm._y > *ity - eps)
					{
						howMany++;
						accumulated += it->_data;
					}
				odv_c.push(2, Coordinates(*itx, *ity, 0), Coordinates(*(itx+1), *(ity+1), 0), howMany? accumulated/howMany:0);
			}
	}
	if(_dim == 1)
	{
		vector<float> linesX = Domains::global()->getDomain(0)->_pMesh->getLines(0);
		for(vector<float>::iterator itx = linesX.begin(); itx != linesX.end()-1; ++itx)
		{
			double howMany = 0;
			double accumulated = 0;
			for(IO::OutDataVectorC<double>::const_iterator it=odv.begin(); it!=odv.end(); ++it)
			{
				if( it->_cm._x < *itx + eps  && it->_cm._x > *itx - eps)
				{
					howMany++;
					accumulated += it->_data;
				}
			}
			odv_c.push(1, Coordinates(*itx, 0, 0), Coordinates(*(itx+1), 0, 0),  howMany? accumulated/howMany:0);
		}
	}
	//dim == 0.
	if(_dim == 0)
	{
		double howMany = 0;
		double accumulated = 0;
		for(IO::OutDataVectorC<double>::const_iterator it=odv.begin(); it!=odv.end(); ++it)
		{
			howMany++;
			accumulated += it->_data;
		}
		odv_c.push(0, Coordinates(), Coordinates(), howMany? accumulated/howMany:0);
	}
	return odv_c;
}

//takes a 3D particle field and "collapses" it into 3,2,1 or 0D creating concentrations.
IO::OutDataVectorC<double> SimData::collapse(const IO::OutDataVectorP<double> &odv_temp) const
{
	IO::OutDataVectorC<double> odv_out;
	IO::OutDataVectorP<double> odv;

	unsigned n[3];
	n[0] = round((_max._x - _min._x) / _odelta[0]);
	n[1] = round((_max._y - _min._y) / _odelta[1]);
	n[2] = round((_max._z - _min._z) / _odelta[2]);
	float delta[3] = {  (_max._x - _min._x)/n[0], (_max._y - _min._y)/n[1], (_max._z - _min._z)/n[2] };
	double factor = 0;

	for(IO::OutDataVectorP<double>::const_iterator it = odv_temp.begin(); it!=odv_temp.end(); ++it)
			if(it->_coords.isIntoEqual(_min, _max))
				odv.push_back(*it);

	if(_dim == 3)
	{
		factor = 1e21/(delta[0]*delta[1]*delta[2]);
		for(unsigned nx=0; nx < n[0]; ++nx)
			for(unsigned ny=0; ny < n[1]; ++ny)
				for(unsigned nz=0; nz < n[2]; ++nz)
					odv_out.push(_dim,
							Coordinates(_min._x + nx    *delta[0], _min._y + ny    *delta[1], _min._z + nz    *delta[2]),
							Coordinates(_min._x + (nx+1)*delta[0], _min._y + (ny+1)*delta[1], _min._z + (nz+1)*delta[2]), 0.);
	}
	if(_dim == 2)
	{
		factor = 1e21/(delta[0]*delta[1]*(_max._z - _min._z));
		for(unsigned nx=0; nx < n[0]; ++nx)
			for(unsigned ny=0; ny < n[1]; ++ny)
				odv_out.push(_dim, Coordinates(_min._x + nx    *delta[0], _min._y + ny    *delta[1], 0),
						           Coordinates(_min._x + (nx+1)*delta[0], _min._y + (ny+1)*delta[1], 0), 0.);
	}
	if(_dim == 1)
	{
		factor = 1e21/(delta[0]*(_max._y - _min._y)*(_max._z - _min._z));
		for(unsigned nx=0; nx < n[0]; ++nx)
			odv_out.push(_dim, Coordinates(_min._x + nx    *delta[0], 0, 0),
							   Coordinates(_min._x + (nx+1)*delta[0], 0, 0), 0.);
	}
	if(_dim == 0)
	{
		factor = 1e21/((_max._x - _min._x)*(_max._y - _min._y)*(_max._z - _min._z));
		odv_out.push(_dim, Coordinates(), Coordinates(), 0.);
	}
	for(IO::OutDataVectorP<double>::const_iterator it=odv.begin(); it!=odv.end(); ++it)
		odv_out[toODVC(delta, n, it->_coords)]._data += it->_data * factor;

	return odv_out;
}

unsigned SimData::toODVC(float delta[3], unsigned n[3], const Coordinates &c) const
{
	unsigned ix = (c._x - _min._x) / delta[0];
	unsigned iy = (c._y - _min._y) / delta[1];
	unsigned iz = (c._z - _min._z) / delta[2];
	switch(_dim)
	{
	case 3:
		return ix*n[2]*n[1] + iy*n[2] + iz;
	case 2:
		return ix*n[1] + iy;
	case 1:
		return ix;
	case 0:
		return 0;
	default:
		ERRORMSG("Incorrect dimension!");
		return 0;
	}
}

IO::OutDataVectorC<double> SimData::getProfileMobile(const string &name, const string &st, const string &mate) const
{

	IO::ParameterManager *PM=Domains::global()->PM();
	double time = Domains::global()->data()->getDeltaTime();
	IO::OutDataVectorC<double> odv;
	bool bAll = st.empty();
	Kernel::M_TYPE mt = Domains::global()->PM()->getMaterialNumber(mate);
	bool bParticleDefined = false;
	Kernel::P_TYPE names[Kernel::MAX_MATERIALS];
	for(unsigned mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		names[mt] = Domains::global()->PM()->getParticleNumber(mt, name);
	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		if (mate.empty()) mt = it->getMaterial();
		if(mt != it->getMaterial())
			continue;
		Kernel::P_TYPE pt = names[it->getMaterial()];
		if(pt == Kernel::UNDEFINED_TYPE || Domains::global()->PM()->isGasLike(mt)) //in gas or similar..
			continue;
		unsigned state = bAll? 0: PM->getStateNumber(mt, pt, st);
		if(state == Kernel::UNDEFINED_STATE)
			continue;
		bParticleDefined = true;
		Coordinates coordm ,coordM;
		it->getCorners(coordm, coordM);
		unsigned endDoWhile = (bAll? PM->getStates(mt, pt): state);
		double temp  = 0;
		do
		{
			double freq = Domains::global()->PM()->getMigrationRate(*it, pt, state);
			unsigned jumps = it->getHops(pt, state);
			temp += jumps? double(jumps)/(time*freq*it->getVolume()*1e-21) : 0.;
			state++;
		} while(state < endDoWhile);
		odv.push(3, coordm, coordM, temp);
	}
	if(bParticleDefined == false)
		ERRORMSG("Particle " << name << " does not seem to be defined in this materials...");
	return collapse(odv);
}

IO::OutDataVectorC<double> SimData::getParticleProfile(const string &needle, const string &def,
		const string &st, const string &mate) const
{
	IO::OutDataVectorP<double> odvp;
	IO::OutDataVectorC<double> odvc;

	std::string name = needle;
	enum class Mode:uint8_t {
		Individual   = 0u,
		AllInActive  = 1u,
		AllInCluster = 2u
	} mode = Mode::Individual;
	if(needle.back() == '*') {
		mode = Mode::AllInActive;
		name.pop_back();
	}
	else if(needle.back() == '@') {
		mode = Mode::AllInCluster;
		name.pop_back();
	}
	else {}

	Kernel::M_TYPE mt = Domains::global()->PM()->getMaterialNumber(mate);
	Kernel::P_TYPE names[Kernel::MAX_MATERIALS];
	for(unsigned mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		names[mt] = Domains::global()->PM()->getParticleNumber(mt, name);

	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		if (!mate.size()) mt = it->getMaterial();
		if(mt != it->getMaterial())
			continue;
		Kernel::P_TYPE pt = names[it->getMaterial()];
		if(mode == Mode::AllInCluster) {
			OKMC::Particle *pPart = it->getFirstPart();
			IO::ParameterManager const * const pm = Domains::global()->PM();
			Kernel::P_TYPE const family = pm->getFamily(pt);
			while(pPart)
			{
				OKMC::Cluster const* pCluster = dynamic_cast<OKMC::Cluster*>(pPart->getDefect());
				OKMC::Interface *pIf = dynamic_cast<OKMC::Interface *>(pPart->getDefect());
				if(pIf != nullptr) {       // We consider the trapped amount as inactive.
					std::vector<OKMC::Particle *> const particles = pIf->getParticles();
					for(OKMC::Particle const * const pPartInIf : particles) {
						OKMC::Cluster *pMcInIf = dynamic_cast<OKMC::Cluster *>(pPartInIf->getDefect());
						if(family == pm->getFamily(pPartInIf->getPType()) && pMcInIf == nullptr) {
							odvp.push(3, pPartInIf->getCoordinates(), 1.);
						}
					}
				}
				else if(pCluster != nullptr) {
					Kernel::ID const theMap = pCluster->getID();
					for(auto const item : theMap._pt) {
						if(family == pm->getFamily(item.first)) {
std::cout << pm->getIDName(theMap) << ' ' << pPart->getCoordinates() << '\n';
							odvp.push(3, pPart->getCoordinates(), item.second);
						}
					}
					odvp.push(3, pPart->getCoordinates(), 1.);
				}
				pPart = pPart->getNext();
			}
		}
		else if(pt == Kernel::UNDEFINED_TYPE) {   // amount of a specific cluster, not the dopants in it
			OKMC::Particle *pPart = it->getFirstPart();
			while(pPart)
			{
				OKMC::Cluster const* pCluster = dynamic_cast<OKMC::Cluster*>(pPart->getDefect());
				if(pCluster != nullptr &&
	    			Domains::global()->PM()->getIDName(pCluster->getID()) == name &&
					it->getDomain()->_pClPar->getParams(mt, pCluster->getEDType())->_name == def) {
					odvp.push(3, pPart->getCoordinates(), 1.);
				}
				pPart = pPart->getNext();
			}
		}
		else {   // name refers a valid mobile particle
			if (!Domains::global()->PM()->isImpurity(mt,pt) && Domains::global()->PM()->isAmorphous((*it)->getMaterial())) //No counting Is nor Vs for amorhpus MEs
				continue;
			Kernel::P_TYPE alloy_pt = Domains::global()->PM()->getMaterial((it->getMaterial()))._alloy;
			if(Domains::global()->PM()->getMaterial(mt)._binary == false && alloy_pt != Kernel::UNDEFINED_TYPE)
				alloy_pt = Domains::global()->PM()->getParticle(mt, alloy_pt, Kernel::POS_0);
			if(alloy_pt == pt)
			{
				Coordinates m, M;
				it->getCorners(m, M);
				odvc.push(3, m, M, float(it->getBAtoms()) / float(it->getVolume()) * 1e21);
	//			odvc.push(3, m, M, it->getAlloyFraction() * Domains::global()->PM()->getMaterial(it->getMaterial())._densityAlloyCm3);
			}
			else if(name.empty())
			{
				Coordinates m, M;
				it->getCorners(m, M);
				odvc.push(3, m, M, float(it->getAAtoms()) / float(it->getVolume()) * 1e21);
			}
			else
			{
				OKMC::Particle *pPart = it->getFirstPart();

				if(mode == Mode::Individual) {
					unsigned stN = Domains::global()->PM()->getStateNumber(it->getMaterial(), pt, st);
					while(pPart)
					{
						Kernel::Event::E_TYPE ev = pPart->getDefect()->getEType();
						OKMC::Cluster *pMC = dynamic_cast<OKMC::Cluster *>(pPart->getDefect());
						OKMC::Interface *pIF = dynamic_cast<OKMC::Interface *>(pPart->getDefect());
						Kernel::P_TYPE thisPt = pPart->getPType();
						if(Domains::global()->PM()->getMaterial(mt)._binary == false && (pMC || pIF))
							thisPt = Domains::global()->PM()->getParticleForCluster(mt, thisPt, Kernel::POS_0);
						bool isPart = thisPt == pt;
						bool isGenericDefect  = Kernel::Event::getEName(ev) == def;
						bool isExtendedDefect = pMC && it->getDomain()->_pClPar->getParams(mt, pMC->getEDType())->_name == def;
						bool isSt = pPart->getDefect()->getState() == stN;
						if(isPart &&
						(def.empty() || isGenericDefect || isExtendedDefect) &&
						(st.empty()  || isSt)) {
							odvp.push(3, pPart->getCoordinates(), 1.);
						}
						pPart = pPart->getNext();
					}
				}
				else {   // Mode::AllInActive
					IO::ParameterManager const * const pm = Domains::global()->PM();
					Kernel::P_TYPE const family = pm->getFamily(pt);
					if(family != 0) {
						while(pPart)
						{
							OKMC::Cluster *pMC = dynamic_cast<OKMC::Cluster *>(pPart->getDefect());
							OKMC::Interface *pIf = dynamic_cast<OKMC::Interface *>(pPart->getDefect());
							if(family == pm->getFamily(pPart->getPType()) && pMC == nullptr && pIf == nullptr) {
								odvp.push(3, pPart->getCoordinates(), 1.);
							}
							pPart = pPart->getNext();
						}
					}
				}
			}
		}
	}
	return odvc.size()? collapse(odvc) : collapse(odvp);
}

IO::OutDataVectorC<double> SimData::getAtomProfile(const string &name, const string &mate) const
{
	IO::OutDataVectorC<double> odv;

	Kernel::M_TYPE mt = Domains::global()->PM()->getMaterialNumber(mate);

	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		if (!mate.size()) mt = it->getMaterial();
		if(mt != it->getMaterial())
			continue;
		Coordinates m, M;
		it->getCorners(m, M);
		if(name == "A.atoms")
			odv.push(3, m, M, it->getAAtoms());
		else
			odv.push(3, m, M, it->getBAtoms());
	}
	return collapse(odv);
}

IO::OutDataVectorC<double> SimData::getBalance(const string &name, const string &mate) const
{
	IO::OutDataVectorC<double> odv;

	Kernel::M_TYPE mt = Domains::global()->PM()->getMaterialNumber(mate);

	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		if (!mate.size()) mt = it->getMaterial();
		if(mt != it->getMaterial())
			continue;
		Coordinates m, M;
		it->getCorners(m, M);
		if(name == "A.atoms")
			odv.push(3, m, M, it->_ABalance);
		else
			odv.push(3, m, M, it->_BBalance);
	}
	return collapse(odv);
}

IO::OutDataVectorC<double> SimData::getLKMCProfile(const string &name) const {
	IO::OutDataVectorP<double> odv;

	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		LKMC::LatticeAtom *pLA = static_cast<LKMC::LatticeAtom *>(it->getFirstLS());

		while (pLA) {
			/*if (name == "lkmc.defects" && pLA->getDefective()) {
				odv.push(3, pLA->getCoordinates(), 1.);
			}*/
			if (name == "lkmc.ac" && !pLA->getPerformed()) {
				for (unsigned i = 0; i < pLA->first(); ++i) {
					if (pLA->getNeighbor(i) && pLA->getNeighbor(i)->getPerformed()) {
						odv.push(3, pLA->getCoordinates(), 1.);
						break ;
					}
				}
			}
			pLA = static_cast<LKMC::LatticeAtom *>(pLA->getNext());
		}
	}
	return collapse(odv);
}

IO::OutDataVectorC<double> SimData::getStrain(const string &name) const
{
	IO::OutDataVectorC<double> odv;
	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		Coordinates m, M;
		it->getCorners(m, M);
		if(name == "xx")
			odv.push(3, m, M, it->strain_xx());
		else if(name == "yy")
			odv.push(3, m, M, it->strain_yy());
		else if(name == "zz")
			odv.push(3, m, M, it->strain_zz());
		else if(name == "xy")
			odv.push(3, m, M, it->strain_xy());
		else if(name == "xz")
			odv.push(3, m, M, it->strain_xz());
		else if(name == "yz")
			odv.push(3, m, M, it->strain_yz());
	}
	return collapse(odv);
}

IO::OutDataVectorC<double> SimData::getStress(const string &name) const
{
	IO::OutDataVectorC<double> odv;
	for(Domains::MeshElementIterator it = Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
	{
		Coordinates m, M;
		it->getCorners(m, M);
		if(name == "xx")
			odv.push(3, m, M, it->stress_xx());
		else if(name == "yy")
			odv.push(3, m, M, it->stress_yy());
		else if(name == "zz")
			odv.push(3, m, M, it->stress_zz());
		else if(name == "xy")
			odv.push(3, m, M,  it->stress_xy());
		else if(name == "xz")
			odv.push(3, m, M, it->stress_xz());
		else if(name == "yz")
			odv.push(3, m, M,  it->stress_yz());
	}
	return collapse(odv);
}

IO::OutDataVectorC<double> SimData::getElectrostatic(const string &name) const
{
	IO::OutDataVectorC<double> odv;
	for(Domains::MeshElementIterator it = Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
	{
		Coordinates m, M;
		it->getCorners(m, M);
		if(name == "potential") {
			odv.push(3, m, M, it->electrostaticPotential());
		}
		else if (name == "electron.density") {
			Kernel::MeshElement *pME = const_cast<Kernel::MeshElement *>(*it);
			std::set<Electrostatics::MeshNode *> s;

			Kernel::Mesh *pMesh = it->getDomain()->_pMesh;
			pMesh->getNodesFromElement(pME, s);

			if (s.empty()) {
				ERRORMSG("SimData::getElectrostatic: cannot find nodes of element " << it->getIndex());
			}

			double dat  = 0;
			for (std::set<Electrostatics::MeshNode *>::const_iterator i = s.begin(); i != s.end(); ++i) {
				dat  += (*i)->getElectronDensity();
			}
			odv.push(3, m, M, dat / static_cast<double>(s.size()));
		}
		else if (name == "hole.density") {
			Kernel::MeshElement *pME = const_cast<Kernel::MeshElement *>(*it);
			std::set<Electrostatics::MeshNode *> s;

			Kernel::Mesh *pMesh = it->getDomain()->_pMesh;
			pMesh->getNodesFromElement(pME, s);

			if (s.empty()) {
				ERRORMSG("SimData::getElectrostatic: cannot find nodes of element " << it->getIndex());
			}

			double dat  = 0;
			for (std::set<Electrostatics::MeshNode *>::const_iterator i = s.begin(); i != s.end(); ++i) {
				dat  += (*i)->getHoleDensity();
			}
			odv.push(3, m, M, dat / static_cast<double>(s.size()));
		}
		else if (name == "bandgap") {
			odv.push(3, m, M, it->bandGap());
		}
		else if (name == "conduction.bandedge") {
			// 0.5 * _Eg(i) - _psi(i)
			odv.push(3, m, M, 0.5 * it->bandGap() - it->electrostaticPotential());
		}
		else if (name == "valence.bandedge") {
			// - 0.5 * _Eg(i) - _psi(i)
			odv.push(3, m, M, - 0.5 * it->bandGap() - it->electrostaticPotential());
		}
	}
	return collapse(odv);
}

std::vector<float> SimData::getMixEnthalpy(const std::string &name) const
{
	std::vector<float> ov(1001);
	IO::FileParameters * pPar =  Domains::global()->getFileParameters();
	IO::Polynomial out = pPar->getPolynomial(name + "/Models/mixing.enthalpy");
	for(unsigned i = 0; i < ov.size(); i++)
	{
		float x = float(i) / 1000.;
		ov[i] = out.getValue(x);
	}
	return ov;
}

string SimData::getCoordination(double r, const Kernel::Coordinates &c, int typ) const
{
	Kernel::Coordinates m(-r-.5, -r-.5, -r-.5), M(r+.5,r+.5,r+.5), rel(0,0,0);
	m += c ;
	M += c;
	std::vector<LKMC::Lattice::LatticeInformation> neis; //type, coordinates
	if(!Domains::global()->client()->isInCell(c))
		ERRORMSG("extract: Coordinate " << c << " is not in this cell");
	Kernel::Domain *pDomain = Domains::global()->getDomain(c);
	Kernel::MeshElement *pEle = pDomain->_pMesh->getElement(pDomain->_pMesh->getIndexFromCoordinates(c));
	Kernel::M_TYPE mt = pEle->getMaterial();
	if(mt == Kernel::MAX_MATERIALS)
		ERRORMSG("Material " << Domains::global()->PM()->getMaterialName(pEle->getMaterial()) << " does not define a material to be transformed into");
	if(!pDomain->_pLat[mt])
		ERRORMSG("Material " << Domains::global()->PM()->getMaterialName(mt) << " does not have a defined lattice!");
	pDomain->_pLat[mt]->fill(m, M, neis, false);
	LOWMSG("Found " << neis.size() << " atoms\n");
	std::multimap<double, int> myMap;  //dist, atom_type
	double minDist = 1e57;
	int atomType = 0;
	for(std::vector<LKMC::Lattice::LatticeInformation>::iterator it=neis.begin(); it<neis.end(); ++it)
	{
		double dist = std::sqrt(
				(it->_coords._x - c._x)*(it->_coords._x - c._x) +
				(it->_coords._y - c._y)*(it->_coords._y - c._y) +
				(it->_coords._z - c._z)*(it->_coords._z - c._z));
		if(typ == -1 || typ == int(it->_type))
			if(dist < minDist)
			{
				minDist = dist;
				atomType = it->_type;
				rel = it->_coords;
			}
	}
	LOWMSG("Taking origin at atom " << rel << " of type " << atomType);
	LOWMSG("Dist   Typ H0  H1  H2");
	for(std::vector<LKMC::Lattice::LatticeInformation>::iterator it=neis.begin(); it<neis.end(); ++it)
	{
		double dist = std::sqrt(
				(it->_coords._x - rel._x)*(it->_coords._x - rel._x) +
				(it->_coords._y - rel._y)*(it->_coords._y - rel._y) +
				(it->_coords._z - rel._z)*(it->_coords._z - rel._z));
		if(dist > r)
			continue;
		myMap.insert(std::pair<double, int>(dist, it->_type));  //dist atom_type
	}
	int howMany[3] = { 0, 0, 0 };
	std::stringstream ss;
	for(std::multimap<double, int>::iterator it=myMap.begin(); it!=myMap.end(); ++it)
	{
		howMany[it->second]++;
		ss << std::right << setw(6) << std::setprecision(3) << it->first << " " << it->second << " " << setw(3) <<
				howMany[0] << " " << setw(3) << howMany[1] << " " << setw(3) << howMany[2] << std::endl;
	}
	return ss.str();
}

unsigned SimData::countDisplaced(const std::string &mate, const bool bAmorph, const float threshold) const
{
	IO::ParameterManager *PM = Domains::global()->PM();
	unsigned counter = 0, meCounter = 0;
	Kernel::M_TYPE mt = PM->getMaterialNumber(mate);
	string name = mate + "/Models/atomic.density";
	double amorph = Domains::global()->getFileParameters()->getFloat(name);
	float thrs = threshold > 0 ? threshold : 1.e24;

	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		if(PM->getMaterial(it->getMaterial())._basicMaterial == mt)
		{
			if(PM->isAmorphous(it->getMaterial()) )
				meCounter += unsigned(amorph/1.e21*it->getVolume());
			else
			{
				OKMC::Particle *pPart = it->getFirstPart();
				while(pPart) {

					Kernel::P_TYPE pt = pPart->getPType();
					if( PM->getFamily(pt) == PM->getMaterial(mt)._pt[0] || PM->getFamily(pt) == PM->getMaterial(mt)._pt[1] )
						meCounter++;
					pPart = pPart->getNext();
				}
			}
		}
		if (!bAmorph && double(meCounter)/it->getVolume()*1.e21 >= thrs)
			counter += unsigned(amorph/1.e21*it->getVolume());
		else
			counter += meCounter;
		meCounter = 0;
	}
	return counter;
}

IO::OutDataVectorC<double> SimData::amorphousFraction() const
{
	IO::ParameterManager *PM = Domains::global()->PM();
	IO::OutDataVectorC<double> odvc;
	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
		Coordinates m, M;
		it->getCorners(m, M);
		if( PM->isAmorphous(it->getMaterial()) )
			odvc.push(3,m,M,1.);
		else
			odvc.push(3,m,M,0.);
	}

	return collapse(odvc);
}

IO::OutDataVectorC<double> SimData::getDamageProfile() const
{
	IO::OutDataVectorC<double> odvc;
	IO::ParameterManager *PM = Domains::global()->PM();

	for(Domains::MeshElementIterator it=Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
	{
			Coordinates m, M;
			it->getCorners(m, M);
			if(Domains::global()->PM()->isAmorphous((*it)->getMaterial()))
				odvc.push(3,m,M,PM->getMaterial(PM->getMaterial(it->getMaterial())._basicMaterial)._amorphThreshold);
			else
			{
				OKMC::Particle *pPart = it->getFirstPart();
				unsigned nParts = 0;
				while(pPart) {
					Kernel::P_TYPE pt = pPart->getPType();
					if( !Domains::global()->PM()->isImpurity(it->getMaterial(),pt) )
						nParts++;
					pPart = pPart->getNext();
				}
				odvc.push(3,m,M,double(nParts)/it->getVolume()*1.e21);
			}
	}

	return collapse(odvc);
}

double SimData::getDeltaTime() const
{
	return Domains::global()->getTime() - _lastTime;
}

void SimData::reset()
{
	_lastTime = Domains::global()->getTime();
	// now reset the hops
	for(Domains::MeshElementIterator it = Domains::global()->beginMEI(); it!=Domains::global()->endMEI(); ++it)
		it.modify()->resetHops();
	for(Domains::ParticleIterator it = Domains::global()->beginPI(); it!=Domains::global()->endPI(); ++it)
		it.modify()->setOrig(it->getCoordinates());
}

//returns an y,z map with the free interface position.
string SimData::getFuzz() const
{
	std::stringstream ss;
	std::map<std::pair<float, float>, float> interface;

	Kernel::Coordinates cell, Cell;
	for(Domains::MeshElementIterator mi = Domains::global()->beginMEI(); mi != Domains::global()->endMEI(); ++mi)
	{
		mi->getCorners(cell, Cell);
		std::pair<float, float> yz(cell._y, cell._z);
		if( (interface.find(yz) == interface.end() || interface[yz] > cell._x) && mi->getInterfaces().size())
			interface[yz] = cell._x;
	}

	for(std::map<std::pair<float, float>, float>::iterator it=interface.begin(); it!=interface.end(); ++it)
		ss << it->first.first << " " << it->first.second << " " << it->second << std::endl;
	return ss.str();
}

}

