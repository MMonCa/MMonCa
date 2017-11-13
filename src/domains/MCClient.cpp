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

#include "MCClient.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include "io/FileParameters.h"
#include "okmc/MobileParticle.h"
#include "okmc/Cluster.h"
#include "okmc/Interface.h"
#include "kernel/Domain.h"
#include "kernel/RateManager.h"
#include "kernel/Mesh.h"
#include "lkmc/LKMCModel.h"

using std::string;
using std::stringstream;
using Kernel::Coordinates;
using Kernel::M_TYPE;

namespace Domains {

MCClient::MCClient(Tcl_Interp *p, const Coordinates &m, const Coordinates &M, const std::string &proc) :
_min(m), _max(M),
_rng2(Domains::global()->getFileParameters()->getInt("MC/General/random.seed")),
 _isFromStream(false)
{
	_getMaterial = new IO::GetMaterial(p, proc);
	beginInsert();
}

MCClient::MCClient(std::istream &is) :
		_rng2(Domains::global()->getFileParameters()->getInt("MC/General/random.seed")),
		_isFromStream(true)
{
	_getMaterial = new IO::GetMaterial(is);
	is >> _min >> _max;
	beginInsert();
}

MCClient::MCClient(const IO::MeshParser *mesh) :
	_rng2(Domains::global()->getFileParameters()->getInt("MC/General/random.seed")),
	_isFromStream(false)
{
	_getMaterial = new IO::GetMaterial(mesh);
	mesh->getCellSize(_min, _max);
	beginInsert();
}

void MCClient::restart(std::ostream &os) const
{
	IO::GetMaterial::restart(os);
	os << _min << " " << _max;
}

void MCClient::resetGetMaterial()
{
	delete _getMaterial;
	_getMaterial = 0;
}

MCClient::~MCClient()
{
}

//Since this is created from global domains, it uses the SubDomain associated with the coordinate
OKMC::Defect * MCClient::createMP(const std::string &name, const std::string &stat, const Kernel::Coordinates &c, bool bReact)
{
	if(c.isInto(_min, _max))
	{
		Kernel::Domain *pDomain = Domains::global()->getDomain(c);
		Kernel::MeshElement * pME = pDomain->_pMesh->getElement(pDomain->_pMesh->getIndexFromCoordinates(c));
		Kernel::M_TYPE mt = pME->getMaterial();
		Kernel::P_TYPE type = Domains::global()->PM()->getParticleNumber(mt, name);
		if(type == Kernel::UNDEFINED_TYPE)
		{
			WARNINGMSG("Discarded invalid particle. " << name);
			_discarded[name]++;
			return 0;
		}
		if(Domains::global()->PM()->getPPos(type) == Kernel::NO_POS && type != Domains::global()->PM()->getMaterial(mt)._alloy)
		{
			WARNINGMSG("Discarded particle with no valid position. " << name);
			_discarded[name]++;
			return 0;
		}
		Kernel::SubDomain *pSub = pDomain->_pRM->getSubDomain(pME->getSubDomainIdx());
		unsigned state = stat.empty()? 0 : Domains::global()->PM()->getStateNumber(mt, type, stat);
		if(state == Kernel::UNDEFINED_STATE)
		{
			_discarded[name]++;
			return 0;
		}
		if(!Domains::global()->PM()->isParticleDefined(type, mt) && (Domains::global()->PM()->getMaterial(mt)._alloy == Kernel::UNDEFINED_TYPE))
		{
			_discarded[name]++;
			return 0;
		}
		if(Domains::global()->PM()->isAmorphous(mt) && !Domains::global()->PM()->isImpurity(mt,type))
		{
			_discarded[name]++;
			return 0;
		}

		Kernel::P_TYPE fam = Domains::global()->PM()->getFamily(type);

		//this code contains many problems here, one of them is that the number in getAtoms is usually different than
		//the number or LKMC atoms. When calling profile, this is fixed later at the "profile" level and atom types
		//are overwritten with more accurate values, in the meanwhile, this is better than nothing...
		if(fam != Kernel::UNDEFINED_TYPE && fam == Domains::global()->PM()->getMaterial(mt)._alloy) //Alloy
		{
			if(pME->getAAtoms() == 0) // Not more B atoms allowed in this cell
			{
				_discarded[name]++;
				return 0;
			}
			_created[mt][name]++;
			pME->decAAtoms();
			pME->incBAtoms();
			char iv = Domains::global()->PM()->getIorV(mt, type);
			type = (iv == Kernel::IV_NONE? Kernel::UNDEFINED_TYPE : Domains::global()->PM()->iorv2pt(mt, iv, Kernel::POS_0)); //TODO: Binaries
			std::vector<LKMC::LatticeSite *> sites;
			LKMC::LatticeSite *pLS = pME->getFirstLS();
			while(pLS)
			{
				if(pLS->getPType() != fam)
					sites.push_back(pLS);
				pLS = pLS->getNext();
			}
			if(sites.size()) //change the type of one
			{
				pLS = sites[Domains::global()->client()->rand()*sites.size()];
				pLS->setPType(fam);
				LKMC::LatticeAtom *pLA = dynamic_cast<LKMC::LatticeAtom *>(pLS);
				if(pLA)
				{
					pLA->getDomain()->_pRM->update(pLA, pLA->getElement());
					for(unsigned i=0; i < pLA->third(); ++i)
					{
						LKMC::LatticeAtom *poLA = pLA->getNeighbor(i);
						if(poLA)
							poLA->getDomain()->_pRM->update(poLA, poLA->getElement());
					}
				}
			}

		}
		if(type != Kernel::UNDEFINED_TYPE)
		{
			_created[mt][name]++; // TODO: Change type into "substitutional" if creating interstitial dopant in amorphous
			OKMC::MobileParticle *mp= new OKMC::MobileParticle(pSub, type, state, pDomain, c, c);
			OKMC::Defect * resDef = mp;
			if(bReact)
			{
				std::vector<OKMC::Particle *> parts, partsRemove;
				pDomain->_pMesh->getInteractions(pSub, mp, parts);
				if(parts.size())
				{
					unsigned idxPart = int(pDomain->_rng_dom.rand() * parts.size());
					OKMC::Particle *intPar = parts[idxPart];
					parts[idxPart] = parts.back();
					parts.pop_back();
					resDef = intPar->getDefect()->interact(pSub, mp, intPar, partsRemove);
				}
			}
			if( !Domains::global()->PM()->isImpurity(mt,type)) // checking amorphization only for Is or Vs
				checkAmorphization(pSub, c);
			return resDef;
		}
	}
	else
		_out[name]++;

	return 0;
}

//Since this is created from global domains, it uses the SubDomain associated with the coordinate
OKMC::Cluster * MCClient::createMC(const std::string &name, const std::string &defect, const Kernel::Coordinates &c)
{
	if(c.isInto(_min, _max))
	{
		Kernel::Domain *pDomain = Domains::global()->getDomain(c);
		Kernel::MeshElement *pEle = pDomain->_pMesh->getElement(pDomain->_pMesh->getIndexFromCoordinates(c));
		Kernel::SubDomain *pSub = pDomain->_pRM->getSubDomain(pEle->getSubDomainIdx());
		M_TYPE mt = pEle->getMaterial();
		unsigned defNumber = pDomain->_pClPar->getDefectNumber(mt, defect);
		Kernel::ID theMap = Domains::global()->PM()->getID(mt, name);
		if(defNumber == pDomain->_pClPar->defectSize(mt) || theMap._pt.empty())
		{
			if(defNumber == pDomain->_pClPar->defectSize(mt))
				WARNINGMSG("Defect " << defect << " seems not to be defined in " << Domains::global()->PM()->getMaterialName(mt));
			else
				WARNINGMSG("Defect " << defect << " ID=" << name << " is not recognized in material " << Domains::global()->PM()->getMaterialName(mt));
			_discarded[name]++;
		}
		else
		{
			const OKMC::EDType::CLType * data = pDomain->_pClPar->getParams(defNumber, theMap);
			if(data)
			{
				OKMC::Cluster *pCl = new OKMC::Cluster(pSub, pDomain, defNumber, c);
				pCl = pCl->growDefect(pSub, theMap);
				if(pCl)
				{
					pDomain->_pRM->insert(pCl, pCl->getElement());
					_created[mt][name]++;
					return pCl;
				}
			}
			else
			{
				WARNINGMSG("Defect " << defect << " ID=" << name << " seems not to be defined in " << Domains::global()->PM()->getMaterialName(mt));
				_discarded[name]++;
			}
		}
	}
	else
		_out[name]++;
	return 0;
}

//Since this is created from global domains, it uses the SubDomain associated with the coordinate
OKMC::Interface * MCClient::createIF(const std::string &name, const Kernel::Coordinates &c)
{
	if(c.isInto(_min, _max))
	{
		Kernel::Domain *pDomain = Domains::global()->getDomain(c);
		Kernel::MeshElement * pME = pDomain->_pMesh->getElement(pDomain->_pMesh->getIndexFromCoordinates(c));
		Kernel::M_TYPE mt = pME->getMaterial();
		Kernel::P_TYPE type = Domains::global()->PM()->getFamily(Domains::global()->PM()->getParticleNumber(mt, name));
		if(type == Kernel::UNDEFINED_TYPE)
		{
			WARNINGMSG("Discarded invalid particle. " << name);
			_discarded[name]++;
			return 0;
		}

		OKMC::Interface *pInt = 0;
		for(std::vector<OKMC::Interface *>::iterator it=pME->getInterfaces().begin(); it!=pME->getInterfaces().end(); ++it)
		{
			pInt = *it;
			OKMC::Particle *pPart = new OKMC::Particle(type, c, pInt, c);
			pInt->insertParticle(pPart);
			_created[mt][name]++; // TODO: Change type into "substitutional" if creating interstitial dopant in amorphous
			break;
		}
		if(!pInt)
			_discarded[name]++;
		return pInt;
	}
	else
		_out[name]++;

	return 0;
}

void MCClient::beginInsert()
{
	_out.clear();
	_discarded.clear();
	for(Kernel::M_TYPE mt=0; mt < Kernel::MAX_MATERIALS; ++mt)
		_created[mt].clear();
}

void MCClient::endInsert()
{
	bool bPrint = false;

	stringstream out;
	out << "Out ";
	for(std::map<std::string, unsigned>::iterator it=_out.begin(); it!=_out.end(); ++it)
	{
		out << it->first << "(" << it->second << ") ";
		bPrint = true;
	}
	if(bPrint)
		LOWMSG(out.str());

	stringstream disc;
	bPrint = false;
	disc << "Discarded ";
	for(std::map<std::string, unsigned>::iterator it= _discarded.begin(); it!=_discarded.end(); ++it)
	{
		disc << it->first << "(" << it->second << ") ";
		bPrint = true;
	}
	if(bPrint)
		LOWMSG(disc.str());

	stringstream creat;
	bPrint = false;
	creat << "Created ";
	for(Kernel::M_TYPE mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		for(std::map<std::string, unsigned>::iterator it=_created[mt].begin(); it!=_created[mt].end(); ++it)
		{
			creat << Domains::global()->PM()->getMaterialName(mt) << "," <<
					it->first << "(" << it->second << ") ";
			bPrint = true;
		}
	LOWMSG(creat.str());
}

//checks if, and amorphizes (if needed), pME and MEs within a radius "dist" from c. It also amorphizes "lonely c-MEs" surrounded by a-MEs.
void MCClient::checkAmorphization(Kernel::SubDomain *pSub, const Coordinates &c)
{
	Kernel::Domain *pDomain = Domains::global()->getDomain(c);
	Kernel::Mesh *pMesh = pDomain->_pMesh;
	unsigned myIndex = pMesh->getIndexFromCoordinates(c);
	Kernel::MeshElement *pEle = pMesh->getElement(myIndex);
	const IO::MaterialProperties &mat = Domains::global()->PM()->getMaterial(pEle->getMaterial());
	IO::ParameterManager *pPM = Domains::global()->PM();

	float threshold = mat._amorphThreshold;
	if(mat._bAmorphous || threshold < 0 ) // amorphous materials cannot amorphize. threshold value is automatically assigned to -1 if it is not defined
		return;
	if(mat._amorphMaterial == Kernel::MAX_MATERIALS )
	{
		WARNINGMSG("It is inconsistent to define amorphous.threshold for " << Domains::global()->PM()->getMaterialName(pEle->getMaterial())
				<< " but not to define Amorphous" << Domains::global()->PM()->getMaterialName(pEle->getMaterial()) << ". Ignoring amorphization.");
		return;
	}

	float dist = 1., vol = 0., nParts = 0;
	unsigned crystElems = 0;
	std::vector<Kernel::MeshElement *> elems;
	pMesh->fillNeighborsAllMat(c,dist,elems);  //IMB: I'm not sure this is the best function to be called...
	for(std::vector<Kernel::MeshElement *>::iterator it=elems.begin(); it!=elems.end(); ++it)
	{
		vol += (*it)->getVolume();
		if(pPM->isAmorphous((*it)->getMaterial()) || pPM->isGasLike((*it)->getMaterial()))
			nParts += (*it)->getAmorphParts();
		else if(!Domains::global()->PM()->isAmorphous((*it)->getMaterial()))
		{
			crystElems++;
			OKMC::Particle *pPart = (*it)->getFirstPart();
			while(pPart)
			{
				Kernel::P_TYPE pt = pPart->getPType();
				if( !Domains::global()->PM()->isImpurity((*it)->getMaterial(),pt) )
					nParts++;
				pPart = pPart->getNext();
			}
		}
	}
	std::vector<Kernel::MeshElement *> toAmor;
	if(nParts/vol*1.e21 >= threshold)
	{
		for(std::vector<Kernel::MeshElement *>::iterator it=elems.begin(); it!=elems.end(); ++it)
			if(!pPM->isAmorphous((*it)->getMaterial()) && (*it)->getMaterial() == mat._basicMaterial)
			{
				float boxParts = 0;
				OKMC::Particle *pPart = (*it)->getFirstPart();
				while(pPart)
				{
					Kernel::P_TYPE pt = pPart->getPType();
					if( !Domains::global()->PM()->isImpurity((*it)->getMaterial(),pt) )
						boxParts+=1.;
					pPart = pPart->getNext();
				}
				assert(pPM->isAmorphous((*it)->getToMat(LKMC::MOD_SPER)) == true);
				(*it)->amorphizeME(pSub);
				(*it)->setAmorphParts(boxParts);
				toAmor.push_back(*it);
			}
		HIGHMSG("Amorphizing... Volume: " << vol << " Boxes: " << elems.size() << " threshold is "
				<< threshold << " nParts is " << nParts	<< " Coord x " << c._x);
	}
	else if(float(crystElems) < 0.12*float(elems.size()))
	{
		for(std::vector<Kernel::MeshElement *>::iterator it=elems.begin(); it!=elems.end(); ++it)
			if(!pPM->isAmorphous((*it)->getMaterial()) && (*it)->getMaterial() == mat._basicMaterial)
			{
				float boxParts = 0;
				OKMC::Particle *pPart = (*it)->getFirstPart();
				while(pPart)
				{
					Kernel::P_TYPE pt = pPart->getPType();
					if( !Domains::global()->PM()->isImpurity((*it)->getMaterial(),pt) )
						boxParts+=1.;
					pPart = pPart->getNext();
				}
				(*it)->setAmorphParts(boxParts);
				assert(pPM->isAmorphous((*it)->getToMat(LKMC::MOD_SPER)) == true);
				(*it)->amorphizeME(pSub);
				toAmor.push_back(*it);
			}
		HIGHMSG("Amorphizing " << crystElems << " 'lonely' elements. Because there were " << elems.size() - crystElems
				<< " surrounding amorphous elements");
	}
	if(toAmor.size())
		pDomain->_pLKMC->localSPER(pSub, toAmor);

}

}
