/*
 * ParameterManager.cpp
 *
 *  Created on: Feb 24, 2011
 *
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

#include "ParameterManager.h"
#include "io/FileParameters.h"
#include "kernel/MeshElement.h"
#include "kernel/RateManager.h"
#include "kernel/Event.h"
#include "kernel/Mesh.h"
#include "okmc/ClusterParam.h"

#include <boost/lexical_cast.hpp>

#include <cstring>
#include <cassert>
#include <cstdlib>

using std::map;
using std::vector;
using std::string;
using std::pair;
using std::stringstream;
using Kernel::Event;
using Kernel::P_TYPE;
using Kernel::P_POS;
using Kernel::UNDEFINED_TYPE;
using Kernel::M_TYPE;

namespace IO {

ParameterManager::FamilyProperties::FamilyProperties(const string &f, const string &d)
: _file(f), _s_family(d)
{
	for(unsigned i=0; i<5; ++i)
		_type[i] = Kernel::UNDEFINED_TYPE;
}

ParameterManager::ParameterManager(Tcl_Interp *pTcl, const IO::FileParameters *p)
{
	for(unsigned mt=0; mt < Kernel::MAX_MATERIALS; ++mt)
		_maxCapRad[mt] = 0;
	readMaterials(pTcl, p);
	readElements(p);
	readDefects  (pTcl, p);
	setAlloy(p);
}

ParameterManager::~ParameterManager()
{
}

void ParameterManager::readElements(const IO::FileParameters *p)
{
	map<string, string> elements = p->getStringMap("MC/Particles/elements");
	_families.push_back(FamilyProperties("Vacancy", "V"));
	_families.back()._type[Kernel::NO_POS] = 0;
	_particles.push_back(ParticleProperties(0, Kernel::NO_POS));
	for(map<string, string>::iterator it=elements.begin(); it != elements.end(); ++it)
	{
		vector<string> tokens;
		string value = it->second;
		getTokens(value, ',', tokens);
		if(tokens.size() != 3)
			ERRORMSG(value <<  " incorrect in MC/Particles/elements. 3 tokens needed!");
		unsigned element = atoi(tokens[1].c_str());
		if(element > Kernel::MAX_ELEMENTS)
			ERRORMSG("Wrong element number " << tokens[2] << " in MC/Particles/elements");
		float mass = atof(tokens[2].c_str());
		_elements[element]._element = element;
		_elements[element]._mass = mass*1.66053892e-27; //kg.
		_elements[element]._name = tokens[0];
		_elements[element]._mt = Kernel::MAX_MATERIALS;

		//build family table.
		_families.push_back(FamilyProperties(tokens[0], it->first));
		_families.back()._type[Kernel::NO_POS] = _families.size() -1;
		_particles.push_back(ParticleProperties(_particles.size(), Kernel::NO_POS));
		if(_families.size() > Kernel::MAX_IMPURITIES)
			ERRORMSG("Maximum number of particles " << int(Kernel::MAX_IMPURITIES) << " exceeded in MC/Particles/families");
	}
	//insert Vs
	_families[0]._type[Kernel::POS_0] = _particles.size();
	_particles.push_back(ParticleProperties(Kernel::V_TYPE, Kernel::POS_0));
	_families[0]._type[Kernel::POS_1] = _particles.size();
	_particles.push_back(ParticleProperties(Kernel::V_TYPE, Kernel::POS_1));
	_families[0]._type[Kernel::POS_I] = _particles.size();  //especial case for clusters
	_particles.push_back(ParticleProperties(Kernel::V_TYPE, Kernel::POS_I));
	//insert particles defining related impurities
	for(vector<FamilyProperties>::iterator it=_families.begin()+1; it!=_families.end(); ++it)
	{
		unsigned fam = it->_type[Kernel::NO_POS];
		for(unsigned pos=0; pos < Kernel::NO_POS; ++pos)
		{
			it->_type[pos] = _particles.size();
			_particles.push_back(ParticleProperties(fam, P_POS(pos)));
		}
	}
	//read which particles and states belong to which materials
	//init all to false
	for(std::vector<ParticleProperties>::iterator it=_particles.begin(); it!=_particles.end(); ++it)
	{
		it->_materials = std::vector<bool>(getNMaterials(), false);
		it->_states = std::vector<std::vector<std::string> >(getNMaterials());
	}
	//read them
	for(M_TYPE mt=0; mt<getNMaterials(); ++mt)
	{
		// read binary materials
		string binary = p->getString(getMaterialName(mt) + "/Models/material.composition");
		vector<string> composition;
		getTokens(binary, ',', composition);
		if(composition.size() == 1) //here because we need the materials definition
		{
			_materials[mt]._binary = false;
			_materials[mt]._s_pt[0] = composition[0];
			_materials[mt]._pt[0] = Kernel::UNDEFINED_TYPE;
			for(unsigned i=0; i < _families.size(); ++i)
				if(composition[0] == _families[i]._s_family)
					_materials[mt]._pt[0] = i;
			if(_materials[mt]._pt[0] == Kernel::UNDEFINED_TYPE && composition[0] != "none")
				ERRORMSG(getMaterialName(mt) + "/Models/material.composition. Cannot understand composition " << composition[0]);
			_materials[mt]._pt[1] = Kernel::UNDEFINED_TYPE;
		}
		else if(composition.size() == 2)
		{
			_materials[mt]._binary = true;
			for(int i=0; i<2; ++i)
			{
				_materials[mt]._s_pt[i] = composition[i];
				_materials[mt]._pt[i] = Kernel::UNDEFINED_TYPE;
				for(unsigned j=0; j < _families.size(); ++j)
					if(composition[i] == _families[j]._s_family)
						_materials[mt]._pt[i] = j;
				if(_materials[mt]._pt[i] == Kernel::UNDEFINED_TYPE)
					ERRORMSG(getMaterialName(mt) + "Models/material.composition. Cannot understand composition " << composition[i]);
			}
		}
		else
			ERRORMSG(getMaterialName(mt) + "/Models/material.composition. Incorrect number of compositions. " << binary);

		map<string,bool> parts = p->getBoolMap(getMaterialName(mt) + "/Models/particles");
		for(map<string,bool>::iterator it=parts.begin(); it != parts.end(); ++it)
		{
			if(!it->second)
				continue;
			std::vector<std::string> tokens,tokens2;
			getTokens(it->first, '_', tokens);
			P_TYPE pt = getParticleNumber(mt, tokens[0]);
			if(pt == UNDEFINED_TYPE)
				ERRORMSG(getMaterialName(mt) + "/particles" << ": Particle not recognized: " << tokens[0]);
			P_POS pos = getPPos(pt);
			if(pos == Kernel::NO_POS)
				ERRORMSG(getMaterialName(mt) + "/particles" << ": Position not recognized: " << tokens[0]);
			_particles[pt]._materials[mt] = true; //set this one.
			_particles[getFamily(pt)]._materials[mt] = true;
			if(tokens.size() == 1)   //add empty state
				_particles[pt]._states[mt].push_back("");
			else if(tokens.size()-1 > Kernel::MAX_STATES)
				ERRORMSG(getMaterialName(mt) + "/particles" << ": Too many states");
			else
			{
				getTokens(tokens[1],',',tokens2);
				for(unsigned i=0; i<tokens2.size(); ++i)
					_particles[pt]._states[mt].push_back(tokens2[i]);
			}
		}
	}
}

void ParameterManager::readMaterials(Tcl_Interp *pTcl, const IO::FileParameters *p)
{
	map<string, string> names = p->getStringMap("MC/General/materials");
	for(map<string, string>::iterator it=names.begin(); it!=names.end(); ++it)
	{
		MaterialProperties mp(it->first, it->second);
		_materials.push_back(mp);
		_materials.back()._bAmorphous = mp._name.find_first_of("Amorphous") == 0;
		MEDMSG("Material " << it->first << " (" << it->second << ") is number " << _materials.size()-1 << " " <<
				(_materials.back()._bAmorphous? " amorphous" : " crystalline"));
		_materials.back()._amorphMaterial = Kernel::MAX_MATERIALS;		
		_materials.back()._permittivity = p->getFloat(it->first + "/Models/permittivity");
		_materials.back()._molarVolume  = p->getFloat(it->first + "/Models/molar.volume");
		_materials.back()._bEpitaxy     = p->getBool (it->first + "/Models/epitaxy");
		_materials.back()._alloy 		= Kernel::UNDEFINED_TYPE;   //initialized here because setAlloy comes after,resulting in segmentation fault if not done.

                if(p->getFloat(it->first + "/Models/atomic.density") == 0)
                        ERRORMSG("Parameter 'atomic.density' in material " << it->first << " must be non-zero.");
                _materials.back()._densityCm3   = p->getFloat(it->first + "/Models/atomic.density");                                
		if(!_materials.back()._bAmorphous && p->specified(it->first + "/Models/amorphization.threshold"))
			_materials.back()._amorphThreshold = p->getFloat(it->first + "/Models/amorphization.threshold");
		else
			_materials.back()._amorphThreshold	= -1;

		_materials.back()._element = 0;
		for(unsigned i=0; i<Kernel::MAX_ELEMENTS; ++i)
		{
			if(_elements[i]._name == it->first)
			{
				MEDMSG("Recognizing element " << i << " for material " << it->first);
				_materials.back()._element = i;
				_elements[i]._mt = _materials.size()-1;
				break;
			}
		}

		// Electronic paramater: bandgap @ 300K
		double Eg0 = 0;
		if (p->specified(it->first + "/ElectronicStructure")) {
			std::string EgPath =  it->first + "/ElectronicStructure/Bandgap";
			std::string T0     = boost::lexical_cast<std::string>(300);

			p->loadProcedure(pTcl, EgPath, 1);
			Eg0 = p->getFloatProc(pTcl, EgPath, T0);
		}
		_materials.back()._Eg0 = Eg0;
		//
		_materials.back().E_                   = p->getFloat(it->first + "/Mechanics/E"); //Young's modulus
		_materials.back().nu_                  = p->getFloat(it->first + "/Mechanics/nu"); //Poisson's ratio
		_materials.back().thermalExp_          = p->getFloat(it->first + "/Mechanics/thermal.expansion"); //thermal expansion
		_materials.back().tempRef_             = p->getFloat(it->first + "/Mechanics/temp.ref"); //reference temperature
		if(_materials.back()._bAmorphous)
		{
			_materials.back().amorphousExpansion_  = p->getFloat(it->first + "/Mechanics/amorphous.expansion");
		}
		else
		{
			map<string, string> theMap = p->getStringMap(it->first + "/Mechanics/eigenstrains");
			for(map<string, string>::iterator sit=theMap.begin(); sit!=theMap.end(); ++sit)
			{
				vector<string> tokens;
				getTokens(sit->second, ',', tokens);
				if(tokens.size() != 2)
					ERRORMSG("Incorrect number of eigenstrain values in " << it->first + "/Mechanics/eigenstrains");
				_materials.back().eigenstrains_[sit->first] = pair<double, double>(atof(tokens[0].c_str()), atof(tokens[1].c_str()));
			}
		}
	}
	if(_materials.size() > Kernel::MAX_MATERIALS)
		ERRORMSG("Too many materials defined! Maximum is " << Kernel::MAX_MATERIALS);
	for(vector<MaterialProperties>::iterator it=_materials.begin(); it!=_materials.end(); ++it)
	{
		it->_basicMaterial  = (it->_bAmorphous? getMaterialNumber(it->_name.c_str() + strlen("Amorphous")): getMaterialNumber(it->_name));
		it->_amorphMaterial = (it->_bAmorphous? getMaterialNumber(it->_name) : getMaterialNumber(string("Amorphous") + it->_name));
	}

}

void ParameterManager::readDefects(Tcl_Interp *pTcl,  const IO::FileParameters *p)
{
	LOWMSG2("Reading defects: ");
	for(M_TYPE mt=0; mt < getNMaterials(); ++mt)
	{
		_defects.push_back(std::map<string, Kernel::Event::E_TYPE>());
		for(M_TYPE mt2=0; mt2 < getNMaterials(); ++mt2)
			if(mt2 != mt)
				_defects[mt][getMaterialName(mt2)] = Kernel::Event::INTERFACE;
		//first particles
		for(P_TYPE pt=0; pt < getNParticles(); ++pt)
			if(isParticleDefined(pt, mt) || pt == _materials[mt]._alloy)
				_defects[mt][getParticleName(mt, pt)] = Kernel::Event::MOBILEPARTICLE;
		//now extended defects
		map<string,bool> theArray = p->getBoolMap(getMaterialName(mt) + '/' + "Models/defined");
		LOWMSG2(getMaterialName(mt) + '(');
		for(map<string,bool>::iterator it=theArray.begin(); it!=theArray.end(); ++it)
		{
			if(!it->second)
				continue;
			if(it->first == "0")
				_defects[mt][it->first] = Kernel::Event::EMPTY;
			else
			{
				_defects[mt][it->first] = Kernel::Event::CLUSTER;
				_edDefectName[mt].push_back(it->first);
				p->loadProcedure(pTcl, getMaterialName(mt) + '/' + it->first + "/formation", 0);
				LOWMSG2(it->first << " ");
				map<string,IO::ArrheniusAlloys> potMap = p->getArrheniusAlloysProc(pTcl,
						getMaterialName(mt) + '/' + it->first + "/formation", Parameters::AA_ENERGY);
				for(map<string,IO::ArrheniusAlloys>::iterator sIt=potMap.begin(); sIt!=potMap.end(); ++sIt)
					_defects[mt][sIt->first] = Kernel::Event::CLUSTER;
			}
		}
		LOWMSG2(") ");
	}
	LOWMSG("");
}

void ParameterManager::setAlloy(const IO::FileParameters *p)
{
	for(vector<MaterialProperties>::iterator it=_materials.begin(); it!=_materials.end(); ++it)
	{
		M_TYPE mt = it - _materials.begin();
		string matName =  it->_name + "/Models/alloy";
		string matDen = it->_name + "/Models/alloy.density";
		string selfDiff = it-> _name + "/Models/self.diffusion";
		string alloy = p->getString(matName);
		it->_alloy = Kernel::UNDEFINED_TYPE;
		for(unsigned i=0; i < _families.size(); ++i)
			if(alloy == _families[i]._file || alloy == _families[i]._s_family)
				it->_alloy = i;
		if(alloy != "none")
		{
			if(it->_binary)
				ERRORMSG("Binary material " << it->_name << " cannot have alloys");
			if(it->_alloy == UNDEFINED_TYPE)
				ERRORMSG("Cannot recognize alloy name " << alloy);
			it->_alloyName = alloy;
			it->_densityAlloyCm3 = p->getFloat(matDen);
			it->_selfDiffusion = p->getBool(selfDiff);
			_defects[it-_materials.begin()][getParticleName(mt, it->_alloy)] = Kernel::Event::MOBILEPARTICLE;
			it->alloyE_          = p->getFloat(it->_name + "/Mechanics/alloy.E");
			it->alloynu_         = p->getFloat(it->_name + "/Mechanics/alloy.nu");
			it->alloythermalExp_ = p->getFloat(it->_name + "/Mechanics/alloy.thermal.expansion");
			it->alloytempRef_    = p->getFloat(it->_name + "/Mechanics/alloy.temp.ref");
			if(it->_alloy == Kernel::UNDEFINED_TYPE)
				ERRORMSG(matName << ": Cannot find " << alloy << " as a particle");
		}
		else
		{
			it->_densityAlloyCm3 = 0;
			it->_selfDiffusion =  false;
			it->_alloyName = "none";
		}
		assert(it->_alloy == Kernel::UNDEFINED_TYPE || it->_alloy < getNFamilies());
	}
}

void ParameterManager::getTokens(const string &key, char delimiter, vector<string> &tokens)
{
	tokens.clear();
	for(string::const_iterator it=key.begin(); it !=key.end(); ++it)
	{
		if(tokens.size() == 0)
			tokens.push_back("");
		if(*it == delimiter)
			tokens.push_back("");
		else
			tokens.back() += *it;
	}
}

M_TYPE ParameterManager::getMaterialNumber(const string &name) const
{
	for(vector<MaterialProperties>::const_iterator it=_materials.begin(); it!=_materials.end(); ++it)
		if(it->_name == name)
			return it - _materials.begin();
	return Kernel::MAX_MATERIALS;
}

//given a particle as a string provides the particle as a type
//accepts different formats for binary and non binary materials
P_TYPE ParameterManager::getParticleNumber(M_TYPE mt, const string &name) const
{
	//thinking, fast...
	map<string, P_TYPE>::const_iterator it = _partNameCache[mt].find(name);
	if(it != _partNameCache[mt].end())
		return it->second;

	//and slow...
	const bool binary = _materials[mt]._binary;
	if(binary)
		for(P_TYPE pt=0; pt < getNParticles(); ++pt)
		{
			P_TYPE fam = getFamily(pt);
			P_POS pos = getPPos(pt);
			if(fam == Kernel::V_TYPE && pos == Kernel::POS_I)
				continue;
			if (name == _families[fam]._s_family + getPositionName(mt, pos))
			{
				_partNameCache[mt][name] = pt;
				return pt;
			}
		}
	else
	{
		if(_materials[mt]._s_pt[0] == name)  //Fe in Fe, Si in Si...
		{
			_partNameCache[mt][name] = getParticleForCluster(mt, _materials[mt]._pt[0], Kernel::POS_0);
			return _partNameCache[mt][name];
		}
		for(P_TYPE pt=getNFamilies(); pt < getNParticles(); ++pt)
		{
			P_TYPE fam = getFamily(pt);
			P_POS pos = getPPos(pt);
			if(pos == Kernel::POS_1 || (fam == Kernel::V_TYPE && pos == Kernel::POS_I) )
				continue;
			if(fam == _materials[mt]._pt[0] && name == getPositionName(mt, pos)) //I
			{
				_partNameCache[mt][name] = pt;
				return pt;
			}
			if (name == _families[fam]._s_family + getPositionName(mt, pos))
			{
				_partNameCache[mt][name] = pt;
				return pt;
			}
		}                
	}
        return UNDEFINED_TYPE;
}

P_TYPE ParameterManager::getElementNumber(const std::string &text) const
{
	for(vector<FamilyProperties>::const_iterator it = _families.begin(); it!=_families.end(); ++it)
		if(it->_s_family == text)
			return it - _families.begin();
	return UNDEFINED_TYPE;
}

//given a string (like I, C or similar) provides the position number
//uses different formats for binary and non binary materials
P_POS ParameterManager::getPositionNumber(M_TYPE mt, const string &name) const
{
	bool binary = _materials[mt]._binary;
	if(name.empty() && binary)
		return Kernel::NO_POS;
	if(name.empty() && !binary)
		return Kernel::POS_0;
	if(name == "I")
		return Kernel::POS_I;
	if(name == _materials[mt]._s_pt[0])
		return Kernel::POS_0;
	if(name == _materials[mt]._s_pt[1])
		return Kernel::POS_1;
	if(name == "V")
		return Kernel::POS_V;
	return Kernel::UNDEFINED_POS;
}

//returns a particle name in readable format.
//format is SiI, CSi, BI, AsV, etc...
//different format for binary and non binary materials
std::string ParameterManager::getParticleName(M_TYPE mt, P_TYPE pt) const
{
	P_POS pos = getPPos(pt);
	P_TYPE family = getFamily(pt);
	bool binary = _materials[mt]._binary;
	if(!binary && pos != Kernel::NO_POS && family == _materials[mt]._pt[0])
		return getPositionName(mt, pos);
	return _families[family]._s_family + getPositionName(mt, pos);
}

P_TYPE ParameterManager::getParticle(M_TYPE mt, P_TYPE fam, P_POS pos) const
{
	P_TYPE pt0 = _materials[mt]._pt[0];
	P_TYPE pt1 = _materials[mt]._pt[1];
	if(pt1 == UNDEFINED_TYPE && pos == Kernel::POS_1)
		return UNDEFINED_TYPE;
	if( (fam == Kernel::V_TYPE && pos == Kernel::POS_I) ||
		(fam == pt0 && (pos == Kernel::POS_0 || pos == Kernel::POS_V)) ||  //Si^Si, SiV
	    (fam == pt1 && (pos == Kernel::POS_1 || pos == Kernel::POS_V))   ) //C^C, CV
		return UNDEFINED_TYPE;
	return _families[fam]._type[pos];
}

//returns particles that are not available as mobile particles (V^I, Si^Si...)
P_TYPE ParameterManager::getParticleForCluster(M_TYPE mt, P_TYPE fam, P_POS pos) const
{
	P_TYPE pt1 = _materials[mt]._pt[1];
	if(pt1 == UNDEFINED_TYPE && pos == Kernel::POS_1)
		return UNDEFINED_TYPE;
	return _families[fam]._type[pos];
}


// Antisite Si is Si^C
P_TYPE ParameterManager::getAntisite(M_TYPE mt, P_POS pos) const
{
	if(_materials[mt]._binary == false)
		return UNDEFINED_TYPE;
	unsigned i = (pos == Kernel::POS_0 ? 0:1);
	P_POS ipos = (pos == Kernel::POS_0 ? Kernel::POS_1: Kernel::POS_0);
	P_TYPE pt = _materials[mt]._pt[i];
	return _families[pt]._type[ipos];
}

P_TYPE ParameterManager::getPositionAsFamily(M_TYPE mt, P_POS pos) const
{
	switch(pos)
	{
	case Kernel::POS_0:
		return _materials[mt]._pt[0];
	case Kernel::POS_1:
		return _materials[mt]._pt[1];
	default:
		ERRORMSG("Incorrect position in getPositionAsFamily");
		return Kernel::UNDEFINED_TYPE;
	}
}

//different names for binary and non-binary materials
string ParameterManager::getPositionName(M_TYPE mt, P_POS pos) const
{
	bool binary = _materials[mt]._binary;
	switch(pos)
	{
	case Kernel::POS_0:
		return (binary ?(_families[_materials[mt]._pt[0]]._s_family) : "");
	case Kernel::POS_1:
		return (_families[_materials[mt]._pt[1]]._s_family);
	case Kernel::POS_I:
		return ("I");
	case Kernel::POS_V:
		return ("V");
	case Kernel::NO_POS:
		return ("");
	default:
		ERRORMSG("Incorrect position in getPositionName");
		return "";
	}

}

string ParameterManager::getRawPositionName(M_TYPE mt, P_POS pos) const
{
	switch(pos)
	{
	case Kernel::POS_0:
		return _families[_materials[mt]._pt[0]]._s_family;
	case Kernel::POS_1:
		return _families[_materials[mt]._pt[1]]._s_family;
	case Kernel::POS_I:
		return ("I");
	case Kernel::POS_V:
		return ("V");
	case Kernel::NO_POS:
		return ("");
	default:
		ERRORMSG("Incorrect position in getRawPositionName");
		return "";
	}

}

//return 1 for vacancies or vacancy pairs, and 0 for particles in interstitial position.
// UNDEFINED_TYPE for everything else.
char ParameterManager::getIorV(M_TYPE mt, P_TYPE pt) const
{
	P_TYPE fam = getFamily(pt);
	P_POS  pos = getPPos(pt);
	if(fam == Kernel::V_TYPE || pos == Kernel::POS_V)
		return Kernel::IV_V;
	if(getPPos(pt) == Kernel::POS_I)
		return Kernel::IV_I;
	return Kernel::IV_NONE;
}

//returns 0 or 1 only for SiI, CI or VSi VC
char ParameterManager::isIorV(M_TYPE mt, P_TYPE pt) const
{
	P_TYPE fam = getFamily(pt);
	if(fam == Kernel::V_TYPE)
		return Kernel::IV_V;
	P_POS  pos = getPPos(pt);
	if(pos == Kernel::POS_I && (fam == _materials[mt]._pt[0] || fam == _materials[mt]._pt[1]) )
		return Kernel::IV_I;
	return Kernel::IV_NONE;
}

//B, As, etc..., i.e., NO_POS and not a constituent material
bool ParameterManager::isPureImpurity(M_TYPE mt, P_TYPE pt) const
{
	if(pt == Kernel::V_TYPE || pt >= _families.size())
		return false;
	if(pt == _materials[mt]._pt[0] || pt == _materials[mt]._pt[1])
		return false;
	return true;
}

//B, As, BI, AsV etc..., i.e., not a constituent material
bool ParameterManager::isImpurity(M_TYPE mt, P_TYPE pt) const
{
	if(getFamily(pt) == Kernel::V_TYPE || getFamily(pt) == _materials[mt]._pt[0] || getFamily(pt) == _materials[mt]._pt[1] )
		return false;
	return true;
}


//returns the correct I or V
P_TYPE ParameterManager::iorv2pt(M_TYPE mt, unsigned char iorv, P_POS pos) const
{
	assert(pos == Kernel::POS_0 || pos == Kernel::POS_1);
	if(iorv == Kernel::IV_V)
		return _families[Kernel::V_TYPE]._type[pos];
	if(iorv == Kernel::IV_I)
		return _families[_materials[mt]._pt[pos]]._type[Kernel::POS_I];
	ERRORMSG("Incorrect petition in iorv2pt " << int(iorv) << " " << int(pos));
	return Kernel::UNDEFINED_TYPE;
}


unsigned ParameterManager::getStates(M_TYPE mt, P_TYPE pt) const
{
	return _particles[pt]._states[mt].size();
}

std::string ParameterManager::getStateName(M_TYPE mt, P_TYPE pt, unsigned state) const
{
	if(state==Kernel::UNDEFINED_STATE)
		return "UNDEFINED_STATE";
	return _particles[pt]._states[mt][state];
}

std::string ParameterManager::getParticleStateName(M_TYPE mt, P_TYPE pt, unsigned st) const
{
	if(st==Kernel::UNDEFINED_STATE)
		return "UNDEFINED_STATE";
	if(_particles[pt]._states[mt].size() > 1)
		return getParticleName(mt, pt) + "_" + _particles[pt]._states[mt][st];
	if(_particles[pt]._states[mt].size() == 1 && !_particles[pt]._states[mt][0].empty())
		return getParticleName(mt, pt) + "_" + _particles[pt]._states[mt][st];
	return getParticleName(mt, pt);
}

unsigned ParameterManager::getStateNumber(M_TYPE mt, P_TYPE pt, const std::string &str) const
{
	for(unsigned i=0; i < _particles[pt]._states[mt].size(); ++i)
		if(_particles[pt]._states[mt][i] == str)
			return i;
	return Kernel::UNDEFINED_STATE;
}

bool ParameterManager::isStateDefined(M_TYPE mt, P_TYPE pt,  unsigned state) const
{
	if (state < _particles[pt]._states[mt].size())
		return true;
	else
		return false;
}

//given a string like HeV2, returns the result as a map
//accepts inputs in different formats for binary and non binary materials
//the output ID is correct (same number of pt and pos)
Kernel::ID ParameterManager::getID(M_TYPE mt, const std::string &txt) const
{
	Kernel::ID result;
	result._mt = mt;
	const bool binary = getMaterial(mt)._binary;
	{
		P_TYPE pt = getParticleNumber(mt, txt);
		if(pt != UNDEFINED_TYPE)
		{
			result._pt [getFamily(pt)] = 1;
			result._pos[getPPos(pt)] = 1;
			return result;
		}
	}
	vector<string> temps[2];
	unsigned t = 0;
	for(unsigned s=0; s < txt.size(); ++s)
	{
		char c = txt[s];
		if(c >= 'A' && c <= 'Z')
			temps[t].push_back(string()+c);
		else if(c >= '0' && c <= '9')
		{
			if(s == 0)
				return result;
			if(txt[s-1] >='0' && txt[s-1] <= '9') //another one...
				temps[t].back() += c;
			else
				temps[t].push_back(string()+c);
		}
		else if(c>='a' && c<='z' && temps[t].size())
			temps[t].back() += c;
		else if(c == '^')
		{
			if(binary == false || ++t > 1)
				return result;
		}
		else
			return result;
	}
	if(binary)
	{
		map<P_TYPE, unsigned>::iterator itPt = result._pt.end();
		for(vector<string>::iterator it=temps[0].begin(); it!=temps[0].end(); ++it)
		{
			if(isdigit((*it)[0]))
			{
				if(itPt == result._pt.end())
					return Kernel::ID();
				itPt->second = atoi(it->c_str());
			}
			else
			{
				P_TYPE pt = getParticleNumber(mt, *it); //literal
				if(pt == UNDEFINED_TYPE || result._pt.find(pt) != result._pt.end())
					return Kernel::ID();
				result._pt[pt] = 1;
				itPt = result._pt.find(pt);
			}
		}
	}
	else
	{
		unsigned nI = 0, nV = 0;
		string previousIt;
		for(vector<string>::iterator it=temps[0].begin(); it!=temps[0].end(); ++it)
		{
			if(isdigit((*it)[0]) == false)
			{
				if(*it == "I")
				{
					result._pos[Kernel::POS_I] = 1;
					nI = 1;
				}
				else if(*it == "V")
				{
					result._pt[Kernel::V_TYPE] = 1;
					nV = 1;
				}
				else
				{
					P_TYPE pt = getParticleNumber(mt, *it); //literal
					if(pt == UNDEFINED_TYPE || result._pt.find(pt) != result._pt.end() || pt == getMaterial(mt)._pt[0])
						return Kernel::ID();
					result._pt[getFamily(pt)] = 1;
				}
			}
			else
			{
				if(previousIt == "I")
				{
					map<P_POS, unsigned>::iterator itPos = result._pos.find(Kernel::POS_I);
					if(itPos == result._pos.end())
						return Kernel::ID();
					itPos->second = atoi(it->c_str());
					nI = itPos->second;
				}
				else
				{
					P_TYPE saved_pt = getParticleNumber(mt, previousIt);
					if(saved_pt == Kernel::UNDEFINED_TYPE)
						return Kernel::ID();
					saved_pt = getFamily(saved_pt);
					if(saved_pt == Kernel::UNDEFINED_TYPE)
						return Kernel::ID();
					map<P_TYPE, unsigned>::iterator itPt = result._pt.find(saved_pt);
					if(itPt == result._pt.end())
						return Kernel::ID();
					itPt->second = atoi(it->c_str());
					if(previousIt == "V")
						nV = itPt->second;
				}
			}
			previousIt = *it;
		}

		if(nI && nV) //IV model, add extra particles.
		{
			result._pt[getMaterial(mt)._pt[0]] = nI;
			result._pos[Kernel::POS_0] = nV;
		}

	}
	{  //fill positions
		map<P_POS, unsigned>::iterator itPos = result._pos.end();
		for(vector<string>::iterator it=temps[1].begin(); it!=temps[1].end(); ++it)
		{
			if(isdigit((*it)[0]))
			{
				if(itPos == result._pos.end())
					return Kernel::ID();
				itPos->second = atoi(it->c_str());
			}
			else
			{
				P_POS pos = getPositionNumber(mt, *it);
				if(pos == Kernel::NO_POS || result._pos.find(pos) != result._pos.end())
					return Kernel::ID();
				result._pos[pos] = 1;
				itPos = result._pos.find(pos);
			}
		}
	}
	//check that it is OK
	unsigned nPts = 0, nPos = 0;
	for(map<P_POS, unsigned>::iterator it = result._pos.begin(); it != result._pos.end(); ++it)
		nPos += it->second;
	for(map<P_TYPE, unsigned>::iterator it = result._pt.begin(); it != result._pt.end(); ++it)
		nPts += it->second;

	if(nPts > nPos) //more particles than positions, fill positions
	{   //check for position 1
		if(result._pt.find(getPositionAsFamily(mt, Kernel::POS_0)) != result._pt.end())
			result._pos[(binary ? Kernel::POS_1 : Kernel::POS_0)] = nPts - nPos;
		else
			result._pos[Kernel::POS_0] = nPts - nPos;
	}
	if(nPts < nPos) //more positions than particles, fill particles
	{
		P_TYPE pos2 = getPositionAsFamily(mt, (binary? Kernel::POS_1 : Kernel::POS_0));
		if(result._pos.find(Kernel::POS_0) != result._pos.end())
			result._pt[pos2] = nPos - nPts;
		else
			result._pt[getPositionAsFamily(mt, Kernel::POS_0)] = nPos - nPts;
	}
	//Checks
	for(map<P_TYPE, unsigned>::const_iterator it=result._pt.begin(); it!=result._pt.end(); ++it)
		if(it->first >= Domains::global()->PM()->getNFamilies())
			return Kernel::ID();
	for(map<P_POS, unsigned>::const_iterator it=result._pos.begin(); it!=result._pos.end(); ++it)
		if(it->first == Kernel::POS_V || it->first == Kernel::NO_POS)
			return Kernel::ID();

	return result;
}

Kernel::ID ParameterManager::getEpiID(M_TYPE mt, const std::string &txt) const
{
	Kernel::ID result;
	result._mt = mt;
	vector<string> temps;
	for(unsigned s=0; s < txt.size(); ++s)
	{
		char c = txt[s];
		if(c >= 'A' && c <= 'Z')
			temps.push_back(string()+c);
		else if(c >= '0' && c <= '9')
		{
			if(s == 0)
				return result;
			if(txt[s-1] >='0' && txt[s-1] <= '9') //another one...
				temps.back() += c;
			else
				temps.push_back(string()+c);
		}
		else if(c>='a' && c<='z' && temps.size())
			temps.back() += c;
		else
			return result;
	}
	map<P_TYPE, unsigned>::iterator itPt = result._pt.end();
	for(vector<string>::iterator it=temps.begin(); it!=temps.end(); ++it)
	{
		if(isdigit((*it)[0]))
		{
			if(itPt == result._pt.end())
				return Kernel::ID();
			itPt->second = atoi(it->c_str());
		}
		else
		{
			P_TYPE pt = getFamily(getParticleNumber(mt, *it)); //literal
			if(pt == UNDEFINED_TYPE || result._pt.find(pt) != result._pt.end())
				return Kernel::ID();
			result._pt[pt] = 1;
			itPt = result._pt.find(pt);
		}
	}
	return result;
}


//given a map returns the name. The output to print in << theMap
//it compacts the output when possible
string ParameterManager::getIDName(const Kernel::ID &theMap, bool printOne) const
{
	stringstream ss1, ss2;
	if(theMap._mt >= _materials.size())
		ERRORMSG("Wrong ID in getIDName");
	if(OKMC::ClusterParam::ID2pt(theMap) != Kernel::UNDEFINED_TYPE)
		return getParticleName(theMap._mt, OKMC::ClusterParam::ID2pt(theMap));
	if(_materials[theMap._mt]._binary == false)
	{
		for(map<P_TYPE,unsigned>::const_reverse_iterator it = theMap._pt.rbegin(); it != theMap._pt.rend(); ++it)
		{
			if(it->first != _materials[theMap._mt]._pt[0])
			{
				ss1 << getParticleName(theMap._mt, it->first);
				if(it->second != 1 || printOne)
					ss1 << it->second;
			}

		}
		for(map<P_POS, unsigned>::const_iterator it = theMap._pos.begin(); it != theMap._pos.end(); ++it)
		{
			if(it->first == Kernel::POS_I)
			{
				ss2 << getPositionName(theMap._mt, it->first);
				if(it->second != 1 || printOne)
					ss2 << it->second;
			}
		}
	}
	else  //binary materials
	{
		for(map<P_TYPE,unsigned>::const_reverse_iterator it = theMap._pt.rbegin(); it != theMap._pt.rend(); ++it)
		{
			ss1 << getParticleName(theMap._mt, it->first);
			if(it->second != 1 || printOne)
				ss1 << it->second;
		}
		for(map<P_POS, unsigned>::const_iterator it = theMap._pos.begin(); it != theMap._pos.end(); ++it)
		{
			ss2 << getPositionName(theMap._mt, it->first);
			if(it->second != 1 || printOne)
				ss2 << it->second;
		}
	}
	stringstream ret;
	if(ss1.str().size())
		ret << ss1.str();
	if(ss2.str().size() && _materials[theMap._mt]._binary)
		ret << "^";
	ret << ss2.str();
	return ret.str();
}

string ParameterManager::getEpiIDName(const Kernel::ID &theMap) const
{
	stringstream ss1;
	if(theMap._mt >= _materials.size())
		ERRORMSG("Wrong ID in getEpiIDName");
	if(OKMC::ClusterParam::ID2pt(theMap) != Kernel::UNDEFINED_TYPE)
		return getParticleName(theMap._mt, OKMC::ClusterParam::ID2pt(theMap));

	for(map<P_TYPE,unsigned>::const_reverse_iterator it = theMap._pt.rbegin(); it != theMap._pt.rend(); ++it)
	{
		ss1 << getParticleName(theMap._mt, it->first);
		if(it->second != 1)
			ss1 << it->second;
	}
	assert(theMap._pos.size() == 0);
	stringstream ret;
	if(ss1.str().size())
		ret << ss1.str();
	return ret.str();
}

string ParameterManager::getRawIDName(const Kernel::ID &theMap) const
{
	stringstream ret;
	if(theMap._mt >= _materials.size())
		ERRORMSG("Wrong ID in getRawIDName");
	for(map<P_TYPE,unsigned>::const_reverse_iterator it = theMap._pt.rbegin(); it != theMap._pt.rend(); ++it)
		ret << getParticleName(theMap._mt, it->first) << it->second;
	ret << '^';
	for(map<P_POS, unsigned>::const_iterator it = theMap._pos.begin(); it != theMap._pos.end(); ++it)
		ret << getRawPositionName(theMap._mt, it->first) << it->second;
	return ret.str();
}

string ParameterManager::getEDName(M_TYPE mt, unsigned edType) const
{
	return _edDefectName[mt][edType];
}

//-1 minds not found
int ParameterManager::getEDType(M_TYPE mt, const string &def) const
{
	for(unsigned i=0; i<_edDefectName[mt].size(); ++i)
		if(_edDefectName[mt][i] == def)
			return i;
	return -1;
}

string ParameterManager::getDefName(M_TYPE mt, unsigned ev, unsigned def) const
{
	if(ev == Event::MOBILEPARTICLE)
		return getParticleName(mt, def);
	if(ev == Event::CLUSTER)
		return getEDName(mt, def);
	return Event::getEName(Event::E_TYPE(ev));
}

bool ParameterManager::isGasLike(M_TYPE mt) const
{
	for(P_TYPE pt=0; pt < getNParticles(); ++pt)
		if(isParticleDefined(pt, mt))
			return false;
	if(_materials[mt]._bAmorphous || _materials[mt]._bEpitaxy || _materials[mt]._binary)
		return false;

	return true;
}

bool ParameterManager::isAmorphous(M_TYPE mt) const
{
	return _materials[mt]._bAmorphous;
}

double ParameterManager::getMigrationRate(const Kernel::MeshElement *pEle, P_TYPE pt, unsigned st) const
{
	M_TYPE mt = pEle->getMaterial();
	return _migrations[mt][pt][st](pEle).getRate(pEle->getDomain()->_pRM->getkT());
}

void ParameterManager::getBreakComponents(Kernel::M_TYPE mt, P_TYPE pt, P_POS pos, P_TYPE &emit1, P_TYPE &emit2) const
{
	const bool binary = getMaterial(mt)._binary;
	P_TYPE fam  = getFamily(pt);
	char iorv   = getIorV(mt, pt);
	emit1 = getParticle(mt, fam, P_POS(pos));
	P_POS otherPos = (pos == Kernel::POS_0 && binary)? Kernel::POS_1 : Kernel::POS_0;
	if(fam == Kernel::V_TYPE || iorv == Kernel::IV_NONE)
		emit2 = getAntisite(mt, P_POS(pos)); //antisite
	else
	{
		if(getPPos(pt) != Kernel::POS_V)
			emit2 = iorv2pt(mt, iorv, P_POS(pos));
		else
			emit2 = iorv2pt(mt, iorv, otherPos);
	}
}

}
