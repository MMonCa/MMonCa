/*
 * ParameterManager.h
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

#ifndef PARAMETERMANAGER_H_
#define PARAMETERMANAGER_H_
#include <string>

#include "kernel/Material.h"
#include "kernel/Event.h"
#include "kernel/Coordinates.h"
#include "kernel/ParticleType.h"
#include "io/ArrheniusAlloys.h"
#include "MaterialProperties.h"
#include "Element.h"
#include <string>
#include <vector>
#include <map>

namespace Kernel { class SubDomain; class MeshElement; }
struct Tcl_Interp;
namespace IO {

class FileParameters;

class ParameterManager
{
public:
	ParameterManager(Tcl_Interp *, const IO::FileParameters *);
	~ParameterManager();

	const std::string &getMaterialName(Kernel::M_TYPE m) const { return _materials[m]._name; }
	const std::string &getAlloyName(Kernel::M_TYPE m) const { return _materials[m]._alloyName; }
	      std::string  getParticleName(Kernel::M_TYPE m, Kernel::P_TYPE d) const;
	const std::string &getElementName(Kernel::P_TYPE d) const { return _families.at(d)._s_family; }
	const std::string &getFamilyName(Kernel::P_TYPE i) const { return _families.at(i)._file; }
	std::string        getPositionName(Kernel::M_TYPE mt, Kernel::P_POS) const;
	std::string        getRawPositionName(Kernel::M_TYPE mt, Kernel::P_POS) const;
	Kernel::P_TYPE       getPositionAsFamily(Kernel::M_TYPE mt, Kernel::P_POS) const;

	bool  isParticleDefined(Kernel::P_TYPE pt, Kernel::M_TYPE mt) const { return pt == Kernel::UNDEFINED_TYPE ? false :_particles[pt]._materials[mt]; }
	bool  isGasLike(Kernel::M_TYPE mt) const;
	bool  isAmorphous(Kernel::M_TYPE mt) const;
	bool  isStateDefined(Kernel::M_TYPE, Kernel::P_TYPE, unsigned state) const;

	Kernel::Event::E_TYPE getDefectType(const std::string &s, Kernel::M_TYPE mt) const
	{ std::map<std::string, Kernel::Event::E_TYPE>::const_iterator it = _defects[mt].find(s);
	  return (it == _defects[mt].end()? Kernel::Event::UNDEFINED_EVENT : it->second); }

	Kernel::M_TYPE     getMaterialNumber(const std::string &) const;
	const MaterialProperties & getMaterial(Kernel::M_TYPE m) const { return _materials[m]; }

	Kernel::P_TYPE	   getParticleNumber(Kernel::M_TYPE mt, const std::string &) const;
	Kernel::P_POS        getPositionNumber(Kernel::M_TYPE mt, const std::string &) const;
	Kernel::M_TYPE     getNMaterials() const { return _materials.size(); }
	Kernel::P_TYPE	   getNParticles() const { return _particles.size(); }
	Kernel::P_TYPE       getNFamilies() const  { return _families.size(); }
	Kernel::P_TYPE       getFamily(Kernel::P_TYPE d) const { return _particles[d]._family; }
	Kernel::P_POS        getPPos(Kernel::P_TYPE d) const   { return _particles[d]._pos; }
	Kernel::P_TYPE       getParticle(Kernel::M_TYPE mt, Kernel::P_TYPE d, Kernel::P_POS pos) const;
	Kernel::P_TYPE       getParticleForCluster(Kernel::M_TYPE mt, Kernel::P_TYPE d, Kernel::P_POS pos) const;
	Kernel::P_TYPE       getAntisite(Kernel::M_TYPE mt, Kernel::P_POS pos) const; //Antisite Si is Si^C
	Kernel::P_TYPE       iorv2pt(Kernel::M_TYPE mt, unsigned char iorv, Kernel::P_POS) const;
	char               getIorV(Kernel::M_TYPE mt, Kernel::P_TYPE p) const; //I or V for SiI, AsV, etc... returns 0 or 1
	char               isIorV(Kernel::M_TYPE mt, Kernel::P_TYPE p) const; //only SiI, CI, V
	bool               isPureImpurity(Kernel::M_TYPE mt, Kernel::P_TYPE p) const;
	bool               isImpurity(Kernel::M_TYPE mt, Kernel::P_TYPE p) const;
	unsigned           getStates(Kernel::M_TYPE, Kernel::P_TYPE) const;
	std::string        getStateName(Kernel::M_TYPE, Kernel::P_TYPE, unsigned) const;
	unsigned           getStateNumber(Kernel::M_TYPE, Kernel::P_TYPE, const std::string &) const;
	std::string        getParticleStateName(Kernel::M_TYPE, Kernel::P_TYPE, unsigned) const;
	const Element &    getElement(unsigned e) const { return _elements[e];  }
	Kernel::P_TYPE     getElementNumber(const std::string &) const;
	static void        getTokens(const std::string &key, char delimiter, std::vector<std::string> &tokens);
	void			   getBreakComponents(Kernel::M_TYPE, Kernel::P_TYPE, Kernel::P_POS, Kernel::P_TYPE &emit1, Kernel::P_TYPE &emit2) const;

	Kernel::ID  getID(Kernel::M_TYPE mt, const std::string &) const;
	Kernel::ID  getEpiID(Kernel::M_TYPE mt, const std::string &) const;
	std::string getIDName(const Kernel::ID &, bool printOne = false) const;
	std::string getEpiIDName(const Kernel::ID &) const;
	std::string getRawIDName(const Kernel::ID &) const;
	std::string getEDName(Kernel::M_TYPE mt, unsigned edType) const;
	int         getEDType(Kernel::M_TYPE mt, const std::string &def) const;
	std::string getDefName(Kernel::M_TYPE mt, unsigned ev, unsigned def) const; //generic for "all" defects

	double getMigrationRate(const Kernel::MeshElement *, Kernel::P_TYPE pt, unsigned st) const;
	void   setMigrationRate(Kernel::M_TYPE mt, unsigned st, const Kernel::P_TYPE &pt, const IO::ArrheniusAlloys &arr)
		{ _migrations[mt][pt][st] = arr; }
	float getMaximumCaptureRadius(Kernel::M_TYPE mt) const  { return _maxCapRad[mt]; }
	void  setMaximumCaptureRadius(Kernel::M_TYPE mt, float r) { if(r > _maxCapRad[mt]) _maxCapRad[mt] = r; }

private:
	void readElements (const IO::FileParameters *);
	void readMaterials(Tcl_Interp *, const IO::FileParameters *);
	void readDefects  (Tcl_Interp *, const IO::FileParameters *);
	void setAlloy     (const IO::FileParameters *);
	struct ParticleProperties
	{
		ParticleProperties(Kernel::P_TYPE fam, Kernel::P_POS pos) : _family(fam), _pos(pos) { }
		Kernel::P_TYPE _family;
		Kernel::P_POS _pos;
		std::vector<bool> _materials;
		std::vector<std::vector<std::string> > _states; //one array per material
	};
	struct FamilyProperties
	{
		FamilyProperties(const std::string &f, const std::string &d);
		std::string _file, _s_family;
		Kernel::P_TYPE _type[5];  //1, 2, I, V, NO_POS
	};
	std::vector<MaterialProperties> _materials;
	std::vector<ParticleProperties> _particles;
	std::vector<FamilyProperties>   _families;
	std::vector<std::map<std::string, Kernel::Event::E_TYPE> > _defects;
	IO::ArrheniusAlloys _migrations[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_STATES];
	std::vector<std::string> _edDefectName[Kernel::MAX_MATERIALS];
	float _maxCapRad[Kernel::MAX_MATERIALS];
	Element _elements[Kernel::MAX_ELEMENTS];

	mutable std::map<std::string, Kernel::P_TYPE> _partNameCache[Kernel::MAX_MATERIALS]; //to speed up get ParticleNumber
};

}

#endif /* PARAMETERMANAGER_H_ */
