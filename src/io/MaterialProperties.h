/*
 * MaterialProperties.h
 *
 *  Created on: May 7, 2013
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

#ifndef MATERIALPROPERTIES_H_
#define MATERIALPROPERTIES_H_

#include "kernel/ParticleType.h"
namespace IO
{

struct MaterialProperties
{
	MaterialProperties(const std::string &n, const std::string &s) : _name(n), _shortName(s),
			amorphousExpansion_(0)
			{ _pt[0] = _pt[1] = Kernel::UNDEFINED_TYPE; }

	std::string    _name, _shortName;
	unsigned       _element;
	bool           _bAmorphous;
	bool		   _bEpitaxy;
	bool           _selfDiffusion;
	Kernel::M_TYPE _basicMaterial;
	Kernel::M_TYPE _amorphMaterial;
	Kernel::P_TYPE   _alloy;
	std::string    _alloyName;
	Kernel::P_TYPE   _pt[2];
	std::string    _s_pt[2];
	bool           _binary;
	float          _densityCm3;
	float		   _amorphThreshold;
	float          _densityAlloyCm3;
	double         _permittivity; // unitless
	double         _molarVolume;  // cm3/mol
	double         _Eg0;          // Bandgap @ 300K

	//mechanics
	double         E_; //Young's modulus
	double         nu_; //Poisson's ratio
	double         thermalExp_; //thermal expansion
	double         tempRef_; //reference temperature
	float		   amorphousExpansion_;
	std::map<std::string, std::pair<double, double> > eigenstrains_;
	double         alloyE_;
	double         alloynu_;
	double         alloythermalExp_;
	double         alloytempRef_;
};

}

#endif /* MATERIALPROPERTIES_H_ */

