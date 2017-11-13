/* Finite Element Method Module
 *
 * Author: ignacio.romero@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain,
 *      and       Technical University of Madrid (UPM), Madrid, Spain
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
/* material.h
 *
 * iro
 * november 1999,
 * converted to c++ june 2006, iro
 * completely redesigned and mp added november 2006, iro
 */


/*
 The material class now serves two purposes:
 
 1) The main tasks of a material is to act as a blueprint for the constructions of new
 material points. This is implemented in the fucntion createMaterialPoint.
 2) Each material acts as a memory space where the constants of each material are stored
 and are referenced by the material points. Since they don't change, it is not worth
 it to store a copy of each one within each material point.
 
 The material point is responsible for the constitutive response, for storing its internal variables, etc,
 but the user needs to create a run-time template that indicates the type of material. The class
 material serves this purpose. When scanning the input file, feliks creates a table of materials, one
 per user-defined material 1 for elastic with e=1, nu=0.2, 2 for elastic with e=2, nu=0.3, material 3
 for elastoplastic with e=1, nu=0.2, H=9, etc. Later, when the elements are created, each individual
 element creates as many material points as quadrature points it has. Each material point must
 be created according to the templated provided by the initially allocated table.
 */



#pragma once
#ifndef _material_h
#define _material_h


#include <map>
#include <iostream>
#include <string>


class commandLine;
class materialPoint;



enum propertyName
{ 
		PR_LAMBDA, PR_MU,	// Lame coefficients
		PR_YOUNG,			// Young's modulus
		PR_POISSON,			// Poisson's ratio
		PR_BULK,			// Bulk modulus
		PR_CP,				// P waves velocity
		PR_CS,				// S waves velocity
		PR_SMAX,			// maximum normal stress
		PR_GF,				// fracture energy
		PR_PLSTRESS_C,		// plane stress wave velocity
		PR_ALPHA,			// thermal conductivity
		PR_TREF,			// reference temperature
		PR_C0,				// thermal capacitance
		PR_CONDUCTIVITY,	
		PR_HEAT_SUPPLY,		// heat supply per unit volume
        PR_HISO,            // isotropic hardening
        PR_HKINE,           // kinematic hardening
        PR_NU               // dynamic viscosity
	};





class material
{
	
public:
                                material();
                                material(const commandLine &cl);
	virtual                     ~material();
	static std::map<int, material*>   allMaterials;


	// these are the function to work with the list of AllMaterials
	static void                 freeMaterials();
	static material&            getMaterial(const int label);
	static material&            getMaterial(const std::string& name);
	static size_t               getNMaterials();
	static void                 listMaterials(std::ostream &of=std::cout);
	static void                 scan(const commandLine &cl);	
	
    
	// operations on single material types
    virtual bool                check() const=0;
	virtual double              density() const=0;
	virtual void                print(std::ostream &of=std::cout) const=0;
    virtual void                setRandom()=0;

    
    // information
    const int&                  label() const       {return _label;}
	const std::string&          name() const        {return _name;}

	
    
private:
	
	static int                  largestMaterialLabel;
	static bool                 addToList(material *m);
	
	std::string                 _name;
    int                         _label;
	
    friend class                element;
    friend class                materialPoint;
};
  





/* 
 The materialPoint objects are things that embody the constitutive 
 response of a single point in a body. As such, they must know how to compute a stress
 from a strain, for example, but are also responsible for book keeping their own 
 internal variables, if they have any, and updating them correspondingly.
*/
class materialPoint
{	
    public:
                                materialPoint();
                                materialPoint(material &m);
    virtual                     ~materialPoint(){};	

    virtual void                commitCurrentState()=0;
    virtual void                resetCurrentState()=0;
    //                          updateCurrentState(some gradient)=0;

};




#endif
