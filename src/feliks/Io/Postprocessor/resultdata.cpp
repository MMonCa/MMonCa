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
/*
 *  resultdata.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 7/21/08.
 *
 */

#include "Io/Postprocessor/resultdata.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <typeinfo>

using namespace blue;


static int sizeofDictionary = 19;

static std::pair<resultCode,std::string> resultDictionary[]={
	std::pair<resultCode,std::string>(RESULT_PLASTIC_SLIP,	"plastic_slip"), 
	std::pair<resultCode,std::string>(RESULT_HEAT_FLUX,		"heat_flux"),
	std::pair<resultCode,std::string>(RESULT_STRAIN,		"strain"),
	std::pair<resultCode,std::string>(RESULT_STRESS,		"stress"),
	std::pair<resultCode,std::string>(RESULT_GENERATION,	"generation"), //5
	std::pair<resultCode,std::string>(RESULT_AXIAL_FORCE,	"axial_force"),
	std::pair<resultCode,std::string>(RESULT_AXIAL_STRAIN,	"axial_strain"),
	std::pair<resultCode,std::string>(RESULT_STRAIN_GAMMA,	"rod_strain_gamma"),
	std::pair<resultCode,std::string>(RESULT_STRAIN_OMEGA,	"rod_strain_omega"),
	std::pair<resultCode,std::string>(RESULT_STRESS_N,		"rod_stress_N"), //10
	std::pair<resultCode,std::string>(RESULT_STRESS_M,		"rod_stress_M"),
	std::pair<resultCode,std::string>(RESULT_TEMPERATURE,	"temperature"),
	std::pair<resultCode,std::string>(RESULT_ENTROPY,		"entropy"),
	std::pair<resultCode,std::string>(RESULT_ENERGY,		"energy"),
	std::pair<resultCode,std::string>(RESULT_PRESSURE,		"pressure"), //15
	std::pair<resultCode,std::string>(RESULT_CELL_NUMBER,	"cell_number"), 
	std::pair<resultCode,std::string>(RESULT_AGE,			"age"),
	std::pair<resultCode,std::string>(RESULT_PHASE,			"phase_field"),
	std::pair<resultCode,std::string>(RESULT_DAMAGE,		"damage")
};


const static resultType resultDictionaryType[]=
{
	RESULT_TYPE_SCALAR,
	RESULT_TYPE_VECTOR,
	RESULT_TYPE_TENSOR,
	RESULT_TYPE_TENSOR,
	RESULT_TYPE_SCALAR,		// generation
	RESULT_TYPE_SCALAR,		// axial force
	RESULT_TYPE_SCALAR,		// axial strain
	RESULT_TYPE_VECTOR,	    // strain gamma
	RESULT_TYPE_VECTOR,		// strain omega
	RESULT_TYPE_VECTOR,		// stress N
	RESULT_TYPE_VECTOR,		// stress M
	RESULT_TYPE_SCALAR,		// temperature
	RESULT_TYPE_SCALAR,		// entropy
	RESULT_TYPE_SCALAR,		// energy
	RESULT_TYPE_SCALAR,		// pressure
	RESULT_TYPE_SCALAR,		// cell_number
	RESULT_TYPE_SCALAR,		// age
	RESULT_TYPE_SCALAR,		// phase_field
	RESULT_TYPE_SCALAR		// damage
};


const static resultDomain resultDictionaryDomain[]={
	RESULT_DOMAIN_NODE,
	RESULT_DOMAIN_NODE,
	RESULT_DOMAIN_NODE,
	RESULT_DOMAIN_NODE,
	RESULT_DOMAIN_BODY,		// generation
	RESULT_DOMAIN_NODE,		// axial force
	RESULT_DOMAIN_NODE,		// axial strain
	RESULT_DOMAIN_NODE,	    // strain gamma
	RESULT_DOMAIN_NODE,		// strain omega
	RESULT_DOMAIN_NODE,		// stress N
	RESULT_DOMAIN_NODE,		// stress M
	RESULT_DOMAIN_NODE,		// temperature
	RESULT_DOMAIN_NODE,		// entropy
	RESULT_DOMAIN_NODE,		// energy
	RESULT_DOMAIN_NODE,		// pressure
	RESULT_DOMAIN_BODY,		// cell_number
	RESULT_DOMAIN_BODY,		// age
	RESULT_DOMAIN_NODE,		// phase_field
    RESULT_DOMAIN_ELEMENT   // damage
};



resultdata :: resultdata() :
	name(),
	code(RESULT_NONE),
	type(RESULT_TYPE_NONE),
	domain(RESULT_DOMAIN_NONE)
{}

resultdata :: resultdata(const resultCode theCode) :
	name(),
	code(RESULT_NONE),
	type(RESULT_TYPE_NONE),
	domain(RESULT_DOMAIN_NONE)
{
	int i(0);
	while (i< sizeofDictionary )
	{
		if (resultDictionary[i].first == theCode) 
		{
			name = resultDictionary[i].second;
			code = resultDictionary[i].first;
			type = resultDictionaryType[i];
			domain=resultDictionaryDomain[i];
			break;
		}
		i++;
	}
}



resultdata :: resultdata(const std::string& rname) :
	name(),
	code(RESULT_NONE),
	type(RESULT_TYPE_NONE),
	domain(RESULT_DOMAIN_NONE)
{
	int i(0);
	while (i< sizeofDictionary )
	{
		if (resultDictionary[i].second == rname) 
		{
			name = resultDictionary[i].second;
			code = resultDictionary[i].first;
			type = resultDictionaryType[i];
			domain=resultDictionaryDomain[i];
			break;
		}
		i++;
	}
}



resultdata :: resultdata(const resultdata& r) :
	name(r.name),
	code(r.code),
	type(r.type),
	data(r.data),
	domain(r.domain)
{
}




int resultdata :: getNComponents() const
{
	int ret(0);
	
	if      (type == RESULT_TYPE_SCALAR) ret = 1;
	else if (type == RESULT_TYPE_VECTOR) ret = 4;
	else if (type == RESULT_TYPE_TENSOR) ret = 9;
	return ret;
}




resultCode resultdata :: getResultCode(const std::string& name)
{
	resultCode r(RESULT_NONE);
	
	int i(0);
	while (i< sizeofDictionary )
	{
		if (resultDictionary[i].second == name) 
		{
			r = resultDictionary[i].first;
			break;
		}
		i++;
	}
	return r;
}




resultDomain resultdata :: getResultDomain(const std::string& name)
{
	resultDomain r(RESULT_DOMAIN_NONE);
	
	int i(0);
	while (i< sizeofDictionary )
	{
		if (resultDictionary[i].second == name) 
		{
			r = resultDictionaryDomain[i];
			break;
		}
		i++;
	}
	return r;
}




std::string resultdata :: getResultName(const resultCode theCode)
{
	std::string name;
	
	int i(0);
	while (i< sizeofDictionary )
	{
		if (resultDictionary[i].first == theCode) 
		{
			name = resultDictionary[i].second;
			break;
		}
		i++;
	}
	return name;
}




int resultdata 	::	getResultNComponents(const resultCode theCode)
{
	int ret(0);
	
	switch ( getResultType(theCode) )
	{
		case RESULT_TYPE_SCALAR: ret = 1; break;
		case RESULT_TYPE_VECTOR: ret = 3; break;
		case RESULT_TYPE_TENSOR: ret = 9; break;
        default: ret = 0;
	}
	return ret;
}




resultType resultdata :: getResultType(const std::string& name)
{
	resultType  t(RESULT_TYPE_NONE);
	int i(0);
	while (i< sizeofDictionary )
	{
		if (resultDictionary[i].second == name) 
		{
			t = resultDictionaryType[i];
			break;
		}
		i++;
	}
	return t;
}




resultType resultdata :: getResultType(const resultCode theCode)
{
	resultType  t(RESULT_TYPE_NONE);
	
	int i(0);
	while (i< sizeofDictionary )
	{
		if (resultDictionary[i].first == theCode) 
		{
			t = resultDictionaryType[i];
			break;
		}
		i++;
	}
	return t;
}




void resultdata :: setZero()
{
	if			(type == RESULT_TYPE_SCALAR) data(0,0) = 0.0;
	else if		(type == RESULT_TYPE_VECTOR) data(0,0) = data(0,1) = data(0,2) = 0.0;
	else if		(type == RESULT_TYPE_TENSOR) data.setZero();
}




double& resultdata :: operator[](const size_t i)
{
	return data((int)i/3,i%3);
}




resultdata&  resultdata ::  operator=(const double &t)
{
	if (this->type != RESULT_TYPE_SCALAR) throw std::runtime_error("Setting a scalar to a non-scalar");
	data(0,0) = t;
	return *this;
}




resultdata&  resultdata ::  operator=(const ivector &t)
{
	if (this->type != RESULT_TYPE_VECTOR) throw std::runtime_error("Adding a vector to a non-vector");
	data(0,0) = t(0);
	data(0,1) = t(1);
	data(0,2) = t(2);
	return *this;
}




resultdata&  resultdata ::  operator=(const itensor &t)
{
	if (this->type != RESULT_TYPE_TENSOR) throw std::runtime_error("Adding a tensor to a non-tensor");
	data = t;
	return *this;
}




resultdata&  resultdata ::  operator+=(const double &t)
{
	if (this->type != RESULT_TYPE_SCALAR) throw std::runtime_error("Adding a scalar to a non-scalar");
	data(0,0) += t;
	return *this;
}




resultdata&  resultdata ::  operator+=(const ivector &t)
{
	if (this->type != RESULT_TYPE_VECTOR) throw std::runtime_error("Adding a vector to a non-vector");
	data(0,0) += t(0);
	data(0,1) += t(1);
	data(0,2) += t(2);
	return *this;
}




resultdata&  resultdata ::  operator+=(const itensor &t)
{
	if (this->type != RESULT_TYPE_TENSOR) throw std::runtime_error("Adding a tensor to a non-tensor");
	data += t;
	return *this;
}




resultdata&  resultdata :: operator+=(const resultdata &t)
{
	if (this-> type != t.type)
		throw std::runtime_error("Non-matching types in result sum.");
	
	if		(this->type == RESULT_TYPE_SCALAR) *this += t.data(0,0);
	else if	(this->type == RESULT_TYPE_VECTOR) *this += ivector(t.data(0,0), t.data(0,1), t.data(0,2) );
	else if	(this->type == RESULT_TYPE_TENSOR) *this += t.data;
	
	return *this;
}




resultdata&  resultdata :: operator*=(const double  &f)
{
	if (this->type == RESULT_TYPE_SCALAR) data(0,0) *= f;
	else
		data *= f;
	return *this;
}




resultdata   operator*(const resultdata &t , const double   a)
{
	resultdata ret(t);
	ret *= a;
	return ret;
}




resultdata&  resultdata :: operator=(const resultdata   &r)
{
	name = r.name ; 
	code = r.code ; 
	type = r.type ; 
	domain = r.domain ; 
	return *this;
}




std::ostream& operator<<(std::ostream &os, const resultdata &t)
{
	if       (t.type == RESULT_TYPE_SCALAR) os << t.data(0,0);
	else if  (t.type == RESULT_TYPE_VECTOR) 
		os << std::setw(9) << std::fixed << std::showpoint << std::right << t.data(0,0) << " " << t.data(0,1) << " " << t.data(0,2);
	else if  (t.type == RESULT_TYPE_TENSOR)
	{
		os << std::showpos << std::setprecision(5) << std::scientific << std::showpoint << std::right;
		for (int i=0; i<3; i++) for (int j=0; j<3; j++) os << t.data(i,j) << " ";
	}
	return os;
}
