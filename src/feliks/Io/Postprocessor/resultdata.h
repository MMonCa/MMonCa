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
 *  resultdata.h
 *  feliks
 *
 *
 *  Created by Ignacio Romero on 7/21/08.
 *
 *  when a finite element is asked for a "result" at a Gauss point or a node, it can be of (at least) three
 *  types: scalar, vector, or tensor. The resultdata class provides a convenient way of passing
 *  arguments to functions that compute such quantities.
 */


#ifndef _resultdata_h
#define _resultdata_h


#include "Math/tensor.h"
#include <string>
#include <utility>
#include <vector>


enum resultCode
{
	RESULT_NONE,
	RESULT_STRAIN,
	RESULT_STRESS,
	RESULT_HEAT_FLUX,
	RESULT_F_RESULTANTS,
	RESULT_M_RESULTANTS,
	RESULT_PLASTIC_SLIP,
	RESULT_GENERATION,
	RESULT_AXIAL_FORCE,
	RESULT_AXIAL_STRAIN,
	RESULT_STRAIN_GAMMA,
	RESULT_STRAIN_OMEGA,
	RESULT_STRESS_N,
	RESULT_STRESS_M,
	RESULT_TEMPERATURE,
	RESULT_ENTROPY,
	RESULT_ENERGY,
	RESULT_PRESSURE,
	RESULT_CELL_NUMBER,
	RESULT_AGE,
	RESULT_PHASE,
    RESULT_DAMAGE
};


enum resultType
{
	RESULT_TYPE_NONE,
	RESULT_TYPE_SCALAR,
	RESULT_TYPE_VECTOR,
	RESULT_TYPE_TENSOR
};



enum resultDomain
{
	RESULT_DOMAIN_NONE,
	RESULT_DOMAIN_NODE,
	RESULT_DOMAIN_ELEMENT,
	RESULT_DOMAIN_BODY
};




class resultdata
{
	
public:
	resultdata();
	resultdata(const resultCode c);
	resultdata(const std::string& rname);
	resultdata(const resultdata& r);
	
	void                setZero();
	
	
	// information about an instance
	int					getNComponents() const;
	resultCode			getCode() const;
	resultType			getType() const;
	resultDomain		getDomain() const;
	std::string			getName() const;
	
	
	// the only math operations are to add and set the data
	resultdata&         operator=(const double  &t);
	resultdata&         operator=(const blue::ivector &t);
	resultdata&         operator=(const blue::itensor &t);
	resultdata&         operator=(const resultdata &r);
	double&             operator[](const size_t i);

	resultdata&         operator+=(const double  &t);
	resultdata&         operator+=(const blue::ivector &t);
	resultdata&         operator+=(const blue::itensor &t);
	resultdata&         operator+=(const resultdata &t);
	resultdata&         operator*=(const double  &f);
	friend resultdata   operator*(const resultdata &t, const double a);

	
	// io
	friend std::ostream& operator<<(std::ostream &os, const resultdata &t);
	
	
	// global lookup
	static resultCode 	getResultCode(const std::string& name);
	static resultDomain getResultDomain(const std::string& name);
	static resultType	getResultType(const std::string& name);
	static resultType	getResultType(const resultCode theCode);
	static std::string	getResultName(const resultCode theCode);
	static int			getResultNComponents(const resultCode theCode);


private:	
	std::string	name;
	resultCode	code;
	resultType	type;
	resultDomain domain;
	blue::itensor		data;
};



// inlines
inline resultCode	resultdata :: getCode() const { return code;}
inline resultDomain resultdata :: getDomain() const {return domain;}
inline resultType	resultdata :: getType() const { return type;}

#endif

