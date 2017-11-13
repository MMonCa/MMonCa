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
 *  energies.h
 *  feliks
 *
 *  Created by Ignacio Romero on Fri Oct 24 2003.
 *  Copyright (c) 2003 Universidad Politecnica de Madrid. All rights reserved.
 *
 *
 *  a energies is just an array of "doubles" each of which has a specific name,
 *  such as "energy "angularMomentum", "errorEstimate", "user". In particular
 *  the "user" array is left so that users can program elements that employ
 *  their own type of "energy".
 *  Since this type of vectors are used both for the global analysis and for each element
 *  it is handy to have a type defined just for this purpose.
 *
 *  If energies are dumped to a file, these are the corresponding values:
 *      (1)  the time
 *      (2)  Kinetic energy
 *      (3)  Potential energy
 *      (4)  Total energy
 *      (5)  Internal
 *      (6)  Free energy
 *      (7)  Entropy
 *      (8)  Angular momentum 1
 *      (9)  Angular momentum 2
 *      (10) Angular momentum 3
 *      (11) Linear momentum 1
 *      (12) Linear momentum 2
 *      (13) Linear momentum 3
 *      (14) Main error estimate (sum of element contributions)
 *      (15) Aux. error estimate 1 (sum)
 *      (16) Aux. error estimate 2 (sum)
 *      (17) Aux. error estimate 3 (sum)
 *      (18) Aux. error estimate 4 (sum)
 *      (19) Time integration total error
 *      (20) Time integration error in last time step
 *      (21) User integral 1
 *      (22) User integral 2 ...
 *      
 */

#pragma once
#ifndef _feliks_energies_h
#define _feliks_energies_h

#include "Math/tensor.h"
#include <cstdio>
#include <iostream>
#include <ostream>


class energies
{
    	
public:
	double energy;
    double kinetic;
    double potential;
	double entropy;
	double helmholtz;
	double internale;
	double mass;
	double volume;
    double dissipation;

    blue::ivector angularMomentum;
    blue::ivector linearMomentum;
    double errorEstimate;
    double error[4];
    double integrationError;
    double stepError;
    int    nuserEnergies;
    double user[10];
	
	energies();
	
	void            add(const energies &e2);  // e1 += e2
    void            chop(const double tol);
	void            initialize(int nuserener=0);
	void            info();
	void            print(std::ostream& of=std::cout) const;
	void            printNames(std::ostream& of=std::cout) const;
	void            setZero();
	
	inline double   getStepError()			{ return stepError;}
	inline void     increaseIntegrationError(double e){ integrationError += e;}
	

	energies& operator+=(const energies &rhs);
};
	

#endif
