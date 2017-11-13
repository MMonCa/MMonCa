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
 *  energies.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on Fri Oct 24 2003.
 *  Copyright (c) 2003 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include "Analysis/energies.h"


#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "Io/message.h"

using namespace std;

double energy;
double kinetic;
double potential;
double entropy;
double helmholtz;
double internale;
double mass;
double volume;
double dissipation;



energies :: energies() :
	energy(0.0),
	kinetic(0.0),
	potential(0.0),
    entropy(0.0),
    helmholtz(0.0),
    internale(0.0),
	mass(0.0),
	volume(0.0),
    dissipation(0.0),
	angularMomentum(),
	linearMomentum()
{
    angularMomentum.setZero();
    linearMomentum.setZero();
}



/*  this function is usually called to accumulate the energy values of e2
    to e1, which is usually a global measure of the energy (and other quantities)
*/
void energies :: add(const energies &e2)
{
    // these are element contributions and need to be accumulated
    energy             += e2.energy;
    kinetic            += e2.kinetic;
    potential          += e2.potential;
	internale		   += e2.internale;
	helmholtz		   += e2.helmholtz;
	entropy			   += e2.entropy;
	mass			   += e2.mass;
    volume	           += e2.volume;
    dissipation        += e2.dissipation;
		
    angularMomentum    += e2.angularMomentum;
    linearMomentum     += e2.linearMomentum;
    errorEstimate      += e2.errorEstimate;
	
    // extra errors
    //for (k=0; k<4; k++)    error[k] += e2.error[k];

    // user defined quantities
    //for (k=0; k<10; k++) user[k]  += e2.user[k];
}




void energies :: chop(const double tol)
{
    if (fabs(energy)    < tol) energy = 0.0;
    if (fabs(kinetic)   < tol) kinetic = 0.0;
    if (fabs(potential) < tol) potential = 0.0;
    if (fabs(internale) < tol) internale = 0.0;
    if (fabs(helmholtz) < tol) helmholtz = 0.0;
    if (fabs(entropy)   < tol) entropy = 0.0;
    if (fabs(mass)      < tol) mass = 0.0;
    if (fabs(volume)    < tol) volume = 0.0;
    if (fabs(dissipation)< tol) dissipation = 0.0;
    
    if (fabs(angularMomentum(0)) < tol) angularMomentum(0) = 0.0;
    if (fabs(angularMomentum(1)) < tol) angularMomentum(1) = 0.0;
    if (fabs(angularMomentum(2)) < tol) angularMomentum(2) = 0.0;
    if (fabs(linearMomentum(0) ) < tol) linearMomentum(0)  = 0.0;
    if (fabs(linearMomentum(1) ) < tol) linearMomentum(1)  = 0.0;
    if (fabs(linearMomentum(2) ) < tol) linearMomentum(2)  = 0.0;
}



void energies ::  initialize(int n)
{
	setZero();
	errorEstimate    = 0.0;
    integrationError = 0.0;
    nuserEnergies    = n;
}




void energies :: info()
{
    Message("\n\n");
    Message("\n Kinetic energy    = % 10.5e", kinetic);
    Message("\n Potential energy  = % 10.5e", potential);
    Message("\n Total energy      = % 10.5e", energy);
    Message("\n Error estimate    = % 10.5e", errorEstimate);
    Message("\n Integration error = % 10.5e", integrationError);
    Message("\n Step error        = % 10.5e", stepError);
    Message("\n Dissipation       = % 10.5e", dissipation);
}




energies& energies :: operator+=(const energies &rhs)
{
	this->add(rhs);
	return *this;
}



/*  this function is usually called when the time counter advances, so some
    manipulations can be done at this stage.
    This goes in a line of the energy output file, whatever its name (user defined or default).
    The first entry in this line is the time, and must be dumped before calling PrintEnergiesInFile.
*/
void energies :: print(std::ostream &of) const
{    
	of << setw(12) << scientific << kinetic;
	of << setw(12) << scientific << potential;
	of << setw(12) << scientific << energy;
	of << setw(12) << scientific << internale;
	of << setw(12) << scientific << helmholtz;
	of << setw(12) << scientific << entropy;
    of << setw(12) << scientific << dissipation;
	of << setw(12) << scientific << linearMomentum[0];
	of << setw(12) << scientific << linearMomentum[1];
	of << setw(12) << scientific << linearMomentum[2];
	of << setw(12) << scientific << angularMomentum[0];
	of << setw(12) << scientific << angularMomentum[1];
	of << setw(12) << scientific << angularMomentum[2];
	
	/*
    // the error quantities 11, 12, 13, 14, 15
    fprintf(fp, "% 10.5e % 10.5e % 10.5e % 10.5e % 10.5e",
            errorEstimate, error[0], error[1], error[2], error[3]);
    
    // the time integration error
    fprintf(fp, " % 10.5e % 10.5e", integrationError, stepError);

    // the user integrals
    for(int k=0; k<nuserEnergies; k++) fprintf(fp,"% 10.5e", user[k]);
	 */
}




void energies :: printNames(std::ostream &of) const
{    
	of << setw(12) << "kinetic";
	of << setw(12) << "potential";
	of << setw(12) << "total nrg";
	of << setw(12) << "internal";
	of << setw(12) << "helmholtz";
	of << setw(12) << "entropy";
	of << setw(12) << "dissip.";
	of << setw(36) << "l i n e a r   m o m e n t u m    ";
	of << setw(36) << "a n g u l a r   m o m e n t u m  ";
	of << flush;
}



/* set all the quantities to zero. This should be called before integrating over
   all the elements. Note that the quantities errorEstimate and integrationError
   are not set to zero
   because we are interested in accumulating the previous errors
*/
void energies :: setZero()
{
    int k;
    
    energy    = 0.0;
    kinetic   = 0.0;
    potential = 0.0;
	entropy   = 0.0;
	helmholtz = 0.0;
	internale = 0.0;
	mass      = 0.0;
	volume    = 0.0;
    dissipation = 0.0;
	
    angularMomentum.setZero();
	linearMomentum.setZero();
    
    for (k=0; k<4; k++)  error[k] = 0.0;
    for (k=0; k<10; k++) user[k]  = 0.0;
}

