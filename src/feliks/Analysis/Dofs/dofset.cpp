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
 *  dofset.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 6/13/09.
 *  Copyright 2009 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "Analysis/Dofs/dofset.h"
#include "Analysis/Integrators/integrator.h"
#include <cassert>
#include <iostream>
#include <iomanip>


using namespace blue;


int	dofset ::	getNConstrainedDofs() const
{
	int count(0);
	for (int a=0; a<ndof(); a++) count += isConstrained(a) ? 1 : 0;
	return count;
}



void dofset :: localizeSolutionIncrement(const longvector& bigdelta)
{
	double* Delta = bigdelta.data;
	extern int global_iteration;
	
	for (int j=0; j < ndof() ; j++)
	{
		int gdof = ID(j);

		// transfer the result from the solution of the global system of equations
		// to the variable deltaU stored at each node.
		if (gdof >= 0)  delta(j) = Delta[gdof];
		
		// remove the imposed bc after the first iteration
		else if (global_iteration > 0) delta(j) = 0.0;
	}
}



bool dofset :: setConstrainedDof(const int local)
{
	return ID.constrain(local);
}



directorDofset :: directorDofset()
{
	ID.setSize(2);
	dTheta.setZero(); Theta.setZero();
	
	Omega.setZero(); Omegan.setZero(); Omegan1.setZero();
	Alpha.setZero(); Alphan.setZero(); Alphan1.setZero();

	M.setZero();
}



directorDofset :: directorDofset(const ivector& d)
{
	ID.setSize(2);

	irotation r;
	r.beRotationWithoutDrill(d);
	rot.extractFromRotation(r);
	rot0 = rotn = rotn1 = rot;

	dTheta.setZero(); Theta.setZero();
	
	Omega.setZero(); Omegan.setZero(); Omegan1.setZero();
	Alpha.setZero(); Alphan.setZero(); Alphan1.setZero();
	
	M.setZero();
}



vector2& directorDofset :: bodyAcceleration(const evaluation_time when)
{
	switch ( when )
	{
		case ( dofset::tn) :
			return Alphan;
			break;
			
		case ( dofset::tn1) :
			return Alphan1;
			break;
			
		default:
			return Alpha;
	}	
}


vector2& directorDofset :: bodyVelocity(const evaluation_time when)
{
	switch ( when )
	{
		case ( dofset::tn) :
			return Omegan;
			break;
			
		case ( dofset::tn1) :
			return Omegan1;
			break;
			
		default:
			return Omega;
	}	
}




double&	directorDofset :: delta(const int l)
{
	return dTheta[l];
}



ivector directorDofset :: director(const evaluation_time when)
{
	iquaternion& q = rotation(when);
	return q*ivector(0.0, 0.0, 1.0);
}



void directorDofset :: printDOFS(std::ostream& of) const
{
}



iquaternion& directorDofset :: rotation(const evaluation_time when)
{
	switch ( when )
	{
		case (dofset::t0) :
			return rot0;
			break;
	
        case (dofset::tn) :
			return rotn;
			break;
			
		case (dofset::tn1) :
			return rotn1;
			break;
			
		default:
			return rot;
	}	
}



void directorDofset :: recoverFromLastSolution()
{
	dTheta.setZero();
	Theta.setZero();
	rot   = rotn1   = rotn;
	Omega = Omegan1 = Omegan;
	Alpha = Alphan1 = Alphan;
}



ivector directorDofset :: spatialAcceleration(const evaluation_time when)
{
	vector2&      a = bodyAcceleration(when);
	iquaternion&  q = rotation(when);
	return q*ivector(a(0), a(1), 0.0);
}



ivector directorDofset :: spatialVelocity(const evaluation_time when)
{
	vector2&      a = bodyVelocity(when);
	iquaternion&  q = rotation(when);
	return q*ivector(a(0), a(1), 0.0);
}



vector2& directorDofset ::  bodyStepRotationVector()
{
	return Theta;
}



rotationDofset :: rotationDofset():
	dtheta(0.0, 0.0, 0.0),
	theta(0.0, 0.0, 0.0),
	omega(0.0, 0.0, 0.0),
	omegan(0.0, 0.0, 0.0),
	omegan1(0.0, 0.0, 0.0),
	alpha(0.0, 0.0, 0.0),
	alphan(0.0, 0.0, 0.0),
	alphan1(0.0, 0.0, 0.0)
{
	ID.setSize(3);
}



const ivector rotationDofset :: bodyAcceleration(const evaluation_time when)
{
	switch ( when )
	{
		case ( dofset::tn) :
			return rotn.conjugate()*alphan;
			break;
			
		case ( dofset::tn1) :
			return rotn1.conjugate()*alphan1;
			break;
			
		default:
			return rot.conjugate()*alpha;
	}	
}



const ivector rotationDofset :: bodyVelocity(const evaluation_time when)
{
	switch ( when )
	{
		case ( dofset::tn) :
			return rotn.conjugate()*omegan;
			break;
			
		case ( dofset::tn1) :
			return rotn1.conjugate()*omegan1;
			break;
			
		default:
			return rot.conjugate()*omega;
	}	
}



double&	rotationDofset :: delta(const int l)
{
	return dtheta(l);
}



ivector rotationDofset :: director(const evaluation_time when)
{
	iquaternion& q = rotation(when);
	return q*ivector(0.0, 0.0, 1.0);
}



ivector& rotationDofset ::  rotationVector()
{
	return dtheta;
}



void rotationDofset :: printDOFS(std::ostream& of) const
{
}



void rotationDofset :: recoverFromLastSolution()
{
	dtheta.setZero();
	theta.setZero();
	
	rot   = rotn1   = rotn;
	omega = omegan1 = omegan;
	alpha = alphan1 = alphan;
}



iquaternion& 	rotationDofset :: rotation(const evaluation_time when)
{
	switch ( when )
	{
		case (dofset::t0) :
			return rot0;
			break;

		case (dofset::tn) :
			return rotn;
			break;
			
		case ( dofset::tn1) :
			return rotn1;
			break;
			
		default:
			return rot;
	}	
}



const ivector rotationDofset ::  bodyStepRotationVector()
{
	return rot.conjugate()*theta;
}



ivector& rotationDofset ::  spatialStepRotationVector()
{
	return theta;
}



ivector& rotationDofset :: spatialVelocity(const evaluation_time when)
{
	switch ( when )
	{
		case ( dofset::tn) :
			return omegan;
			break;
			
		case ( dofset::tn1) :
			return omegan1;
			break;
			
		default:
			return omega;
	}	
}



ivector& rotationDofset :: spatialAcceleration(const evaluation_time when)
{
	switch ( when )
	{
		case ( dofset::tn) :
			return alphan;
			break;
			
		case ( dofset::tn1) :
			return alphan1;
			break;
			
		default:
			return alpha;
	}	
}



scalarDofset :: scalarDofset() :
	dp(0.0),
	p(0.0), pn(0.0), pn1(0.0),
    V(0.0), Vn(0.0), Vn1(0.0),
    A(0.0), An(0.0), An1(0.0),
	F(0.0)
{
	ID.setSize(1);
}



double& scalarDofset ::  acceleration(const evaluation_time when)
{
	switch ( when )
	{
		case ( dofset::tn) :
			return An;
			break;
			
		case ( dofset::tn1) :
			return An1;
			break;
			
		default:
			return A;
	}
}



double& scalarDofset ::  displacement(const evaluation_time when)
{
	switch ( when )
	{
		case (dofset::tn):
			return pn;
			break;
			
		case (dofset::tn1):
			return pn1;
			break;
			
		default:
			return p;
	}
}



double& scalarDofset ::  increment()
{
	return dp;
}



double& scalarDofset ::  delta(const int l)
{
	return dp;
}



void scalarDofset :: printDOFS(std::ostream& of) const
{
	extern double global_macheps;
	
	double dof = pn1;
	dof = fabs(dof) > 1000.0*global_macheps ? dof : 0.0;
	of << std::scientific << std::setw((int)of.precision()+8) << dof;
}



void scalarDofset :: recoverFromLastSolution()
{
	dp = 0.0;
	p  = pn1 = pn;
	V  = Vn1 = Vn;
	A  = 0.0; //An1 = An;
}



void scalarDofset :: setDof(const double v)
{
	p = v;
}



double& scalarDofset ::  velocity(const evaluation_time when)
{
	switch ( when )
	{
		case (dofset::tn):
			return Vn;
			break;
		
		case (dofset::tn1):
			return Vn1;
			break;
			
		default:
			return V;
	}
}



vectorDofset :: vectorDofset() 
{
	ID.setSize(3);
	
	dU.setZero();
	
	U.setZero();
	Un.setZero();
	Un1.setZero();
	
	V.setZero();
	Vn.setZero();
	Vn1.setZero();

	A.setZero();
	An.setZero();
	An1.setZero();
	
	F.setZero();
}



ivector& vectorDofset ::  acceleration(const evaluation_time when)
{
	switch ( when )
	{
		case (dofset::tn):
			return An;
			break;
			
		case (dofset::tn1):
			return An1;
			break;
			
		default:
			return A;
	}
}



ivector& vectorDofset ::  displacement(const evaluation_time when)
{
	switch ( when )
	{
		case (dofset::tn):
			return Un;
			break;
			
		case (dofset::tn1):
			return Un1;
			break;
			
		default:
			return U;
	}
}



double& vectorDofset ::  delta(const int l)
{
	return dU(l);
}



void vectorDofset :: printDOFS(std::ostream& of) const
{
	extern double global_macheps;	
	
	for (int k=0; k<3; k++)
	{
		double dof = Un1(k);
		dof = fabs(dof) > 1000.0*global_macheps ? dof : 0.0;
		of << " " << scientific << std::showpos << right << setw((int)of.precision()+8) << dof;
	}
}



void vectorDofset :: recoverFromLastSolution()
{
	dU.setZero();
	U = Un1 = Un;
	V = Vn1 = Vn;
	A = An1 = An;
}



void vectorDofset :: setDof(const int l, const double v)
{
	assert(l<ndof());
	U(l) = v;
}



ivector& vectorDofset ::  velocity(const evaluation_time when)
{
	if ( when == dofset::tna)
		return V;
	
	else if ( when == dofset::tn)
		return Vn;

	else 
		return Vn1;

	return V;
}



const ivector& vectorDofset ::  velocity(const evaluation_time when) const
{
	return const_cast<vectorDofset&>(*this).velocity(when);
}



std::ostream&  operator<<(std::ostream &os, const directorDofset &ds)
{
	os << ds.ID;
	// Solutions
	os << "\n   Iter. rotation vector: " << ds.dTheta[0] << " " << ds.dTheta[1];
	os << "\n   Rotation at tn       : " << ds.rotn;
	os << "\n   Velocity at tn+1     : " << ds.Omegan1;
	os << "\n   Velocity at tn       : " << ds.Omegan;
	os << "\n   Acceleration at tn+1 : " << ds.Alphan1;
	os << "\n   Acceleration at tn   : " << ds.Alphan;
	
	return os;
}



std::ostream&  operator<<(std::ostream &os, const vectorDofset &ds)
{
	os << ds.ID;
	
	// Solutions
	os << "\n   Solution at tn+1     : " << ds.Un1;
	os << "\n   Solution at tn       : " << ds.Un;
	os << "\n   Last solver solution : " << ds.dU;
	os << "\n   Velocity at tn+1     : " << ds.Vn1;
	os << "\n   Velocity at tn       : " << ds.Vn;
	os << "\n   Acceleration at tn+1 : " << ds.An1;
	os << "\n   Acceleration at tn   : " << ds.An;
	
	return os;
}



std::ostream&  operator<<(std::ostream &os, const scalarDofset &ds)
{
	os << ds.ID;
	
	// Solutions
	os << "\n   Solution at tn+1     : " << ds.pn1;
	os << "\n   Solution at tn       : " << ds.pn;
	os << "\n   Last solver solution : " << ds.dp;
	os << "\n   Velocity at tn+1     : " << ds.Vn1;
	os << "\n   Velocity at tn       : " << ds.Vn;
	os << "\n   Acceleration at tn+1 : " << ds.An1;
	os << "\n   Acceleration at tn   : " << ds.An;
	
	return os;
}



std::ostream&  operator<<(std::ostream &os, const rotationDofset &ds)
{
	os << ds.ID;
	
	// Solutions
	os << "\n   Iter. rotation vector: " << ds.dtheta;
	os << "\n   Step rotation vector : " << ds.theta;
	os << "\n   Rotation at tn       : " << ds.rotn;
	os << "\n   Velocity at tn+1     : " << ds.omegan1;
	os << "\n   Velocity at tn       : " << ds.omegan;
	os << "\n   Acceleration at tn+1 : " << ds.alphan1;
	os << "\n   Acceleration at tn   : " << ds.alphan;
	
	return os;
}




