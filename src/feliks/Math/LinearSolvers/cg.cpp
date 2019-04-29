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
 *  cg.cpp
 *
 *  Preconditioned Conjugate Gradient linear solver
 *  feliks
 *
 *  Created by Ignacio Romero on Wed May 05 2004.
 *  Copyright (c) 2004 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include "Math/LinearSolvers/cg.h"

#include "Io/logger.h"
#include "Math/linearsoe.h"
#include "Math/vector.h"
#include "Math/Sparse/sparse.h"
#include "Math/Sparse/compressed.h"
#include "General/timer.h"
#include "Io/message.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <utility>
#include <limits>

#define RTOL 1.e-10
#define ATOL 1.e-14


cg :: cg(const commandLine &cl) :
	linearSolver(cl),
	preconditioner("diagonal"),
	droptol(1e-2),
	perm(0),
	invperm(0),
	isPermuted(false),
	prec(0)
{
	name = "PCG";
	storage = SCSC;
	
	for (int k=0; k < cl.size(); k++)
	{
		const usercommand &uc = cl[k];		
		if      (uc.keyword() == "preconditioner") preconditioner = uc.option();
		else if (uc.keyword() == "droptol")        droptol        = uc.value();
	}
}


cg :: cg() :
	preconditioner("diagonal"),
	droptol(1e-2),
	perm(0),
	invperm(0),
	isPermuted(false),
	prec(0)
{
	name = "PCG";
	storage = SCSC;
}	



cg :: ~cg()
{
	if (bp      != 0 ) delete bp;
	if (perm    != 0 ) delete [] perm;
	if (invperm != 0 ) delete []invperm;
}



void cg :: initialize(const model& theModel)
{
	A.createProfile(theModel);
	b.resize(A.getNRows());
	theLinearSOE.setLinks(A,b); 
}



void cg :: info(std::ostream &of)
{
	linearSolver::info(of);
	
    of  <<  "\n Preconditioned conjugate gradient solver."
		<<  "\n Iterative solver for symmetric positive definte systems.";
	
	// ILU preconditioner missing
	if (false)
	{
		of  << "\n Incomplete LU preconditioner."
			<< "\n Drop tolerance = " << setprecision(4)  << droptol
			<< "\n Relative convergence tolerance = "     << RTOL
			<< "\n Absolute tolerance for convergence = " << ATOL
			<< "\n";
	}
	
	else
	{
		of  <<  "\n Diagonal preconditioner."
			<<  "\n Relative tolerance for convergence = " << RTOL
			<<	"\n Absolute tolerance for convergence = " << ATOL
			<<  "\n";
	}
}




bool cg :: solve(linearSOE &theSOE, double &mflops)
{
	
	// ILU preconditioner missing
	if  (false)
		return solveILU(theSOE, mflops);
	else
		return solveJacobi(theSOE, mflops);
}



/* Conjugate Gradient with Jacobi (diagonal) preconditioner.
 See Golub & Van Loan, 2nd edition, pg 534,
 See also Handbook of numerical analysis VI, Ferenz & Hughes, pg 31
 The matrix in the linearSOE is assumed to be in SCSC format, with C numbering.
 Since amux() expects a matrix in CSR format, we use atmux_() which multiplies
 the transpose of the tangent matrix (or the matrix in CSC format) with
 the vector.
*/
bool cg :: solveJacobi(linearSOE &theSOE, double &mflops)
{
    int         n;
    int         i, max_iter;
    int         dneg, dzero;
    double      dmax, dmin;
    double      res0, res, rres;
    double      alpha, beta, gamma, gammap=1.0, dtmp;
    double      rtol, atol, tiny;
    timer       t;

    extern double global_macheps;
    
    // tolerances
    rtol = RTOL;					// relative tolerance
    atol = ATOL;					// absolute tolerance
    tiny = 1e4*global_macheps;      // anything smaller, is a zero
    n    = A.getNColumns();

	if (preconditioner == "diagonal")
	{
		// extract the matrix diagonal
		if (diag.length == 0) diag.resize(n);
		A.getDiagonal(diag);
		
		// some statistics of the diagonal
		dneg = dzero = 0;
		dmin = std::numeric_limits<double>::max();
		dmax = 0.0;
		for (int k=0; k<n; k++)
		{
			dtmp = diag[k];
			
			if      (dtmp < 0.0 ) dneg++;
			else if (dtmp < tiny) dzero++;
			
			dtmp  = fabs(dtmp);
			dmax = max(dmax, dtmp);
			dmin = min(dmin, dtmp);
		}
		if (dneg >0 || dzero>0) Message("\n\t\t[ Negative diagonal terms : %3d, zero : %3d\n", dneg, dzero);
		A.setConditionNumber(dmax/dmin);
		
		// invert the diagonal: the diagonal preconditioner
		for (int k=0; k<n; k++) diag[k] = (diag[k] > 0.0) ?  1.0/diag[k] : 0.0;		
	}	
		
	
    // memory allocation
	if (p.length == 0) p.resize(n);		// the search direction
    if (r.length == 0) r.resize(n);		// the residual b-Ax
    if (v.length == 0) v.resize(n);
    if (x.length == 0) x.resize(n);
	if (z.length == 0) z.resize(n);
	    
    // initialization with x = 0
	x.setZero();
	r = b;
    res0 = r.norm();				// the initial residual norm
    res  = res0;					// absolute residual
    rres = 1.0;						// relative residual

	
    // the max number of iterations allowed. n should suffice ... in perfect arithmetic
    max_iter = min((int) (n*1.1), 1000);

    // the CG loop
    i = 0;
	t.start();
    while (i++ < max_iter && rres> rtol && res>atol)
    {
		
		/*	cout << "\n PCG iter : " << i << scientific 
			<< ", res = "  << res  << " (lim = " << atol << ")"
			<< ", rres = " << rres << " (lim = " << rtol << ")" 
			<< flush; */

		if (preconditioner == "diagonal")
			for (int k=0; k<n; k++) 
				z[k] = diag[k]*r[k];		// z <-- M^(-1) r
		else if (preconditioner == "ilu")
			continue;
			
        gamma = r.dot(z);								// gamma = r*z
        if (i == 1)
        {
            beta = 0.0;
			p = z;
        }
        else
        {
            beta = gamma/gammap;								// beta <-- gamma/oldgamma
            for (int k=0; k<n; k++) p[k] = z[k]+beta*p[k];		// p <-- z + beta*p
        }

        A.timesVector(p, v);									// v  = A p
        alpha  = gamma/p.dot(v);								// alpha = gamma/(p*v)
        for (int k=0; k<n; k++) x[k] += alpha*p[k];				// x <-- x + alpha*p
        
        if (i%50 == 0)
        {
            // every 50 iterations recompute the residual in the standard way,
            // to alleviate the accumulation of floating point errors
            A.timesVector(x, v);								// v  = A x
            for (int k=0; k<n; k++) r[k] = b[k]-v[k];			// r <-- b - A x
        }
        else
			for (int k=0; k<n; k++) r[k] -= alpha*v[k];

        res  = r.norm();
        rres = res/res0;

        gammap = gamma;									// old value of gamma
        if (i%1 == 0) 
			logger::mainlog << "\n\t[ PCG Iteration " << setw(5) << i
				<< ", residual: " << scientific << setprecision(4) << res
				<< ", time: " << setw(5) << setprecision(1) << t.seconds();
    }

    // compute the megaflops: my count:
	//   preconditioning : 1 dot prod + 1 vector op = 3 n
	//   per step        : 3 dot prod + 7 vector op + 1 Mv = 13 n + 2*nnz
    double sec      = t.seconds();
    double megaflop = 1.0e-6*static_cast<double>(3*n + i *( 13*n+2*A.getNProfileTerms() ) );

	if (sec > 0.0) 	mflops = megaflop/sec;

    if (i<max_iter-1) DebugMessage("\n\t[ PCG iteration %5d, residual: %5.4e, time: %5.1f s", i, res, sec);
    else              DebugMessage("\n\t\t[ PCG has not converged. Max. iterations (%d) reached", max_iter);

    // the solution, in x, is copied to b
    b = x;			// b <-- x
    return true;
}



// preconditioned conjugate gradient using TAUCS
bool cg :: solveILU(linearSOE &theSOE, double &mflops)
{	
	return false;
}







