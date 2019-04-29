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
 *  direct.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on Wed May 05 2004.
 *  Copyright (c) 2004 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include "Math/LinearSolvers/direct.h"

#include "Math/Sparse/compressed.h"
#include "Math/LinearSolvers/linearsolver.h"
#include "Math/linearsoe.h"
#include "Io/logger.h"
#include "Model/model.h"
#include "Io/message.h"
#include "Math/Sparse/sparse.h"
#include "Math/vector.h"



#include <cstdio>
#include <fstream>
#include <list>
#include <string>
#include <vector>


#include <typeinfo>        

#include <superlu/supermatrix.h>
#include <superlu/slu_ddefs.h>
#include <superlu/slu_util.h>

using namespace std;

 
direct :: direct() :
	L(0),
	U(0),
	colPerm(0),
	rowPerm(0),
	elTree(0),
	isPermuted(false)
{
	name = "unsymmetric direct linear equation solver";
	storage = CSC;	
}


direct :: direct(const commandLine &cl) :
	linearSolver(cl),
	L(0),
	U(0),
	colPerm(0),
	rowPerm(0),
	elTree(0),
	isPermuted(false)
{
	name = "unsymmetric direct linear equation solver";
	storage = CSC;
}



direct :: ~direct()
{
	if (L != 0)  
	{
		Destroy_SuperNode_Matrix(L);
		L = 0;
	}
	if (U != 0) 
	{
		Destroy_CompCol_Matrix(U);
		U = 0;
	}
	
	if (colPerm != NULL) delete [] colPerm;
	if (rowPerm != NULL) delete [] rowPerm;
	if (elTree  != NULL) delete [] elTree;
}



void direct :: initialize(const model& theModel)
{
	A.createProfile(theModel);
	b.resize(A.getNColumns());
	theLinearSOE.setLinks(A,b); 
}




/* Print a bried description of the solver */
void direct :: info(ostream &of)
{
	linearSolver::info(of);
	of  << "\n Direct solver for general unsymmetric systems of equations."
		<< "\n SuperLU 3.0 version. See http://crd.lbl.gov/~xiaoye/SuperLU\n";
}




bool direct ::  solve(linearSOE &theSOE, double &mflops)
{
	bool r;
	
	sparseMatrix &m = theSOE.getA();
	if (typeid(m) != typeid(CSCMatrix))
	{
		Message("\n Error, matrix must be csc to solve with general, unsymmetric, direct solver");
		r = false;
	}
	
	//else if (m.getDatatype() == DATA_DCOMPLEX)	
	//	r = solveComplex(theSOE, mflops);
	else	
		r = solveReal(theSOE, mflops);
	
	return r;
}




/* This routine is hardwired with the CSC storage functions. This should
be separated in the future.
Return 1: if ok, -1 if error
*/
bool direct :: solveReal(linearSOE &theSOE, double &mflops)
{
    SuperMatrix  A, B, X;
    char         equed[1];
    int          nrhs;
    int          info=0;
    int         ndata, ncol,lwork(0);
    void        *work=NULL;
    double      *Bdata=NULL;
    double       ferr[1], berr[1], rpg, rcond;
    bool         ret = false;

    // recover matrix and rhs vector
	// we know that the matrix is of type CSCMatrix
	CSCMatrix  &m = dynamic_cast<CSCMatrix&> (theSOE.getA());
    longvector &v = theSOE.getB();
	
    // recover sparseMatrix structure for CSC storage
    ncol      = m.getNColumns();
    ndata     = m.getNProfileTerms();
    nrhs      = 1;

    // default options, common for all tasks
    superlu_options_t options;
    set_default_options(&options);
	
    options.Trans           = NOTRANS;
    options.DiagPivotThresh = 0.0;      //0.0, no pivoting. Use for SPD matrices
    options.IterRefine      = NOREFINE;
    options.PrintStat       = NO;
    options.SymmetricMode   = YES;
	options.ConditionNumber = YES;
	options.PivotGrowth     = NO;
	options.Equil           = NO; equed[0] = 'N';
    
    // if the LU factors have to be refactorized
    if ( ! m.isFactorized() )
    {
        // if the LU factors were once computed but they are obsolete, free their content, but not the structure
        if (L != 0)
        {
            Destroy_SuperNode_Matrix(L);
            Destroy_CompCol_Matrix(U);
            options.Fact = SamePattern;
        }

        // if the LU factors do not exist, allocate the structure and indicate superlu to do everything
        else
        {
            L = (SuperMatrix *) malloc(sizeof(SuperMatrix));
            U = (SuperMatrix *) malloc(sizeof(SuperMatrix));
            options.Fact = DOFACT;
        }
    }

    // the matrix factorization is correct. Reuse it.
    else
        options.Fact = FACTORED;


    // If the matrix has never been preorderer allocate space & set options
    if ( !isPermuted )
    {
        // allocate space the permutation matrices & elimination tree
        colPerm = new int[ncol];
        rowPerm = new int[ncol];
        elTree  = new int[ncol];

        /* tell what ordering to use
            NATURAL=no permutation,
            MMD_ATA=mininum degree on A'*A,
            MMD_AT_PLUS_A= min. degree on A+A',
            COLAMD=column approx. min. degree */
        if (ncol < 50)
            options.ColPerm = NATURAL;
        else
            options.ColPerm = MMD_AT_PLUS_A;

		isPermuted = true;
    }
    else
    {
        // reuse permutation matrices already computed

        // tell the solver not to compute the permutation matrices, but to reuse
        options.ColPerm = MY_PERMC;
    }


    // Copy of the rhs of the problem 
    Bdata = (double *) malloc(ncol* sizeof(double));
    Dcopy(v.data, Bdata, ncol);
	
    dCreate_CompCol_Matrix(&A, ncol, ncol, ndata, m.ddata, m.prow, m.colstart, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(  &B, ncol,       nrhs,  Bdata,    ncol,              SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(  &X, ncol,       nrhs,  v.data,   ncol,              SLU_DN, SLU_D, SLU_GE);

    // initialize statistics
	SuperLUStat_t   stat;
    StatInit(&stat);

   // TODO: New argument introduced in SUPERLU 5 GlobalLU_t, I added for compiling reasons
   // but it will not be used. This should be reviewed for guarantiee no impact is dervied
   // from this action
   GlobalLU_t to_throw;


    // the main solution function
	mem_usage_t  mem_usage;
    dgssvx(&options, &A, colPerm, rowPerm, elTree, equed, NULL, NULL,
           L, U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &to_throw, &mem_usage, &stat, &info);
	
	m.setFactorized();


	// errors
	if (info != 0)    WarningMessage("Error in SuperLU. Code = %d", info);
	if (info == ncol+1) 
		WarningMessage("Double superLU found matrix to be singular upto machine precision");
	
	
	
    /* Total megaflops, as the average of the factor and solution operations */
    if (stat.utime[FACT]+stat.utime[SOLVE] > 0.0)
        mflops = (stat.ops[FACT]+stat.ops[SOLVE])*1e-6/(stat.utime[FACT]+stat.utime[SOLVE]);
    else
        mflops = -1.0;
    if (options.PrintStat == YES) StatPrint(&stat);
    StatFree(&stat);


    // store matrices for future usage
    // store the new matrices in the sparse matrix data
    m.factorized  = true;
    m.condition   = 1.0/rcond;

    //this only destroys the superLU structure, leaving the data
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&X);
	
	//this deallocates everything
    Destroy_Dense_Matrix(&B);

    return ret;
}





/* This routine is hardwired with the CSC storage functions. This should
be separated in the future.
Return 1: if ok, -1 if error

static int  SolveComplex(linearSOE theSOE, double *mflops)
{
    SuperMatrix  A, B, *L=NULL, *U=NULL, X;
    sparseMatrix m;
    vector      v;
    char         equed[1];
    int         *perm_r, *perm_c, nrhs;
    int          info=0, *etree;
    int         *prow, *colstart, ndata, ncol;
    void        *work=NULL;
    complex double  *thedata=NULL, *Bdata=NULL;
    double       ferr[1], berr[1], rpg, rcond;
    int          ret = -1;
	
    mem_usage_t       mem_usage;
    SuperLUStat_t     stat;
    superlu_options_t options;
	
	// functions without prototype
	void zCreate_CompCol_Matrix();
	void zCreate_Dense_Matrix();
	void zgssvx();
	
	
    // recover matrix and rhs vector
    m = GetAInLinearSOE(theSOE);
    v = GetBInLinearSOE(theSOE);
	
    // recover sparseMatrix structure for CSC storage
    ncol      = m->cols;
    prow      = m->ipointer[0];
    colstart  = m->ipointer[1];
    thedata   = (complex double *)m->ddata;
    ndata     = colstart[ncol];
    nrhs      = 1;
    L         = m->L;
    U         = m->U;
	
    // default options, common for all tasks
    set_default_options(&options);
    options.Trans           = NOTRANS;
    options.DiagPivotThresh = 0.0;      //0.0, no pivoting. Use for SPD matrices
    options.IterRefine      = NOREFINE;
    options.PrintStat       = NO;
    options.SymmetricMode   = YES;
    options.Equil           = NO; equed[0] = 'N';
    options.ConditionNumber = YES;
    options.PivotGrowth     = NO;
	
    // if the LU factors have to be refactorized
    if (m->factorized != 1)
    {
        // if they have been previously computed, free their content, but not the structure
        if (L != NULL)
        {
            Destroy_SuperNode_Matrix(L);
            Destroy_CompCol_Matrix(U);
            options.Fact = SamePattern;
        }
		
        // if they are completely fresh, just allocate the structure
        else
        {
            L = (SuperMatrix *) Allocate(sizeof(SuperMatrix));
            U = (SuperMatrix *) Allocate(sizeof(SuperMatrix));
            options.Fact = DOFACT;
        }
    }
	
    // the matrix factorization is correct. Reuse it.
    else
        options.Fact = FACTORED;
	
	
    // If the matrix has never been preorderer allocate space & set options
    if (m->permuted != 1)
    {
        // allocate space the permutation matrices & elimination tree
        perm_c = intAllocate(ncol);
        perm_r = intAllocate(ncol);
        etree  = intAllocate(ncol);
		
        // tell what ordering to use
        //    NATURAL=no permutation,
        //    MMD_ATA=mininum degree on A'*A,
        //    MMD_AT_PLUS_A= min. degree on A+A',
        //    COLAMD=column approx. min. degree 
        if (ncol < 50)
            options.ColPerm = NATURAL;
        else
            options.ColPerm = MMD_AT_PLUS_A;
    }
    else
    {
        // reuse permutation matrices already computed
        perm_c = m->ipointer[2];
        perm_r = m->ipointer[3];
        etree  = m->ipointer[4];
		
        // tell the solver not to compute the permutation matrices, but to reuse
        options.ColPerm = MY_PERMC;
    }
	
	
    Copy of the rhs of the problem 
    Bdata = dcomplexAllocate(2*ncol);
    Zcopy((double complex *)v->data, Bdata, ncol);
	
    zCreate_CompCol_Matrix(&A, ncol, ncol, ndata, m->ddata, prow, colstart, SLU_NC, SLU_Z, SLU_GE);
    zCreate_Dense_Matrix(  &B, ncol,       nrhs,  Bdata,    ncol,           SLU_DN, SLU_Z, SLU_GE);
    zCreate_Dense_Matrix(  &X, ncol,       nrhs,  v->data,  ncol,           SLU_DN, SLU_Z, SLU_GE);
	
		
    // initialize statistics
    StatInit(&stat);
	
    // the main solution function
    zgssvx(&options, &A, perm_c, perm_r, etree, equed, NULL, NULL,
           L, U, work, 0, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);
		
	
	// what to do with errors
	if (info == ncol+1) 
		WarningMessage("Complex superLU found matrix to be singular upto machine precision");
	
	//zPrint_CompCol_Matrix("tangent", &A);
	// zPrint_Dense_Matrix("rhs", &B);
	// zPrint_Dense_Matrix("sol", &X);
	
	
    // Total megaflops, as the average of the factor and solution operations 
    if (stat.utime[FACT]+stat.utime[SOLVE] > 0.0)
        *mflops = (stat.ops[FACT]+stat.ops[SOLVE])*1e-6/(stat.utime[FACT]+stat.utime[SOLVE]);
    else
        *mflops = -1.0;
    if (options.PrintStat == YES) StatPrint(&stat);
    StatFree(&stat);
	
	
    // store matrices for future usage
    // store the new matrices in the sparse matrix data
    m->permuted    = 1;
    m->factorized  = 1;
    m->L           = (void *) L;
    m->U           = (void *) U;
    m->condition   = 1.0/rcond;
    m->ipointer[2] = perm_c;
    m->ipointer[3] = perm_r;
    m->ipointer[4] = etree;
	
    //this only destroys the superLU structure & leaves the data
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&X);
    Destroy_Dense_Matrix(&B);
	
    return ret;
}

*/
			
		
cholesky :: cholesky(const commandLine &cl) :
	linearSolver(cl),
	L(0),
	perm(0),
	invperm(0),
	isPermuted(false),
	xp(0),
	bp(0)
{
	name = "Cholesky solver for linear systems of equations with symmetric, positive definite matrices.";
	storage = SCSC;

	
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if       ( uc.keyword() == "nprocessors" ) nprocessors = (int) uc.value();
	}
}
		
		
cholesky :: cholesky() :
L(0),
perm(0),
invperm(0),
isPermuted(false),
xp(0),
bp(0)
{
	name        = "Cholesky solver for linear systems of equations with symmetric, positive definite matrices.";
	storage     = SCSC;
	nprocessors = 1;
}

		
cholesky :: ~cholesky()
{
}




/* Print a bried description of the solver */
void cholesky :: info(ostream &of)
{
	linearSolver::info(of);

	of  << "\n Direct Cholesky solver not implemented\n";
}



void cholesky :: initialize(const model& theModel)
{
	A.createProfile(theModel);
	b.resize(A.getNProfileTerms());
	theLinearSOE.setLinks(A,b); 
}


bool cholesky ::  solve(linearSOE &theSOE, double &mflops)
{
	logger::mainlog  << "\n Error, Cholesky solver not implemented";
	logger::warnings << "\n Error, Cholesky solver not implemented";
	
	return false;
}





bool cholesky :: solveReal(linearSOE &theSOE, double &mflops)
{
	return false;
}


