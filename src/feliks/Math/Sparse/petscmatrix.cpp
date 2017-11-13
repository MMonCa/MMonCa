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
 *  petscmatrix.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 11/9/09.
 *  Copyright 2009 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "petscmatrix.h"

#include "Math/Sparse/compressed.h"
#include "Model/model.h"
#include "Model/Parts/modelpart.h"

#include "Io/logger.h"
#include <boost/foreach.hpp>
#include <string>

#ifdef WITHPETSC

petscmatrix :: petscmatrix() :
	colstart(0),
	prow(0),
	A(PETSC_NULL)
{
}


petscmatrix :: petscmatrix(const petscmatrix& pm):
	colstart(0),
	prow(0)
{
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	MatDuplicate(pm.A, MAT_COPY_VALUES, &(this->A));	
	MatSetOption(A, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE); 
	MatGetSize(A, &rows, &cols);	
}



petscmatrix :: petscmatrix(const compressedMatrix& spm)
{
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if (spm.getStorage() == SCSC)
	{
		CSCMatrix tmp(dynamic_cast<const SCSCMatrix&>(spm));
		rows = cols = tmp.getNRows();
		nnz  = tmp.getNProfileTerms();
		
		// go entry by entry filling up the Harwell-Boeing pointers
		colstart = new int[rows+1];
		prow     = new int[nnz];
		ddata    = new double[nnz]();
		
		Icopy( tmp.getColstart(), colstart, rows+1);
		Icopy( tmp.getProw()    , prow    , nnz);
		Dcopy( tmp.getData()    , ddata   , nnz);
	}
	
	if (spm.getStorage() == CSC)
	{
		rows = cols = spm.getNRows();
		nnz  = spm.getNProfileTerms();
	
		// go entry by entry filling up the Harwell-Boeing pointers
		colstart = new int[rows+1];
		prow     = new int[nnz];
		ddata    = new double[nnz];
	
		Icopy( spm.getColstart(), colstart, rows+1);
		Icopy( spm.getProw()    , prow    , nnz);
		Dcopy( spm.getData()    , ddata   , nnz);
	}
	
	else
	{
		// must take the diagonal and upper part and construct the three pointers for the full matrix
		
		
	}	

	MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, rows, cols, colstart, prow, ddata, &A);
	MatSetOption(A, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE); 
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); 
	MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);

	
	if ( spm.getStorage() == SCSC)
		MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);
}


petscmatrix :: petscmatrix(const model &m, const char* datatype) :
	colstart(0),
	prow(0)
{
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetType(A, datatype);
	createProfile(m);
	MatSetOption(A, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
}



// this step is crucial for efficiency. Petsc can allocate an empty matrix, and then
// modify the memory by adding one by one the matrix terms but this takes forever.
// a certain estimate of the number of terms per row must be given
// sparse matrices must follow CSR format, however, because matrix A is profile symmetric.
// CSR format is equivalent to CSC format, which is implemented in this function.
bool petscmatrix :: createProfile(const model &m, const char* datatype)
{
	int psSize;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &psSize);
	
	// quick return
	int ndof = m.getNDofs();
	if (ndof <= 0) return false;
	
	// gather the positions of non-zero entries on the matrix
	// a colection of vectors, one per column
		
	PetscErrorCode ierr;
	if(!strcmp(datatype, MATSEQAIJ))
    {
		MatCreate(PETSC_COMM_WORLD,&A);
		MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, ndof, ndof);
		MatSetType(A,  MATSEQAIJ);
		MatSeqAIJSetPreallocation(A, ndof, PETSC_NULL);
	}
    
	else if(!strcmp(datatype, MATAIJ))
    {				
		MatCreate(PETSC_COMM_WORLD,&A);
		
		if(psSize > 1)
        {
			MatSetType(A, MATAIJ);
			MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m.getNDofs(), m.getNDofs());
		}
		else
        {
			MatSetType(A, MATSEQAIJ);
			MatSetSizes(A,PETSC_DECIDE, PETSC_DECIDE, m.getNDofs(), m.getNDofs());
			MatSeqAIJSetPreallocation(A, m.getNDofs(), PETSC_NULL);
		}		
	}
    
	// now go element by element and fill up nonzero data
	if (rank == 0)
    {
        std::vector<int> theDofs;
        std::list<modelpart*>::const_iterator iter = m.theParts.begin();
        
		while (iter != m.theParts.end()) 
		{                        
            longvector vtmp;
            std::vector<element*>::const_iterator e = (*iter)->elements.begin();
            while (e != (*iter)->elements.end())
            {
				(*e)->getDofNumbers(theDofs);
				int   neldofs = theDofs.size();
				vtmp.resize(neldofs*neldofs);
				ierr = MatSetValues(A, neldofs, &theDofs[0], neldofs, &theDofs[0], vtmp.data, INSERT_VALUES); 
                CHKERRQ(ierr);
                ++e;
			}
			++iter;
		}		
	}
	MatGetSize(A, &rows, &cols);	
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	MatSetOption(A, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
	
    
	if (false)
	{
		PetscViewer pv;
		PetscViewerDrawOpen(PETSC_COMM_WORLD, 0, "Tangent", PETSC_DECIDE,PETSC_DECIDE,PETSC_DRAW_FULL_SIZE,PETSC_DRAW_FULL_SIZE, &pv);
		MatView(A, pv);
		getchar();
	}

	return true;
}



petscmatrix :: ~petscmatrix()
{
	MatDestroy(A);
}



void petscmatrix :: addTo(double x, int row, int col)
{
	MatSetValue(A, row, col, x, ADD_VALUES);
}



void petscmatrix ::	assemble(const matrix& block, const std::vector<int>& list)
{
	MatSetValues(A, list.size(), &list[0], list.size(), &list[0], block.start, ADD_VALUES);
}



void petscmatrix ::	finalizeAssembly()
{
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); 
	MatAssemblyEnd(  A, MAT_FINAL_ASSEMBLY);		
}



void petscmatrix :: getDiagonal(longvector& diag)
{
	Vec tmp;
	
	VecCreateSeq(PETSC_COMM_WORLD, rows, &tmp);
	MatGetDiagonal(A, tmp);
	
	double *array;
	VecGetArray(tmp, &array);
	Dcopy(array, diag.data, rows);
	VecRestoreArray(tmp, &array);
	VecDestroy(tmp);
}



double petscmatrix :: getEntry(int row, int col) const
{
	double tmp;
	MatGetValues(A, 1, &row, 1, &col, &tmp); 
	return tmp;
}



int petscmatrix :: 	getNProfileTerms() const
{
	std::cout << "getNProfileTerms not implemented" << std::endl;
	return 0;
}



int petscmatrix :: getNTermsInColumn(const int col)
{
	std::cout << "getNTermsInColumn not implemented" << std::endl;
	return 0;
}


void petscmatrix ::	info(std::ostream &of)
{
	of << "\n Petsc matrix of dimensions " << rows << " x " << cols << endl;
}



void petscmatrix :: insert(double x, int row, int col)
{
	MatSetValue(A, row, col, x, INSERT_VALUES); 
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); 
	MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);
}



void petscmatrix :: setZero()
{
	MatZeroEntries(A);
}



 void petscmatrix :: print(FILE *fp)
{
	MatView(A, PETSC_VIEWER_STDOUT_WORLD);
}



void petscmatrix :: timesVector(longvector& b, longvector& Axb)
{
	Vec vb, vprod;
	VecCreateSeq(PETSC_COMM_WORLD, rows, &vb);
	VecCreateSeq(PETSC_COMM_WORLD, rows, &vprod);
	
	double *array;
	VecGetArray(vb, &array);
	Dcopy(b.data, array, rows);
	VecRestoreArray(vb, &array);
	
	MatMult(A, vb, vprod);
	
	VecGetArray(vprod, &array);
	Dcopy(array, Axb.data, rows);
	VecRestoreArray(vprod, &array);
	
	VecDestroy(vprod);
	VecDestroy(vb);
}



void petscmatrix ::	timesDoubleArray(double b[], double Axb[])
{
	Vec vb, vprod;
	VecCreateSeq(PETSC_COMM_WORLD, rows, &vb);
	VecCreateSeq(PETSC_COMM_WORLD, rows, &vprod);
	
	double *array;
	VecGetArray(vb, &array);
	Dcopy(b, array, rows);
	VecRestoreArray(vb, &array);
	
	MatMult(A, vb, vprod);
	
	VecGetArray(vprod, &array);
	Dcopy(array, Axb, rows);
	VecRestoreArray(vprod, &array);
	
	VecDestroy(vprod);
	VecDestroy(vb);	
}



void petscmatrix :: transpose()
{
	MatTranspose(A, MAT_REUSE_MATRIX, &A);
}


#endif

