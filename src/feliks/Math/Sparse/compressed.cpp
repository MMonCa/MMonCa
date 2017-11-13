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
/* functions for CSC sparse matrices, symmetric and unsymmetric
 *
 * ignacio romero, december 2002
 *
 * Description of Harwell-Boeing / Compressed Column Storage:
 * 
 * For a NxN sparse matrix we define:
 *
 *   nnz
 *   ddata       = data     : vector of nnz doubles, the nonzero entries in the matrix.
 *   ipointer[0] = prow     : vector of nnz ints, the row index of the elements Aij in ddata.
 *   ipointer[1] = colstart : vector of N+1 ints, the start of each column.
 *   ipointer[2] = an array of ndof ints, the column permutations for SuperLU.
 *   ipointer[3] = an array of ndof ints, the row permutations for SuperLU.
 *   ipointer[4] = a int* of ndof, the elimination tree etree for SuperLU.
 *   ipointer[5] = an array of ndof ints, the positions of the diagonal terms in ddata 
 *
 *
 *   Note: even this is not needed in general, the routines in this file assume that
 *         the nonozero elements in a column are stored with increasing row numbers.
 * 
 *
 *   Example: (taken from superlu manual)
 *                          
 *                        [19  0  21  21  0]
 *                        [12 21   0   0  0]
 *                   A  = [ 0 12  16   0  0]
 *                        [ 0  0   0   5 21]
 *                        [12 12   0   0 18] 
 *	
 *                   n     = 5
 *                   ddata = [19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18]
 *                   nnz   = 12
 *                   prow  = [ 0,  1,  4,  1,  2,  4,  0,  2,  0, 3,  3,  4]
 *                   colstart = [0, 3, 6, 8, 10, 12]
 *
 */

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <iomanip>
#include <cstddef>
#include <fstream>
#include <math.h>
#include <ostream>
#include <stdexcept>
#include <typeinfo>
#include "boost/foreach.hpp"

#include "Elements/evalspot.h"
#include "Elements/element.h"

#include "Math/Sparse/sparse.h"
#include "Math/Sparse/compressed.h"
#include "Math/linearsoe.h"
#include "Math/LinearSolvers/linearsolver.h"
#include "Math/LinearSolvers/direct.h"
#include "Math/LinearSolvers/cg.h"

#include "Math/vector.h"
#include "Math/matrix.h"

#include "iohb.h"

#include "Model/model.h"
#include "Model/Parts/modelpart.h"

#include "Io/message.h"

using namespace std;



struct index_real_value_pair_t
{
	int index;
	double value;
};



/**
 * Compares two (index,value) pairs by their indices.
 */
static int
compare_index_real_value_pairs (const void *pa, const void *pb)
{
	const struct index_real_value_pair_t* a = (const struct index_real_value_pair_t*) pa;
	const struct index_real_value_pair_t* b = (const struct index_real_value_pair_t*) pb;
	
	if (a->index > b->index)
		return 1;
	else if (a->index < b->index)
		return -1;
	else
		return 0;
}





/** 
 * Takes column number j of the given sparse matrix A, and copies all the 
 * indices (and values, if there are any values) of that column into "col".
 *
 * @param A [IN]    Sparse matrix in CSC format
 * @param j [IN]    Index of the column to copy (zero-based indices)
 * @param col [OUT] An array of either index_real_value_pair_t, 
 *                  index_complex_value_pair_t or int, depending on whether 
 *                  the value type of the matrix is REAL, COMPLEX, or PATTERN,
 *                  respectively.
 * @param max_nnz [IN]  The max number of nonzeros in any column of A
 */
void compressedMatrix :: copy_col2pairs (int j, void* col, int max_nnz)
{
	int a = colstart[j];
	int b = colstart[j+1];
	int nnz = b - a;
	int k;
	
	assert (nnz <= max_nnz);
	
	struct index_real_value_pair_t* _col = (struct index_real_value_pair_t*) col;
	const double* const values = (const double* const) (ddata);
		
	for (k = 0; k < nnz; k++)
	{
		assert ((a+k) < b);
		_col[k].index = prow[a+k];
		_col[k].value = values[a+k];
    }
}


/** 
 * Given a sparse matrix A in CSC format, a column number j, and a list of
 * (index,value) pairs, copies the list of (index,value) pairs back into that
 * column of the matrix.
 */
void compressedMatrix :: copy_pairs2col (const void* col, int max_nnz, int j)
{
	int a = colstart[j];
	int b = colstart[j+1];
	int nnz = b - a;
	
	int k;
	
	assert (nnz <= max_nnz);
	
	struct index_real_value_pair_t* _col = (struct index_real_value_pair_t*) col;
	double* values = (double*) (ddata);
		
	for (k = 0; k < nnz; k++)
	{
		assert ((a+k) < b);
		prow[a+k] = _col[k].index;
		values[a+k] = _col[k].value;
	}
}





compressedMatrix :: compressedMatrix(model &m, char datatype) :
	prow(0),
	colstart(0),
	pdiag(0),
	nnz(0)
{
}


compressedMatrix :: compressedMatrix():
	prow(0),
	pdiag(0),
	colstart(0),
	nnz(0)
{
}


// constructor that reads data from file
compressedMatrix :: compressedMatrix (const char *filename, bool symmetric) :
	prow(0),
	colstart(0),
	pdiag(0)
{
	char    title[73];       // Title.
	char    name[9];         // Name.
	
	// Declaring incremental variables for loops.
	int    i, j;
	
	// Variables for the second line
	int    lintot; // total number of lines in the file 
	int    linptr; // number of lines for pointers
	int    linind; // number of lines for row indices
	int    linval; // number of lines for numerical values
	int    linrhs; // number of lines for right-hand sides
	
	// Variables for the third line
	int    npcol;
	int    fpcol;
	int    nirow;
	int    firow;
	int    nval;
	int    fval;
	char   c;
	char   num[81];
	double value;
	
	// Opening file.
	ifstream file(filename);
	if (!file) throw runtime_error("Error opening file");
	
	// Reading the first line.
	file.get((char*)title, 73, '\n');
	file.get((char*)name, 9, '\n');
	do file.get(c); while (c!='\n'); 
	
	// Reading the second line.
	file >> lintot >> linptr >> linind >> linval >> linrhs;
	do file.get(c); while (c!='\n'); 
	if ((linptr < 1) || (linind < 1)) throw runtime_error("Parameter error");
	
	
	// Reading the third line.
	// Consists of 3 chars [R|C|P][S|U|H|Z|R][A|E]
	// real|complex|just pattern, no data
	// symmetric|unsymmetric|hermitian|skew|rectangular
	// assembled|element blocks
	char    type[4];         // Matrix type.
	file.get((char*)type,4,'\n');
	file >> rows >> cols >> nnz;
	do file.get(c); while (c!='\n'); 
	cout << "input format = " << type << "\n";
	if (((type[0] != 'R') && (type[0] != 'C')) || 
		(symmetric && (type[1] != 'S')) || 
		(type[2] != 'A')) 
	{
		cout << "Matrix type in file does not match data structure" << "\n";
		throw runtime_error("WRONG MATRIX TYPE");
	}
	else if ((rows < 1) || (cols < 1) || (nnz < 1)) {
		throw runtime_error("PARAMETER ERROR");
	}
	
	// Reading the fourth line.
	ReadFormat(file, npcol, fpcol);
	ReadFormat(file, nirow, firow);
	ReadFormat(file, nval, fval);
	do file.get(c); while (c!='\n'); 
	if ((fpcol<1) || (firow<1) || (fval<1)) throw runtime_error("WRONG DATA TYPE");
	
	
	// Skipping the fifth line.
	if (linrhs) {
		do file.get(c); while (c!='\n'); 
		runtime_error("RHS IGNORED");
	}
	
	
	// Reading column pointers.
	pdiag    = new int[cols];
	colstart = new int[cols+1];
	fpcol++;
	i = 0;
	while ((i <= cols) && (file.get((char*)num,fpcol,'\n'))) {
		colstart[i++] = atoi((char*)num)-1;
		if (!(i%npcol)) do file.get(c); while (c!='\n'); 
	}
	if (i%npcol) do file.get(c); while (c!='\n'); 
	
	if (i <=cols) {
		throw runtime_error("UNEXPECTED EOF");
	}
	
	
	// Reading row indices.
	prow = new int[nnz];
	firow++;
	i = 0;
	while ((i < nnz) && (file.get((char*)num,firow,'\n'))) {
		prow[i++] = atoi((char*)num)-1;
		if (!(i%nirow)) do file.get(c); while (c!='\n'); 
	}
	if (i%nirow) do file.get(c); while (c!='\n');   
	if (i < nnz) throw runtime_error("UNEXPECTED EOF");
	
	
	
	// Reading matrix elements.
	fval++;
	ddata = new double[nnz];
	i = 0;
	j = 0;
	while ((i < nnz) && (ReadEntry(file, nval, fval, j, value))) {
		ddata[i++] = value;
	}  
	if (j%nval) do file.get(c); while (c!='\n'); 
	
	if (i < nnz) {
		throw runtime_error("UNEXPECTED EOF");
	}
	
	// fill diagonal counter
	for (int i=0; i<cols; i++)
	{
		int j;
		for (j = colstart[i]; j < colstart[i+1]; j++)
			if ( prow[j] == i)
			{
				pdiag[i] = j;
				break;
			}
		
		if (j == colstart[i+1])
			printf("\n Error in computing diagonal. Entry is empty");
	}
	
	
	// Closing file and reporting success.
	file.close();
	
	// expand the matrix if it must be non-symmetric
	if (!symmetric && type[1] == 'S')
		expand_symmetric_storage();
	
}







compressedMatrix :: ~compressedMatrix()
{
	if (prow != 0)		delete [] prow;
	if (colstart != 0)	delete [] colstart;
	if (pdiag != 0)		delete [] pdiag;	
}



void compressedMatrix :: convertDouble(char* num)
{
	char* pd;
	
	pd = strchr(num,'D');
	if (pd) *pd = 'E';
	pd = strchr(num,'d');
	if (pd) *pd = 'E';
};



// change the numbering of integer pointer arrays to C type 
void compressedMatrix :: convertToC()
{
	if (numbering == NUM_C) return;
	
	int k, n = cols;
	
	for (k=0; k<nnz; k++) prow[k]--;
	for (k=0; k<n+1; k++) colstart[k]--;    
	
	numbering = NUM_C;
}




// change the numbering of integer pointer arrays to fortran type 
void  compressedMatrix :: convertToFortran()
{
	if (numbering == NUM_FORTRAN) return;
	
	int k, n = cols;
	
	for (k=0; k<nnz; k++) prow[k]++;
	for (k=0; k<n+1; k++) colstart[k]++;
	
	numbering = NUM_FORTRAN;
}



void compressedMatrix :: dump(string format)
{    
	if (format == "matrixmarket" || format == "mm")
	{
		ofstream os;
		os.open("feliks.tangent.mm");
		
		os << "%%MatrixMarket matrix coordinate real general" << "\n";
		os << "%% standard establishes 1 based " << "\n";
		os << "%% Nrows = " << rows << ", Ncols = " << cols << ", nnz = " << nnz << "\n";
		os << rows << " " << cols << " " << nnz << "\n";
		
		int cc = 0;
		for (int i=0; i<nnz; i++)
		{
			if (i == colstart[cc+1]) cc++;
			os  << setw(5)  << prow[i]+1 << " " 
			<< setw(5)  << cc+1 << " " 
			<< setw(18) << setprecision(14) << ddata[i] << "\n";
		}
		
		os.close();
	}
	
	else if (format == "harwellboeing" || format == "harwell-boeing" || format == "hb")
	{
		/* This uses the Matrix Market subroutines to write a hwb sparse matrix
		 to a file. It produces a "nicer" output that the one above */
		
		
		int  space;
		char filename[]="feliks.tangent.hwb";
		char title[]="feliks tangent";
		char key[]="feliks key";
		char type[]="RUA";
		char Ptrfmt[]="(13I7)";
		char Indfmt[]="(16I7)";
		char Valfmt[]="(3E26.18)";
		char Rhsfmt[]="(3E26.18)";
		
		// get the size for the row and column pointers and leave one space for readability
		space = (int) log10((double)nnz) + 2;
		Ptrfmt[4] = space + toascii('0');
		Indfmt[4] = space + toascii('0');
		
		writeHB_mat_double(filename, rows, cols, nnz, colstart, prow,
						   ddata, 0, NULL, NULL, NULL,
						   title, key, type, 
						   Ptrfmt, Indfmt, Valfmt, Rhsfmt,
						   NULL);
	}
}




void compressedMatrix :: getDiagonal(longvector& diag)
{
	double *D  = diag.data;
	for (int k=0; k<cols; k++) D[k] = ddata[pdiag[k]];
}



int compressedMatrix :: getNTermsInColumn(const int col)
{
	return colstart[col+1] - colstart[col];
}



/*Function to read the format of the variables in line 4*/
void compressedMatrix :: ReadFormat( ifstream& file, int& n, int& fmt)
{
	char c;
	
	do file.get(c); while ((c != '(') && (c!='\n'));
	file >> n;
	file.get(c);
	while ((c!='I') && (c!='i') && (c!='E') && (c!='e') &&
		   (c!='D') && (c!='d') && (c!='\n')) {
		do file.get(c); while ((c != ',') && (c!='\n'));  
		file >> n;
		file.get(c);
	}
	if ((c==')')||(c=='\n')) { // Reading error!
		fmt = 0;
	}
	else {
		file >> fmt;
	}	
} 




bool compressedMatrix :: ReadEntry( ifstream& file, int nval, int fval, int& j, double& val)
{
	
	char num[81];
	char c;
	
	if (file.get((char*)num,fval,'\n')) 
	{
		convertDouble((char*)num);
		val = atof((char*)num);
		if (!((++j)%nval)) do file.get(c); while (c!='\n'); 
		return true;
	}
	else {
		return false;
	}
	
} 




void compressedMatrix :: timesVector(longvector& b, longvector& Axb)
{	
	// some checks
	if ( cols != b.length)
		printf("matrix and vector dimensions do not match in compressedMatrix::times vector");
	
	else if ( b.length != Axb.length)
		printf("input and output vector have different dimensions");
	
	else if (datatype == DATA_DOUBLE)
		timesVectorDouble(b, Axb);
	
	
	else printf("unknown type in compressedMatrix::times vector");
}



/**
 * Sorts the row indices within each column of the sparse matrix A.
 */
void compressedMatrix :: csc_matrix_sort_rowidx()
{
	int j;
	int max_nnz = 0;  /* Will be the max # of nonzeros in each column */
	
	
	if (cols <= 0) return;
	
	/* Find the max # of nonzeros in each column.  This lets us allocate 
	 * workspace large enough to accomodate any column (otherwise we would
	 * have to allocate it separately for _each_ column, which would waste
	 * a lot of malloc calls). */
	max_nnz = colstart[1] - colstart[0];
	assert (colstart[1] - colstart[0] >= 0);
	for (j = 1; j < cols; j++)
    {
		int nnz = colstart[j+1] - colstart[j];
		if (nnz < 0)
			
		{
			printf("*** csc_matrix_sort_rowidx: at column %d, colstart%d]=" 
					   "%d < colstart[%d]=%d ***\n", 
					   j, j+1, colstart[j+1], j, colstart[j]);
		}
		max_nnz = std::max<double>(nnz, max_nnz);
    }
	

	/* Sort the row indices in each column, reordering the corresponding values accordingly. */
		/* col: Workspace for sorting (index,value) pairs */
		struct index_real_value_pair_t* col = (index_real_value_pair_t*) malloc (max_nnz * sizeof (struct index_real_value_pair_t));
		
		
		for (j = 0; j < cols; j++)
		{
			int nnz = colstart[j+1] - colstart[j];
			
			if (nnz < 0)
			{
				printf("*** csc_matrix_sort_rowidx: At column %d, nnz = %d < 0 ***\n", j, nnz);
				exit(1);
			}
			
			/* Copy the (index,value) pairs into the temp location "col". */
			copy_col2pairs (j, col, max_nnz);
			
			/* Do the sorting in "col". */
			qsort (col, nnz, sizeof (struct index_real_value_pair_t), compare_index_real_value_pairs);
			
			/* Copy the (index,value) pairs back from "col" into the matrix. */
			copy_pairs2col (col, max_nnz, j);
		}
		
		/* Free up the temp storage. */
		free(col);
}



//"borrowed" from BeBOP library
void compressedMatrix :: expand_symmetric_storage () 
{
	/* Code borrowed from Rich Vuduc's Sparsity-0.1 suite (spmv_util.c,
	 * spmv_expand_symmcsr) and adapted for CSC-format matrices (the original was
	 * for CSR-format matrices). */
	
	int *cur_col_nnz = NULL; /* # of nonzeros currently in each col of the matrix A */
	int *new_col_nnz = NULL; /* # of nonzeros in each col of the expanded matrix */
	int  new_nnz;            /* Total number of nonzeros in the expanded matrix */	
	int *colptr = NULL;      /* Will be the new colptr array of A */
	int *rowidx = NULL;      /* Will be the new rowidx array of A */
	

	if (rows != cols)
    {
		cout << "Error, matrix is not square" << endl;
		return;
    }
	
	else if (storage == CSC)
    {
		// Not an error -- it just means that we don't have any work to do
		return ;
    }
		
	cur_col_nnz = new int[cols];
	new_col_nnz = new int[cols];
	
	/* 
	 * Scan A and count how many new nonzeros we will need to create. 
	 * When we are done scanning:
	 *
	 * cur_col_nnz[j] == # of nonzeros in column j of original matrix
	 * new_col_nnz[j] == # of nonzeros to be stored in col j of matrix
	 * new_nnz        == Total # of nonzeros to be stored in matrix
	 */
	new_nnz = 0;
	
	// First, count the number of nonzeros currently in each column.
	//bebop_log (2, "Counting # nnz currently in each column...");
	int j;                   /* Current column index */
	for (j=0; j<cols; j++)
    {
		cur_col_nnz[j] = colstart[j+1] - colstart[j];
		new_col_nnz[j] = cur_col_nnz[j];
		new_nnz += new_col_nnz[j];
    }
	
	
	/* Now figure out how many elements we need to reflect across the diagonal. */
	for (j=0; j<cols; j++)
    {
		int k;
		for (k = colstart[j]; k < colstart[j+1]; k++)
		{
			int ii = prow[k];
			
			/* Reflect off-diagonal elements across the diagonal; don't count
			 * diagonal elements twice.  Element (ii,j) goes into column ii (it
			 * is reflected into element (j,ii)). */
			if (ii != j)
			{
				new_col_nnz[ii]++;
				new_nnz++;
			}
		}
    }
	
	/* Initialize new arrays.  These will hold the expanded version of the matrix. */
	colptr = new int [cols + 1];
	rowidx = new int [new_nnz];
	
	/* Copy old colptr into new colptr */
	memcpy (colptr, colstart, (cols + 1) * sizeof (int));
	
	/* 
	 * Initialize column pointers in A.  After we are done:
	 *
	 * colptr initialized to the correct, final values.
	 * new_col_nnz[j] reset to be equal to cur_col_nnz[j].
	 */
	for (j = 1; j <= cols; j++)
    {
		colptr[j] = colptr[j-1] + new_col_nnz[j-1];
		new_col_nnz[j-1] = cur_col_nnz[j-1];
    }
	colptr[cols] = new_nnz;
	
		
	/* 
	 * Complete expansion of the matrix to full storage.  After we are done:
	 *
	 * (colptr, rowidx, values) is the full-storage equivalent of A.
	 * new_col_nnz[j] == # of nonzeros in col j of the (expanded) matrix.
	 */
		double* values = new double[new_nnz];
		for (j = 0; j < cols; j++)
		{
			int cur_nnz = cur_col_nnz[j]; /* number of nonzeros in current row of old matrix */
			int k_cur = colstart[j];    /* current position in old matrix */
			int k_new = colptr[j];       /* current position in expanded matrix */
			
		
			
			/* Copy current nonzeros from old matrix to new matrix */
			memcpy (&rowidx[k_new], &(prow[k_cur]), cur_nnz * sizeof (int));
			memcpy (values + k_new, ddata + k_cur, cur_nnz * sizeof (double));
			
			/* Fill in the symmetric "missing" values */
			while (k_cur < colstart[j+1])
			{
				/* Nonzero of A */
				int ii = prow[k_cur];
				
				if (ii != j) /* If not a diagonal element */
				{
					/* Get the current element from the old matrix */
					double a = ddata[k_cur];
					
					/* Position of this transposed element in A */
					k_new = colptr[ii] + new_col_nnz[ii];
					
					/* Store the nonzero */
					rowidx[k_new] = j;
					values[k_new] = a;
					
					/* Update so that the next element stored at column ii will
					 * appear at the right place */
					new_col_nnz[ii]++;
				}
				
				k_cur++;
			}
		}
		
		/* Cleanup */
		delete [] cur_col_nnz;
		delete [] new_col_nnz;
		
		/* Free the old arrays and assign the new ones */
		delete [] colstart;
		delete [] prow;
		delete [] ddata;
		
		colstart = colptr;
		prow     = rowidx;
		ddata    = values;
		nnz      = new_nnz;

		
	/* Sort the row indices in the matrix */
	csc_matrix_sort_rowidx();
	
	/* Now A uses unsymmetric storage */
	storage = CSC;
}




/* this function will dump the tangent to a postscript file where the sparsity pattern could 
 be inspected. Use Sparskit function pspltm()
 

extern "C" { void pspltm_(int*, int*, int*, int*, int*, int*, float*, char*, int*, int*, int*); }
void compressedMatrix :: plot()
{
	
	int one=1, zero=0, unit;
	//char title[]="FELIKS tangent";
	char cm[2];
	float size=20.0; 
	
	unit = 20;
	cm[0] = 'c'; cm[1] = 'm';
	
	pspltm_(&rows, &cols, &one, prow, colstart, &zero, &size, cm, &zero, NULL, &unit);
}

*/



void compressedMatrix :: setZero()
{
	if      (datatype == DATA_DOUBLE)    Dzero(ddata, nnz);
	//else if (datatype == DATA_DCOMPLEX)  Zzero((complex double *) ddata, nnz);
}






/* functions for SCSC sparse matrices
 *
 * ignacio romero, december 2002
 *
 * Description of Harwell-Boeing / Symmetric Compressed Column Storage:
 * We store the lower half of the matrix, for faster matrix x vector operations
 * 
 * For a NxN sparse matrix we define:
 *
 *   nnz
 *   ddata       = data     : vector of nnz doubles, the nonzero entries in the matrix.
 *   prow     : vector of nnz ints, the row index of the elements Aij in ddata.
 *   colstart : vector of N+1 ints, the start of each column.
 *
 *   Note: even though it is not needed in general, the routines in this file assume that
 *         the nonozero elements in a column are stored with increasing row numbers.
 * 
 *
 *   Example:
 *                          
 *                        [19  0  21  11  0]
 *                        [ 0 22  12   0  0]
 *                   A  = [21 12  16   7  0]
 *                        [11  0   7   5 21]
 *                        [ 0  0   0  21 18] 
 *
 *                   n     = 5
 *                   ddata = [19, 22, 21, 12, 16, 11, 7, 5, 21, 18]
 *                   nnz   = 10
 *                   prow  = [ 0,  1,  0,  1,  2,  0, 2, 3,  3 , 4]
 *                   colstart = [0, 1, 2, 5, 8, 10]
 *
 */



void compressedMatrix :: transpose()
{
	cout << "funciton transpose is no longer defined, to avoid using sparskit" << endl;
	
	/*
	int* iwk = new int[ getNProfileTerms() ];
	int  ierr(0);
	int  nrows( getNRows() );
	int  ncols( getNColumns() );
	
	if ( numbering == NUM_C) 
	{
		this->convertToFortran();
		transp_( &nrows, &ncols, getData(), getProw(), getColstart(), iwk, &ierr);
		this->convertToC();
	}
	else
		transp_( &nrows, &ncols, getData(), getProw(), getColstart(), iwk, &ierr);
	delete [] iwk;
	 */
}



SCSCMatrix :: SCSCMatrix() :
compressedMatrix()
{
}



/*	function to allocate the memory space for a matrix in the Harwel-Boeing (aka compress sparse column) format.
 *	We might consider matrices for models with 0 dofs. These are not common but sometimes
 *	appear in element checks. In this case, we build the matrix structure with no data.
 *	
 *	This function works for both double and dcomplex data.
 */

SCSCMatrix :: SCSCMatrix (const char* filename) :
	compressedMatrix(filename,true)
{
	storage = SCSC;
}





SCSCMatrix :: SCSCMatrix(model &m, char datatype)
{
	createProfile(m);
}


// Dealllocate the data allocated in NewSCSC but also the L U factors that might have
// been created
SCSCMatrix :: ~SCSCMatrix()
{
}




// create a new SCSC matrix, with the same profile as a model one
SCSCMatrix :: SCSCMatrix(const SCSCMatrix &m)
{
	// get information about the dimensions of the model matrix m
	nnz  = m.nnz;
	cols = m.cols;
	
	// allocate the space for the 4 vectors of the new matrix
	// if data is complex, allocate twice the space
	int fdouble  = (datatype == DATA_DOUBLE) ? 1 : 2;
	ddata    = new double[nnz*fdouble]();
	colstart = new int[cols+1];
	prow     = new int[nnz];
	pdiag    = new int[cols];
	if (ddata == NULL || colstart == NULL || prow == NULL || pdiag == NULL)
		printf("Error in SCSCDuplicate. No memory for data.");
	
	
	// copy the int & double data
	Icopy(m.prow, prow, nnz);
	Icopy(m.colstart, colstart, cols+1);
	Icopy(m.pdiag, pdiag, cols);
	Dcopy(m.ddata, ddata, nnz);
	
	
	// allocate the whole structure and connect pointers
	datatype    = m.datatype;
	numbering   = m.numbering;
	rows        = m.rows;
}




// create an identical SCSC matrix, not factorized
SCSCMatrix& SCSCMatrix :: operator=(const SCSCMatrix &m)
{    
	// get information about the dimensions of the model matrix m
	nnz  = m.nnz;
	cols = m.cols;
	
	// allocate the space for the 4 vectors of the new matrix
	// if data is complex, allocate twice the space
	int fdouble  = (datatype == DATA_DOUBLE) ? 1 : 2;
	ddata    = new double[nnz*fdouble]();
	colstart = new int[cols+1];
	prow     = new int[nnz];
	pdiag    = new int[cols];
	if (ddata == NULL || colstart == NULL || prow == NULL || pdiag == NULL)
		printf("Error in SCSCDuplicate. No memory for data.");
	
	
	// copy the int & double data
	Icopy(m.prow, prow, nnz);
	Icopy(m.colstart, colstart, cols+1);
	Icopy(m.pdiag, pdiag, cols);
	Dcopy(m.ddata, ddata, nnz);
	
	// allocate the whole structure and connect pointers
	datatype    = m.datatype;
	numbering   = m.numbering;
	rows        = m.rows;
	
	return *this;
}




/* add one element in the sparse matrix m given its row and column, to the existing data,
 in the standard way. This is only for matrices with SCSC storage */

void SCSCMatrix ::  addTo(double x, int row, int col)
{
	//quick return
	if (col > row) return;
	
    if (false)
    {
        int  k, kl;
	
        k  = colstart[col];		// pointer to the first data of column col
        kl = colstart[col+1];	// pointer to the first data of following column
	
        while( k < kl)
        {
            if ( prow[k] == row ) 
            {
                ddata[k] += x; 
                break;
            }
            k++;
        }
	
        if (k == kl) printf("\nCannot add data in position row:%d  col:%d", row, col);
    }
    
    else
    {
        int s  = colstart[col];
        int e  = colstart[col+1]-1;
        int rs = prow[s];
        int re = prow[e];
        
        if (rs == row)
        {
            ddata[s] += x;
            return;
        }
        
        if (re == row)
        {
            ddata[e] += x;
            return;
        }
        
        
        while (true)
        {
            int jump = ((row-rs)*(e-s))/(re-rs);
            jump = (jump == 0) ? 1 : jump;
            jump = (jump == e-s) ? e-s-1 : jump;
            int i  = s + jump;
            int ri = prow[i];
            
            if (ri < row)
            {
                s  = i;
                rs = ri;
            }
            
            else if (ri > row)
            {
                e  = i;
                re = ri;
            }
            
            else 
            {
                ddata[i] += x;                
                break;
            }
        }
    }
}






void SCSCMatrix :: assemble(const matrix& block, const std::vector<int>& position) 
{
	if      (datatype == DATA_DOUBLE)   assembleDouble (block, position);
	// else if (datatype == DATA_DCOMPLEX) assembleComplex(block, position);
	else printf("Can not assemble. Datatype unknown");
}




/*  assemble matrix 'block', that has full format in matrix 'big', which is sparse.
 The idata list gives the equation number where the entries go as:
 big( list(i), list(j) ) += block(i,j) 
 Note: the equations numbered with negative numbers need not be assembled
 Note: we assume that in the SCSC format, the nonzero elements in a column are stored
 with growing row numbers.
 */ 

void SCSCMatrix :: assembleDouble(const matrix& block, const std::vector<int>& position) 
{
	int     r, m,c_start, c_end;
	int     globalcol, globalrow, blkdim;
	double  **BLOCK=NULL;
	
	
	if (cols == 0) return;
	
	blkdim    = block.rows;
	BLOCK     = block.data;
	
	for (int k=0; k<blkdim; k++)					// loop over the columns of the matrix 'block'
	{
		globalcol = position[k]; 
		
		if (globalcol > -1)                         // continue only if the equation needs to be assembled
		{
			c_start   = colstart[globalcol];    	// position in sparsedata[] of the first data of column col
			c_end     = colstart[globalcol+1];  	// position in sparsedata[] to the first data of next column
			for (m=0; m<blkdim; m++)                // loop over the rows of the matrix 'block'
			{
				globalrow = position[m]; 
				
				if (globalrow >= globalcol)         // continue only if the equation needs to be assembled
				{
					r = c_start;                    // the row that is being scanned
					while(r < c_end)
					{
						if (prow[r] == globalrow)
						{
							ddata[r] += BLOCK[m][k];
							break;
						}
						r++;
					}
					if (r == c_end) printf("\nCannot add data in row:%d  col:%d", globalrow, globalcol);
				}
			}
		}
	}
} 




/* assemble matrix 'block', that has full format in matrix 'big', which is sparse.
 The idata list gives the equation number where the entries go as:
 big( list(i), list(j) ) += block(i,j) 
 Note: the equations numbered with negative numbers need not be assembled
 Note: we assume that in the SCSC format, the nonzero elements in a column are stored
 with growing row numbers.
 
 void SCSCMatrix :: assembleComplex(matrix block, idata position) 
 {
 int     r, m, c_start, c_end;
 int     globalcol, globalrow, blkdim;
 complex double *sparsedata=NULL, **BLOCK=NULL;
 
 // recover sparseMatrix structure for SCSC storage
 sparsedata= (complex double *) ddata;
 
 if (cols == 0) return;
 
 blkdim    = block->rows;
 BLOCK     = (complex double **)block->data;
 
 for (int k=0; k<blkdim; k++)                        // loop over the columns of the matrix 'block'
 {
 globalcol = position->data[k]; 
 
 if (globalcol > -1)                         // continue only if the equation needs to be assembled
 {
 c_start   = colstart[globalcol];    	// position in sparsedata[] of the first data of column col
 c_end     = colstart[globalcol+1];  	// position in sparsedata[] to the first data of next column
 for (m=0; m<blkdim; m++)                // loop over the rows of the matrix 'block'
 {
 globalrow = position->data[m]; 
 
 if (globalrow > -1)                 // continue only if the equation needs to be assembled
 {
 r = c_start;                    // the row that is being scanned
 while(r < c_end)
 {
 if (prow[r]  == globalrow)
 {
 sparsedata[r] += BLOCK[m][k];
 break;
 }
 r++;
 }
 if (r == c_end) printf("\nCannot add data in row:%d  col:%d", globalrow, globalcol);
 }
 }
 }
 }
 } 
 */


/* we plan to put several tests in here to confirm that the basic data of the SCSC matrix
 is ok */


void SCSCMatrix :: check()
{
	// check that at least every column has 1 element
	for (int k=1; k<cols; k++)
		if (colstart[k-1] >= colstart[k]) printf("Bad sparse matrix check. Empty column");
	
	// check that the numbering is either fortran or c
	if (numbering != NUM_FORTRAN && numbering != NUM_C)
		printf("Bad sparse matrix check. Numbering style unknown.");
	
	// check that the datatype is either double of complex double
	if (datatype != DATA_DOUBLE && datatype != DATA_DCOMPLEX)
		printf("Bad sparse matrix check. Datatype unknown.");
}




void SCSCMatrix :: createProfile(const model &m)
{
	// quick return
	int ndof = m.getNDofs();
	if (ndof <= 0) return;
	
	// gather the positions of non-zero entries on the matrix
	// a colection of vectors, one per column
	vector<int> *global = new vector<int>[ndof]();
	
	// now go element by element and fill up nonzero data
    vector<int> theDofs;
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(element* ele, bd->elements)
        {
			ele->getDofNumbers(theDofs);
			int   neldofs = theDofs.size();
			
			for (int i=0; i<neldofs; i++)
			{
				int col = theDofs[i];             // retrieve degree of freedom
				
				if (col >= 0)							// modify if dof is active, ie, nonnegative
				{
					for (int k=0; k<neldofs; k++)
					{
						int row = theDofs[k];;
						
						// modify profile if dof is active. Only lower part
						if (row >= col) global[col].push_back(row);
					}
				}
			}
		}
        
        BOOST_FOREACH(evalspot* sp, bd->evalspots)
        {
			sp->getDofNumbers(theDofs);
			int   neldofs = theDofs.size();
			
			for (int i=0; i<neldofs; i++)
			{
				int col = theDofs[i];             // retrieve degree of freedom
				
				if (col >= 0)                       // modify if dof is active, ie, nonnegative
				{
					for (int k=0; k<neldofs; k++)
					{
						int row = theDofs[k];
						if (row >= 0) global[col].push_back(row);
					}
				}
			}
		}
	}
	
	
	// go column by column, ordering them and trimming repeated values
	nnz = 0;
	for (int i=0; i<ndof; i++)
	{
		sort( global[i].begin(), global[i].end() );
		vector<int>::iterator end_unique = unique( global[i].begin(), global[i].end() );
		global[i].erase( end_unique, global[i].end() );
		nnz += global[i].size();
	}
	
	
	// go entry by entry filling up the Harwell-Boeing pointers
	colstart = new int[ndof+1];
	prow     = new int[nnz];
	ddata    = new double[nnz]();
	pdiag    = new int[ndof];
	
	int counter = 0;
	for (int i=0; i<ndof; i++)
	{
		colstart[i] = pdiag[i] = counter;
		for (int j=0; j<global[i].size(); j++) prow[counter++] = global[i][j];
	}
	colstart[ndof] = counter;
	
	
	// erase the temporary vectors
	for (int i=0; i<ndof; i++) global[i].clear();
	delete [] global;
	
	//  place everything in a matrix structure
	rows    = ndof;
	cols    = ndof;
	storage = SCSC;
}




void SCSCMatrix :: dump(string format)
{  
	if (format == "matrixmarket" || format == "mm")
	{
		ofstream os;
		os.open("feliks.tangent.mm");
		
		os << "%%MatrixMarket matrix coordinate real general" << "\n";
		os << "%% standard establishes 1 based " << "\n";
		os << "%% Nrows = " << rows << ", Ncols = " << cols << ", nnz = " << nnz << "\n";
		os << rows << " " << cols << " " << 2*nnz-cols << "\n";
		
		int cc = 0;
		for (int i=0; i<nnz; i++)
		{
			if (i == colstart[cc+1]) cc++;
			os  << setw(5)  << prow[i]+1 << " " 
			<< setw(5)  << cc+1 << " " 
			<< setw(18) << setprecision(14) << ddata[i] << "\n";
			if (prow[i] != cc) 
				os  << setw(5)  << cc+1 << " " 
				<< setw(5)  << prow[i]+1 << " " 
				<< setw(18) << setprecision(14) << ddata[i] << "\n";
		}
		
		os.close();
	}
	
	else if (format == "harwellboeing" || format == "harwell-boeing" || format == "hb")
	{
		/* This uses the Matrix Market subroutines to write a hwb sparse matrix
		 to a file. It produces a "nicer" output that the one above */
		
		int  space;
		char filename[]="feliks.tangent.hwb";
		char title[]="feliks tangent";
		char key[]="feliks key";
		char type[]="RSA";
		char Ptrfmt[]="(13I7)";
		char Indfmt[]="(16I7)";
		char Valfmt[]="(3E26.18)";
		char Rhsfmt[]="(3E26.18)";
		
		// get the size for the row and column pointers and leave one space for readability
		space = (int) log10((double)nnz) + 2;
		Ptrfmt[4] = space + toascii('0');
		Indfmt[4] = space + toascii('0');
		
		writeHB_mat_double(filename, rows, cols, nnz, colstart, prow,
						   ddata, 0, NULL, NULL, NULL,
						   title, key, type, 
						   Ptrfmt, Indfmt, Valfmt, Rhsfmt,
						   NULL);
	}
}



double SCSCMatrix :: getEntry(int row, int col)
{
	int k, kl;
	double x=0.0;
	
	if (row >= col)
	{
		
		// recover sparseMatrix structure for SCSC storage
		k  = colstart[col];			// pointer to the first data of column col
		kl = colstart[col+1];		// pointer to the first data of following column
		
		while( k < kl)
		{
			if ( prow[k] == row )
			{
				x = ddata[k];
				break;
			}
			k++;
		}
		
		if (k == kl) printf("Can not extract data from position row:%d  col:%d", row, col);
	}
	else
		x = getEntry(col,row);
	
	return x;
}



void  SCSCMatrix :: info(ostream &of)
{
	sparseMatrix::info(of);
	
	of	<< "\n\tFormat                     : symmetric Compressed Sparse Column\n";
}





void SCSCMatrix :: insert(double x, int row, int col)
{
	int k, kl;
	
	// Make sure we work in C style
	convertToC();
	
	k  = colstart[col];    /* pointer to the first data of column col */
	
	kl = colstart[col+1];  /* pointer to the first data of following column */
	
	
	while( k < kl)
	{
		if ( prow[k] == row )
		{
			ddata[k] = x;
			break;
		}
		k++;
	}
	
	if (k == kl)
	{
		printf("Can not insert data in position row:%d  col:%d", row, col);
	}
}




void  SCSCMatrix :: print(FILE *fp)
{
	matrix  full;
	int     k;
	
	fprintf(fp, "\n Sparse Matrix, symmetric CSC format");
	fprintf(fp, "\n  >Number of data    = %3d", nnz);
	fprintf(fp, "\n  >Number of columns = %3d", cols);
	fprintf(fp, "\n  >Number of rows    = %3d", rows);
	if (datatype == DATA_DOUBLE)
		fprintf(fp, "\n  >Type of data      = DOUBLE");
	else if (datatype == DATA_COMPLEX)
		fprintf(fp, "\n  >Type of data      = COMPLEX");
	fprintf(fp, "\n  >Columns start in data number :");
	fprintf(fp, "\n     "); for(k=0; k<cols+1; k++) fprintf(fp, " %2d", colstart[k]);
	fprintf(fp, "\n >Row pointers :");
	fprintf(fp, "\n     "); for(k=0; k<nnz; k++) fprintf(fp, " %2d", prow[k]);
	fprintf(fp, "\n >Diagonal pointers :");
	fprintf(fp, "\n     "); for(k=0; k<cols; k++) fprintf(fp, " %2d", pdiag[k]);
	
	full = toFull();
	full.print();
	
	longvector  diag( full.diagonal() );
	fprintf(fp, "\n\n>Diagonal:");
	diag.print();
	fprintf(fp, "\n\n");
	
	fflush(stdout);
}



void SCSCMatrix :: timesDoubleArray(double* b, double* Axb)
{
	int               startc, endc, row;
	double            bdatac, Adatak;
	register int      c, k;
	register double   tmp;
	
	// recover sparseMatrix structure for SCSC storage
	int     n       = cols;
	double *Adata   = ddata;
	
	
	Dzero(Axb, n);
	for (c=0; c<n; c++) 
	{
		// pointers to the start and end data of the c-th column below the diagonal
		startc = colstart[c];
		endc   = colstart[c+1];
		
		// multiply diagonal element by x
		bdatac      = b[c];
		Axb[c] += Adata[startc]*bdatac; 
		
		tmp = 0.0;
		// compute inner product of row i with vector x
		for (k=startc+1; k<endc; k++)
		{
			row    = prow[k];
			Adatak = Adata[k]; 
			
			// upper part
			tmp += Adatak*b[row];
			
			// lower part
			Axb[row] += Adatak* bdatac;
		}
		Axb[c] += tmp;
	}
}


void SCSCMatrix :: timesVectorDouble(longvector& b, longvector& Axb)
{
	int               startc, endc, row;
	double            bdatac, Adatak;
	register int      c, k;
	register double   tmp;
	
	// recover sparseMatrix structure for SCSC storage
	int     n       = cols;
	double *Adata   = ddata;
	double *bdata   = b.data;
	double *Axbdata = Axb.data;
	
	Dzero(Axbdata, n);
	for (c=0; c<n; c++) 
	{
		// pointers to the start and end data of the c-th column below the diagonal
		startc = colstart[c];
		endc   = colstart[c+1];
		
		// multiply diagonal element by x
		bdatac      = bdata[c];
		Axbdata[c] += Adata[startc]*bdatac; 
		
		tmp = 0.0;
		// compute inner product of row i with vector x
		for (k=startc+1; k<endc; k++)
		{
			row    = prow[k];
			Adatak = Adata[k]; 
			
			// upper part
			tmp += Adatak*bdata[row];
			
			// lower part
			Axbdata[row] += Adatak* bdatac;
		}
		Axbdata[c] += tmp;
	}	
}




/*
 void  SCSCMatrix :: timesVectorDComplex(longvector b, longvector Axb)
 {
 int                startc, endc, row;
 double             bdatac, Adatak;
 register int       c, k;
 complex double   tmp;
 
 // recover sparseMatrix structure for SCSC storage
 int n         = cols;
 complex double *Adata   = (complex double *)ddata;
 complex double *bdata   = (complex double *)b->data;
 complex double *Axbdata = (complex double *)Axb->data;
 
 
 // Full multiplication 
 if (mathtype == GEN)
 {
 for(c=0; c<n; c++) 
 {
 // Standard implementation
 // pointers to the start and end data of the cth row (or column)
 startc = colstart[c];
 endc   = colstart[c+1];
 
 // compute inner product of row c with vector x
 tmp = 0.0 + 0.0*I;
 for (k=startc; k<endc; k++) tmp += Adata[k] * bdata[prow[k]];
 Axbdata[c] = tmp; 
 }
 }
 
 else if (mathtype == SYMM)
 {
 // first multiply the diagonal of the matrix A by b 
 for (c=0; c<n; c++) Axbdata[c] = Adata[pdiag[c]]*bdata[c];
 
 for(c=0; c<n; c++) 
 {
 // pointers to the start and end data of the ith column above the diagonal
 startc = colstart[c];
 endc   = pdiag[c];
 bdatac = bdata[c];
 tmp    = 0.0;
 
 // compute inner product of row i with vector x
 for (k=startc; k<endc; k++)
 {
 row = prow[k];
 
 // lower part
 Adatak = Adata[k]; 
 tmp   += Adatak*bdata[row];
 
 // upper part
 Axbdata[row] += Adatak* bdatac;
 }
 Axbdata[c] += tmp;
 }
 }
 }
 */




// creates a full matrix by filling up entries from sparse matrix
matrix SCSCMatrix :: toFull()
{
	int     k   ,  j;
	
	matrix full(rows, cols);
	if (datatype == DATA_DOUBLE)
	{
		for (k=0; k< cols; k++)
		{
			full.data[k][k] = ddata[colstart[k]];
			for (j=colstart[k]+1; j<colstart[k+1]; j++)
				full.data[prow[j]][k] = full.data[k][prow[j]] = ddata[j];
		}
	}
	
	else
	{
	}		
	
	return full;
}




CSCMatrix :: CSCMatrix() :
	compressedMatrix()
{
}



CSCMatrix :: CSCMatrix (const char* filename) :
	compressedMatrix(filename,false)
{
	storage = CSC;
}



/*	function to allocate the memory space for a matrix in the Harwel-Boeing (aka compress sparse column) format.
 *	We might consider matrices for models with 0 dofs. These are not common but sometimes
 *	appear in element checks. In this case, we build the matrix structure with no data.
 *	
 *	This function works for both double and dcomplex data.
 */

CSCMatrix :: CSCMatrix(model &m, char datatype)
{
	createProfile(m);
}



void CSCMatrix :: createProfile(const model &m)
{
    // quick return
	int ndof = m.getNDofs();
	if (ndof <= 0) return;
	
	// gather the positions of non-zero entries on the matrix
	// a colection of vectors, one per column
	vector< set<int> > global; global.resize(ndof);
	vector<int> theDofs;
    
	// now go element by element and fill up nonzero data
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(element* ele, bd->elements)
        {
			ele->getDofNumbers(theDofs);
			int   neldofs = theDofs.size();
			
			for (int i=0; i<neldofs; i++)
			{
				int col = theDofs[i];             // retrieve degree of freedom
				
				if (col >= 0)                       // modify if dof is active, ie, nonnegative
				{
					for (int k=0; k<neldofs; k++)
					{
						int row = theDofs[k];
						if (row >= 0) global[col].insert(row);
					}
				}
			}
		}

        BOOST_FOREACH(evalspot* sp, bd->evalspots )
        {
			sp->getDofNumbers(theDofs);
			int   neldofs = theDofs.size();
			
			for (int i=0; i<neldofs; i++)
			{
				int col = theDofs[i];             // retrieve degree of freedom
				
				if (col >= 0)                       // modify if dof is active, ie, nonnegative
				{
					for (int k=0; k<neldofs; k++)
					{
						int row = theDofs[k];
						if (row >= 0) global[col].insert(row);
					}
				}
			}
		}
	}
    
	
	
	// go column by column, ordering them and trimming repeated values
	nnz = 0;
	for (int i=0; i<ndof; i++)
	{
		nnz += global[i].size();
	}
	
	
	// go entry by entry filling up the Harwell-Boeing pointers
	colstart = new int[ndof+1];
	prow     = new int[nnz];
	ddata    = new double[nnz]();
	pdiag    = new int[ndof];
	
	int counter = 0;
	for (int i=0; i<ndof; i++)
	{
		colstart[i] = counter;
        BOOST_FOREACH(int j, global[i])
            prow[counter++] = j;
//		for (int j=0; j<global[i].size(); j++) prow[counter++] = global[i][j];
	}
	colstart[ndof] = counter;
	
	
	// fill diagonal counter
	for (int i=0; i<ndof; i++)
	{
		int j;
		for (j = colstart[i]; j < colstart[i+1]; j++)
			if ( prow[j] == i)
			{
				pdiag[i] = j;
				break;
			}
        
		if (j == colstart[i+1])
			printf("\n Error in computing diagonal. Entry is empty");
	}
	
	
	//  place everything in a matrix structure
	rows    = ndof;
	cols    = ndof;
	storage = CSC;
}




// create a new CSC matrix, with the same profile as a model one, but with entries equal to 0
CSCMatrix :: CSCMatrix(const CSCMatrix &m) 
{
	// get information about the dimensions of the model matrix m
	nnz  = m.nnz;
	cols = m.cols;
	
	// allocate the space for the 4 vectors of the new matrix
	// if data is complex, allocate twice the space
	int fdouble  = (datatype == DATA_DOUBLE) ? 1 : 2;
	ddata    = new double[nnz*fdouble]();
	colstart = new int[cols+1];
	prow     = new int[nnz];
	pdiag    = new int[cols];
	if (ddata == NULL || colstart == NULL || prow == NULL || pdiag == NULL)
		printf("Error in CSCDuplicate. No memory for data.");
	
	
	// copy the int & double data
	Icopy(m.prow, prow, nnz);
	Icopy(m.colstart, colstart, cols+1);
	Icopy(m.pdiag, pdiag, cols);
	
	
	// allocate the whole structure and connect pointers
	storage     = CSC;
	datatype    = m.datatype;
	numbering   = m.numbering;
	rows        = m.rows;
}


CSCMatrix :: CSCMatrix(const SCSCMatrix &m) :
	compressedMatrix()
{
	// get information about the dimensions of the model matrix m
	nnz  = m.getNProfileTerms();
	cols = m.getNColumns();
	
	// allocate the space for the 4 vectors of the new matrix
	// if data is complex, allocate twice the space
	int fdouble  = (datatype == DATA_DOUBLE) ? 1 : 2;
	ddata    = new double[nnz*fdouble]();
	colstart = new int[cols+1];
	prow     = new int[nnz];
	pdiag    = new int[cols];
	if (ddata == NULL || colstart == NULL || prow == NULL || pdiag == NULL)
		printf("Error in CSCDuplicate. No memory for data.");
	
	
	// copy the int & double data
	Icopy(m.getProw(), prow, nnz);
	Icopy(m.getColstart(), colstart, cols+1);
	Icopy(m.getPdiag(), pdiag, cols);
	
	
	// allocate the whole structure and connect pointers
	storage     = SCSC;
	rows        = m.getNRows();
	
	expand_symmetric_storage();
}


// Dealllocate the data allocated in NewCSC but also the L U factors that might have
// been created in the SuperLU solver
CSCMatrix :: ~CSCMatrix()
{
}





// create an identical CSC matrix, not factorized
CSCMatrix& CSCMatrix :: operator=(const CSCMatrix &m)
{    
	// get information about the dimensions of the model matrix m
	nnz  = m.nnz;
	cols = m.cols;
	
	// allocate the space for the 4 vectors of the new matrix
	// if data is complex, allocate twice the space
	int fdouble  = (datatype == DATA_DOUBLE) ? 1 : 2;
	ddata    = new double[nnz*fdouble]();
	colstart = new int[cols+1];
	prow     = new int[nnz];
	pdiag    = new int[cols];
	if (ddata == NULL || colstart == NULL || prow == NULL || pdiag == NULL)
		printf("Error in CSCDuplicate. No memory for data.");
	
	
	// copy the int & double data
	Icopy(m.prow, prow, nnz);
	Icopy(m.colstart, colstart, cols+1);
	Icopy(m.pdiag, pdiag, cols);
	Dcopy(m.ddata, ddata, nnz*fdouble);
	
	// allocate the whole structure and connect pointers
	storage     = CSC;
	datatype    = m.datatype;
	numbering   = m.numbering;
	rows        = m.rows;
	factorized  = -1;
	condition   = 0.0;	
	return *this;
}




/* add one element in the sparse matrix m given its row and column, to the existing data,
 in the standard way. This is only for matrices with CSC storage */

void CSCMatrix ::  addTo(double x, int row, int col)
{
    // this is the old, brute force assember. Leave it for the 
    // future, but deactivate it
    if (false)
    {
        int k  = colstart[col];		// pointer to the first data of column col
        int kl = colstart[col+1];	// pointer to the first data of following column
	
        while( k < kl)
        {
            if ( prow[k] == row ) 
            {
                ddata[k] += x; 
                break;
            }
            k++;
        }
	
        if (k == kl) printf("\nCannot add data in position row:%d  col:%d", row, col);
    }
    
    else
    {
        int s  = colstart[col];
        int e  = colstart[col+1]-1;
        int rs = prow[s];
        int re = prow[e];

        if (rs == row)
        {
            ddata[s] += x;
            return;
        }
        
        if (re == row)
        {
            ddata[e] += x;
            return;
        }

        
        while (true)
        {
            int jump = ((row-rs)*(e-s))/(re-rs);
            jump = (jump == 0) ? 1 : jump;
            jump = (jump == e-s) ? e-s-1 : jump;
            int i  = s + jump;
            int ri = prow[i];
            
            if (ri < row)
            {
                s  = i;
                rs = ri;
            }
            
            else if (ri > row)
            {
                e  = i;
                re = ri;
            }
            
            else 
            {
                ddata[i] += x;                
                break;
            }
        }
    }
    
}




void CSCMatrix :: assemble(const matrix& block, const std::vector<int>& position) 
{
	if      (datatype == DATA_DOUBLE)   assembleDouble (block, position);
	else if (datatype == DATA_DCOMPLEX) assembleComplex(block, position);
	else printf("Can not assemble. Datatype unknown");
}




/*  assemble matrix 'block', that has full format in matrix 'big', which is sparse.
 The idata list gives the equation number where the entries go as:
 big( list(i), list(j) ) += block(i,j) 
 Note: the equations numbered with negative numbers need not be assembled
 Note: we assume that in the CSC format, the nonzero elements in a column are stored
 with growing row numbers.
 */ 

void CSCMatrix :: assembleDouble(const matrix &block, const std::vector<int>& position) 
{
	int     r, m,c_start, c_end;
	int     globalcol, globalrow, blkdim;
	double  **BLOCK=NULL;
	
	if (cols == 0) return;
	
	blkdim = block.rows;
	BLOCK  = block.data;
	
	for (int k=0; k<blkdim; k++)                    // loop over the columns of the matrix 'block'
	{
		globalcol = position[k]; 
		
		if (globalcol > -1)                         // continue only if the equation needs to be assembled
		{
			c_start   = colstart[globalcol];    	// position in sparsedata[] of the first data of column col
			c_end     = colstart[globalcol+1];  	// position in sparsedata[] to the first data of next column
			for (m=0; m<blkdim; m++)                // loop over the rows of the matrix 'block'
			{
				globalrow = position[m]; 
				
				if (globalrow > -1)                 // continue only if the equation needs to be assembled
				{
					r = c_start;                    // the row that is being scanned
					while(r < c_end)
					{
						if (prow[r] == globalrow)
						{
							ddata[r] += BLOCK[m][k];
							break;
						}
						r++;
					}
					if (r == c_end) 
						printf("\nCannot add data in row:%d  col:%d", globalrow, globalcol);
				}
			}
		}
	}
} 




/* assemble matrix 'block', that has full format in matrix 'big', which is sparse.
 The idata list gives the equation number where the entries go as:
 big( list(i), list(j) ) += block(i,j) 
 Note: the equations numbered with negative numbers need not be assembled
 Note: we assume that in the CSC format, the nonzero elements in a column are stored
 with growing row numbers.
 */
void CSCMatrix :: assembleComplex(const matrix &block, const std::vector<int>& position) 
{
} 



/* we plan to put several tests in here to confirm that the basic data of the CSC matrix
 is ok */
void CSCMatrix :: check()
{
	int k;
	
	/* check that at least every column has 1 element */
	for (k=1; k<cols; k++)
		if (colstart[k-1] >= colstart[k]) printf("Bad sparse matrix check. Empty column");
	
	// check that the numbering is either fortran or c
	if (numbering != NUM_FORTRAN && numbering != NUM_C)
		printf("Bad sparse matrix check. Numbering style unknown.");
	
	// check that the datatype is either double of complex double
	if (datatype != DATA_DOUBLE && datatype != DATA_DCOMPLEX)
		printf("Bad sparse matrix check. Datatype unknown.");
}





double CSCMatrix :: getEntry(int row, int col)
{
	int k, kl;
	double x=0.0;
	
	/* recover sparseMatrix structure for CSC storage */	
	k  = colstart[col];			/* pointer to the first data of column col */
	
	kl = colstart[col+1];		/* pointer to the first data of following column */
	
	while( k < kl)
	{
		if ( prow[k] == row )
		{
			x = ddata[k];
			break;
		}
		k++;
	}
	
	if (k == kl) printf("Can not extract data from position row:%d  col:%d", row, col);
	
	return x;
}



void  CSCMatrix :: info(ostream &of)
{
	sparseMatrix::info(of);
	of	<< "\n\tFormat                     : unsymmetric Compressed Sparse Column";
}



void CSCMatrix :: insert(double x, int row, int col)
{
	int k, kl;
	
	// Make sure we work in C style
	
	convertToC();
	
	k  = colstart[col];    /* pointer to the first data of column col */
	
	kl = colstart[col+1];  /* pointer to the first data of following column */
	
	while( k < kl)
	{
		if ( prow[k] == row )
		{
			ddata[k] = x;
			break;
		}
		k++;
	}
	
	if (k == kl)
	{
		printf("Can not insert data in position row:%d  col:%d", row, col);
	}
}



void CSCMatrix :: timesVectorDouble(longvector& b, longvector& Axb)
{
	int               startc, endc;
	register int       c, k;
	register double   tmp;
	
	// recover sparseMatrix structure for CSC storage
	int n           = cols;
	double *Adata   = ddata;
	double *bdata   = b.data;
	double *Axbdata = Axb.data;
	
	for(c=0; c<n; c++) 
	{
		// Standard implementation
		// pointers to the start and end data of the cth column
		startc = colstart[c];
		endc   = colstart[c+1];
		
		// compute inner product of row c with vector x
		tmp = 0.0;
		for (k=startc; k<endc; k++) tmp += Adata[k] * bdata[prow[k]];
		Axbdata[c] = tmp; 
	}
}




/****************************************************************************************/
/*************************** A√ëADIDO MADE IN BORJA***************************************/
/****************************************************************************************/

CSCMatrix :: CSCMatrix(int size, int numbernnz, double *nzelements, int *rowvect, int *colvect)
{
	rows     = size;
	cols     = size;
	nnz 	 = numbernnz;
	ddata 	 = new double[nnz];
	prow	 = new int[nnz];
	colstart = new int[size+1];
	for(int i=0; i<nnz; i++)
	{
		ddata[i] = nzelements[i];
		prow[i]  = rowvect[i];
	}
	for(int i=0; i<=size; i++) colstart [i] = colvect[i];
}






void CSCMatrix ::	timesDoubleArray(double* b, double* Axb)
{
	int       startc, endc;
	int       c, k;
	double   tmp;
	
	// recover sparseMatrix structure for CSC storage
	int n           = cols;
	double *Adata   = ddata;
	
	for(c=0; c<n; c++) 
	{
		// Standard implementation
		// pointers to the start and end data of the cth column
		startc = colstart[c];
		endc   = colstart[c+1];
		
		// compute inner product of row c with vector x
		tmp = 0.0;
		for (k=startc; k<endc; k++) tmp += Adata[k] * b[prow[k]];
		Axb[c] = tmp; 
	}
}





/********************************************************************************************/




/*
 void  CSCMatrix :: timesVectorDComplex(longvector b, longvector Axb)
 {
 int                startc, endc, row;
 double             bdatac, Adatak;
 register int       c, k;
 complex double   tmp;
 
 // recover sparseMatrix structure for CSC storage
 int n         = cols;
 complex double *Adata   = (complex double *)ddata;
 complex double *bdata   = (complex double *)b->data;
 complex double *Axbdata = (complex double *)Axb->data;
 
 
 // Full multiplication 
 if (mathtype == GEN)
 {
 for(c=0; c<n; c++) 
 {
 // Standard implementation
 // pointers to the start and end data of the cth row (or column)
 startc = colstart[c];
 endc   = colstart[c+1];
 
 // compute inner product of row c with vector x
 tmp = 0.0 + 0.0*I;
 for (k=startc; k<endc; k++) tmp += Adata[k] * bdata[prow[k]];
 Axbdata[c] = tmp; 
 }
 }
 
 else if (mathtype == SYMM)
 {
 // first multiply the diagonal of the matrix A by b 
 for (c=0; c<n; c++) Axbdata[c] = Adata[pdiag[c]]*bdata[c];
 
 for(c=0; c<n; c++) 
 {
 // pointers to the start and end data of the ith column above the diagonal
 startc = colstart[c];
 endc   = pdiag[c];
 bdatac = bdata[c];
 tmp    = 0.0;
 
 // compute inner product of row i with vector x
 for (k=startc; k<endc; k++)
 {
 row = prow[k];
 
 // lower part
 Adatak = Adata[k]; 
 tmp   += Adatak*bdata[row];
 
 // upper part
 Axbdata[row] += Adatak* bdatac;
 }
 Axbdata[c] += tmp;
 }
 }
 }
 */




void  CSCMatrix :: print(FILE *fp)
{
	int     k;

	fprintf(fp, "\n\n\n");
	fprintf(fp, "\n CSC Sparse Matrix");
	fprintf(fp, "\n  >Number of columns = %3d", cols);
	fprintf(fp, "\n  >Number of rows    = %3d", rows);
	fprintf(fp, "\n  >Number of data    = %3d", nnz);
	fprintf(fp, "\n  >Type of data      = %3d", datatype);
	fprintf(fp, "\n  >Columns start in data number :");
	fprintf(fp, "\n     "); for(k=0; k<cols+1; k++) fprintf(fp, " %2d", colstart[k]);
	fprintf(fp, "\n >Row pointers :");
	fprintf(fp, "\n     "); for(k=0; k<nnz; k++) fprintf(fp, " %2d", prow[k]);
	fprintf(fp, "\n >Diagonal pointers :");
	fprintf(fp, "\n     "); for(k=0; k<cols; k++) fprintf(fp, " %2d", pdiag[k]);
	
	matrix  full;
	full = toFull();
	full.print();
	
	fprintf(fp, "\n\n>Diagonal:");
	longvector  diag(cols);
	getDiagonal(diag);
	diag.print();
	fprintf(fp, "\n\n");
	
	
	fflush(stdout);
}


// creates a full matrix by filling up entries from sparse matrix
matrix CSCMatrix :: toFull()
{
	matrix full(rows,cols);
	if (datatype == DATA_DOUBLE)
	{
		for (int k=0; k< cols; k++)
			for (int j=colstart[k]; j<colstart[k+1]; j++)
				full.data[prow[j]][k] = ddata[j];
	}
	
	else
	{
	}		
	
	return full;
}


