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
    /* matrix.h
 * matrices
 * Ignacio Romero
 * July 1999  
 *
 * defines  (full) matrix and functions to deal with them. 
 *
 * Despite the fact feliks is written with some "object oriented" philosophy,
 * the full matrix and vector data types is not treated in this way throughout
 * the program. This is because data and vectors are accessed all the time
 * and it is faster to just use the data as m->data[i][j] or things like this.
 * This is different than the sparse structure, which should never be
 * accessed directly
 *
 * this set of functions are only for FULL MATRICES!!!
 */


#ifndef _matrix_h
#define _matrix_h


#include "Math/vector.h"
#include "General/idata.h"




typedef enum { 
	DATA_SINGLE,
	DATA_DOUBLE,
	DATA_COMPLEX,
	DATA_DCOMPLEX
} dataT;


class matrix
	{
	private:
		dataT   datatype;				// double, complex, etc. see vector.h
		
		

	public:
		int     rows;					// number of rows for full matrices
		int     cols;					// number of cols for full matrices
		double **data;
		double *start;					// where the data goes
		
		matrix();
		matrix(const size_t rows, const size_t cols);
		matrix(const matrix& m);
		matrix(const longvector& v);
		~matrix();
		
		matrix&						operator=(const matrix& m);
		inline			double&		operator()(const int i, const int j)		{return data[i][j];}
		inline const	double&		operator()(const int i, const int j) const	{return data[i][j];}
		matrix&						operator+=(const matrix &v);
		matrix&						operator-=(const matrix &v);


		static matrix  IdentityMatrix(int n);
		void   resize(const size_t newrows, const size_t newcols);



/* functions to put elements or blocks in the matrix */
void            addBlock(matrix& block, const int row , const int col);
inline	void	addTo				(const double f, int row, int col){ data[row][col] += f;}
void            AddBlockToMatrix    (matrix bl  , matrix    m , int row , int    col);
void            CopyMatrix          (matrix m   , matrix  mcopy);                   
void            InsertElement       (matrix m   , double    e , int row , int    col);
static void     InsertBlockInMatrix (matrix &bl, matrix &m , const int row, const int col);
static void     InsertRowInMatrix(matrix &big, longvector &row , const int r);
void            setZero();
 

/* this function assembles a small matrix A inside a big one B, with the idata list
   giving the map from the rows and columsn of A to the rows and columsn of B */
static void    AssembleMatrixInMatrix(matrix& block , matrix& big, idata list);

 
/* functions to get elements or blocks of the matrix */
matrix		block          (const int row1, const int row2, const int col1, const int col2);
longvector  MatrixColumn         (matrix m , int   c);
longvector  diagonal() const;
double		MatrixElement        (matrix m , int   r , int   c);
matrix		MatrixMinor          (matrix m,  int   r,  int c);
longvector  MatrixRow            (matrix m , int   r);


/* math operations on matrices */
static void     AddBtCB(matrix& m, const matrix& b, const matrix& c, double f);           /* m += B^t C B *fact  */
void            AddSparseBtCB(matrix& m, const matrix& b, const matrix& c,  double f);      /* m += B^t C B *fact  */
static void     AddBtS (longvector& r, const matrix& b, const longvector& s, double f);	/* r += B^t S *fact */
		
void    AddSparseBtS (longvector r , matrix b , longvector s , double f); /* r += B^t S *fact */
static  void    AddMatrices(matrix& m1, matrix& m2);                          /* stores result in m2 */
void    chsgn();                                        /* change sign */
void    InvertMatrix(matrix m);
double  determinant() const;
		
void    eigendata(matrix& elongvectors, longvector& evaluesR, longvector& evaluesI);
	double  norm(const char normtype='F');
		
void    MatrixTimesMatrix(double a, matrix A, matrix B, double b, matrix C); /* C= a*A*B+b*C */
void    MatrixTrTimesMatrix(double a, matrix A, matrix B, double b, matrix C); /* C=a*At*B+b*C */
void    MatrixTimesVector(double a, matrix A, longvector V, double b , longvector B); // B = b*B+a*A*V 
void    MatrixTrTimesVector(double a, matrix A, longvector V, double b, longvector B);// B = b*B+a*A'*V
		void    round() const;
void    SymmetrizeMatrix(matrix m);
void    VectorTrTimesMatrix(double a, longvector A, matrix B, double b, matrix C); /* C=a*At*B+b*C */


/* function to LU decompose the matrix, with no pivotting */
bool     factorLU();


/* solve Ax=b, overwriting b with x. No pivoting, for general unsymmetric matrices */
int     solveFull(longvector& b);
void    solveLU(longvector& b);



/* functions to display information of the matrix. EV: eigenvalues/longvectors */
void    print();
void    PrintMatrixEV(matrix m);
void    PrintMatrixToFile(matrix m, FILE *fp);

/* functions to transform the matrix */
void    SwapColumns(matrix m , int col1 , int col2);
void    SwapRows   (matrix m , int row1 , int row2);


};	

#endif

