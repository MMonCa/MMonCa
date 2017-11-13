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
/* matrix.c
 * Matrix manipulation.
 * Right now limited to double precission for real numbers.
 * Ignacio Romero (August 1999)
 *
 * Implementation for full matrices.
 */

#define NUM_PER_LINE   5         /* numbers per line, for screen output */

#include <ctype.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "Math/matrix.h"
#include "Math/vector.h"
#include "General/idata.h"
#include "Io/message.h"
#include "Main/feliks.h"

#ifdef MAGERIT
#define dgemv_ dgemv
#define dgemm_ dgemm
#endif


matrix :: matrix() :
	datatype(DATA_DOUBLE),
	rows(0),
	cols(0),
	start(0),
	data(0)
{
}




/* 'data' is a pointer to a list of pointers, each of these pointing
 * a row in the matrix. For compactness, all the data is allocated
 * at once, in a big block and then the pointers suitably directed.
 */
matrix :: matrix(const size_t rows_, const size_t cols_) :
datatype(DATA_DOUBLE),
rows((int) rows_),
cols((int) cols_)
{
	data  =  new double* [rows];
	start =  new double[rows*cols];
	
	// set the pointers in data to each row
	for (int i=0; i<rows; i++)
	{
		data[i] =  &(start[i*cols]) ;
		for (int j=0; j<cols; j++)  data[i][j] = 0.0;
	}
}




matrix :: matrix(const matrix& n) :
rows(n.rows),
cols(n.cols),
datatype(n.datatype)
{
	data  =  new double* [rows];
	start =  new double  [rows*cols];
	
	// set the pointers in data to each row
	for (int i=0; i<rows; i++)
	{
		data[i] =  &(start[i*cols]) ;
		for (int j=0; j<cols; j++)  data[i][j] = n.data[i][j];
	}
}




matrix :: ~matrix()
{
	delete [] start;
	delete [] data;
}




void  matrix :: addBlock(matrix& block, const int row , const int col)
{
	int      i, j , nrows, ncols;
	double **BIG, **BLOCK;
	
	BIG   = this->data;
	BLOCK = block.data;
	nrows = block.rows;
	ncols = block.cols;
	
	if ( (row + nrows > this->rows) || (col + ncols > this->cols))
		printf("Block too big in InsertBlockInMatrix");
	
	else
		for (i=0; i<nrows; i++)
			for(j=0; j<ncols; j++)
				BIG[row+i][col+j] += BLOCK[i][j];
}


extern "C" {void dgemv_(char*, const int*, const int*, double*, const double*, const int*,
                       const double*, int*, double*, double*, int*);}


/* r += B^T * S * fact */
void matrix :: AddBtS (longvector& r, const matrix& b, const longvector& s , double f)
{
	char notr='N';
	int  inc=1;
	double one=1.0;
	
	if (b.cols != r.length || b.rows != s.length) 
    {
		printf("\n Error in AddBtS: dimensions don't match. (r  += Bt S)");
		printf("\n r is %d, B is %d x %d and S is %d",
			   r.length, b.rows, b.cols, s.length);
    }
	else
		dgemv_(&notr   , &(b.cols) , &(b.rows) , &f  , b.start , &(b.cols) ,
			   s.data , &inc       , &one       , r.data , &inc);
}




/* r += B^T * S * f, using the sparsity structure of the B-matrix for standard
 finite elements. For 2D, B has dimensions 4x(2xnnode). For 3D, B is 6x(3xnnode)
 and is of the form 
 
 [ N,1  0 ]       [ N,1 0    0  ]
 [ 0   N,2]       [ 0   N,2  0  ]
 2D =  [ x   x  ]  3D = [ 0   0    N,3]
 [ N,2 N,1]       [ 0   N,3  N,2]
 [ N,3 0    N,1]
 [ N,2 N,1  0  ]
 */
void matrix :: AddSparseBtS (longvector r , matrix b , longvector s , double f)
{
	int     nnode, k, a;
	double **B, *R, *S;
	
	if (b.cols != r.length || b.rows != s.length) 
    {
		printf("\n Error in AddSparseBtS: dimensions don't match. (r  += Bt S)");
		printf("\n r is %d, B is %d x %d and S is %d",
			   r.length, b.rows, b.cols, s.length);
    }
	else if ( b.rows == 4)
    {
		nnode = b.cols/2;
		B = b.data;
		R = r.data;
		S = s.data;
		for (a = 0; a<nnode; a++)
		{
			k = 2*a;
			R[k]   += (B[0][k]  * S[0] + B[2][k]  *S[2] + B[3][k]  *S[3])*f;
			R[k+1] += (B[1][k+1]* S[1] + B[2][k+1]*S[2] + B[3][k+1]*S[3])*f;
		}
    }
	else if (b.rows == 6)
    {
		nnode = b.cols/3;
		B = b.data;
		R = r.data;
		S = s.data;
		for (a = 0; a<nnode; a++)
		{
			k = 3*a;
			R[k]   += (+B[0][k]  *S[0] + B[4][k]  *S[4]  + B[5][k]  *S[5] )*f;
			R[k+1] += (+B[1][k+1]*S[1] + B[3][k+1]*S[3]  + B[5][k+1]*S[5] )*f;
			R[k+2] += (+B[2][k+2]*S[2] + B[3][k+2]*S[3]  + B[4][k+2]*S[4] )*f;
		}
    }
	else
		printf("\n ERROR AddSparseBtS is not defined for this type of matrix");
}



extern "C" {void dgemm_(char*, char*, const int*, const int*, const int*, double*,
                       const double*, const int*, const double*, const int*, double*, 
                       double*, int*);}


/* M += B^T * C *  B fact
 This is one of the most important functions, because it is called very frequently. It
 is important to optimize it as much as possible */
void matrix :: AddBtCB(matrix& m , const matrix& b , const matrix& c, double fact)
{
    int tcols, trows;
    char no='N' , tr='T';
    double *t=NULL, tfix[270];
    double one=1.0 , zero=0.0;
	
    if (c.cols != b.rows)
    {
        printf("\n Error in AddBtCB: dimensions don't match. (M  += Bt C B)");
        printf("\n M is %d x %d, B is %d x %d and C is %d x %d",
               m.rows, m.cols, b.rows, b.cols, c.rows, c.cols);
    }
    else
    {
        /* Compute first C*B */
        tcols = b.cols;
        trows = c.rows;
        if (trows*tcols <= 270)
            t = tfix;
        else
            t = (double *)malloc( trows*tcols*sizeof(double));
		
        /*  this is tricky. We obtain t^T by multiplying B^T*C^T. Now, each of these
		 transposed matrices is passed as a C block, which means that for the Fortran
		 blas they are already transposed, and that's why we say there is no need to
		 transpose them again */
        /* t = CB */
        dgemm_(&no, &no   , &(b.cols) , &(c.rows) , &(b.rows) , &one ,
               b.start   , &(b.cols) ,   c.start , &(c.cols) , &zero  ,
               t          , &tcols);
		
        /* then m += fact* B^T (CB) */
        dgemm_(&no, &tr   , &tcols     , &(b.cols) , &trows     , &fact ,
               t         , &tcols     ,   b.start , &(b.cols) , &one  ,
               m.start   , &(m.cols));
		
        if (trows*tcols > 270) free(t);
    }
}



/* M += B^T * C *  B * fact 
 using the sparsity structure of the B-matrix for standard
 finite elements. For 2D, B has dimensions 3x(2xnnode). For 3D, B is 6x(3xnnode)
 and is of the form
 
 [ N,1  0 ]       [ N,1 0    0  ]
 [ 0   N,2]       [ 0   N,2  0  ]
 2D =  [ x   x  ]  3D = [ 0   0    N,3]
 [ N,2 N,1]       [ 0   N,3  N,2]
 [ N,3 0    N,1]
 [ N,2 N,1  0  ]
 
 No assumption on the structure of C is made. In particular, C may by full
 and unsymmetric
 */
void  matrix :: AddSparseBtCB(matrix& m, const matrix& b, const matrix& c, const double f)
{
	int a, k, i, nnode;
	matrix t;
	double **M, **B, **C, **T;
	
	if (c.cols != b.rows) 
    {
		printf("\n Error in AddBtCB: dimensions don't match. (M  += Bt C B)");
		printf("\n M is %d x %d, B is %d x %d and C is %d x %d",
			   m.rows, m.cols, b.rows, b.cols, c.rows, c.cols);
    }
	else
    {
		t.resize(c.rows , b.cols);
		T = t.data;
		M = m.data;
		B = b.data;
		C = c.data;
		
		if (b.rows == 4) 
		{
			nnode = b.cols/2;
			
			/* Compute first C*B */
			for (i=0; i<4; i++)
				for (a=0; a<nnode; a++)
				{
					k = 2*a;
					T[i][k] = C[i][0]*B[0][k] + C[i][2]*B[2][k] + C[i][3]*B[3][k];
					k ++;
					T[i][k] = C[i][1]*B[1][k] + C[i][2]*B[2][k] + C[i][3]*B[3][k];
				}
			
			/* then B^T CB */
			for (a=0; a<nnode; a++)
			{
				k = 2*a;
				for (i=0; i< b.cols; i++)
				{
					M[k][i]   += f*(+ B[0][k]  *T[0][i] + B[2][k]  *T[2][i] + B[3][k]*T[3][i]);
					M[k+1][i] += f*(+ B[1][k+1]*T[1][i] + B[2][k+1]*T[2][i] + B[3][k+1]*T[3][i]);
				}
			}
			
		}
		else if (b.rows == 6) 
		{
			nnode = b.cols/3;
			
			/* Compute first C*B */
			for (i=0; i<6; i++)
				for (a=0; a<nnode; a++)
				{
					k = 3*a;
					T[i][k]   = C[i][0]*B[0][k] + C[i][4]*B[4][k] + C[i][5]*B[5][k]; 
					k++;
					T[i][k] = C[i][1]*B[1][k]   + C[i][3]*B[3][k] + C[i][5]*B[5][k];
					k++;
					T[i][k] = C[i][2]*B[2][k]   + C[i][3]*B[3][k] + C[i][4]*B[4][k];
				}
			
			/* then B^T CB */
			for (a=0; a<nnode; a++)
			{
				k = 3*a;
				for (i=0; i< b.cols ; i++)
				{
					M[k][i]   += f*(+ B[0][k]  *T[0][i] 
									+ B[4][k]  *T[4][i]  
									+ B[5][k]  *T[5][i]);
					M[k+1][i] += f*(+ B[1][k+1]*T[1][i] 
									+ B[3][k+1]*T[3][i]
									+ B[5][k+1]*T[5][i]);
					M[k+2][i] += f*(+ B[2][k+2]*T[2][i] 
									+ B[3][k+2]*T[3][i]
									+ B[4][k+2]*T[4][i]);
				}
			}
		}
		
    }
}








void matrix :: AddMatrices(matrix& m1, matrix& m2)
{
    int r,c;
    double **M1, **M2;
	
    if ( (m1.cols != m2.cols) || (m1.rows != m2.rows) )
		printf("Error in AddMatrices. Matrices of different size");
	
    else
	{
		M1 = m1.data;
		M2 = m2.data;
		for (r=0; r < m1.rows; r++)
			for (c=0; c < m1.cols; c++)
				M2[r][c] += M1[r][c];
	}
}






/* this function "assembles" a small matrix into a large one. The elements of the
 small matrix are ADDED to the positions indicated by the idata list as
 customary in finite element assembly procedures. No bounds are checked for the
 positions of the assembled block in order to increase speed 
 */
void  matrix ::  AssembleMatrixInMatrix(matrix& block , matrix& big, idata list)
{
	int r,c;
	int globalcol, globalrow;
	int len;
	double **Big, **Block;
	
	Big   = big.data;
	Block = block.data;
	len   = IdataSize(list);
	
	for (r=0; r<len; r++)
    {
		globalrow = GetIntInIdata(r, list);
		for (c=0 ; c<len ; c++)
		{
			globalcol = GetIntInIdata(c, list);
			Big[globalrow][globalcol] += Block[r][c];
		}
    }
}



/* note: the integers row1, row2, col1, col2 are start form 0
 */
matrix matrix :: block(const int row1, const int row2, const int col1, const int col2)
{
	int r,c;
	matrix bl;
	
	if ( (row1 > row2) || (col1 > col2) ||  (row2 > rows-1) || (col2 > cols-1)) 
		printf("Error in MatrixBlock");
	
	else
    {
		bl.resize(row2-row1+1, col2-col1+1);
		for (r=0; r<=row2-row1; r++)
			for (c=0; c<=col2-col1; c++) 
				bl.data[r][c] = data[row1+r][col1+c];
    }
	return bl;
}



void matrix :: chsgn()
{
	int r , tot;
	matrix& m=*this;
	tot = m.rows*m.cols;
    double *M   = m.start;
    for (r=0; r<tot; r++) M[r] *= -1.0;
}



matrix matrix ::  IdentityMatrix(int n)
{
	int i;
	matrix m(n,n);
	
	for (i=0; i<n; i++) m.data[i][i] = 1.0;
	
	return m;
}



/*---------------------------------------------------------------------*/
/*                             InsertBlockInMatrix                     */
/*---------------------------------------------------------------------*/
void  matrix ::  InsertBlockInMatrix(matrix& block, matrix& big, const int row, const int col)
{
    int i,j;
	
    if ( (row + block.rows > big.rows) || (col + block.cols > big.cols))
        printf("Block too big in InsertBlockInMatrix");
	
    else
	{
        for (i=0; i<block.rows; i++)
            for(j=0; j<block.cols; j++)
                big.data[row+i][col+j] = block.data[i][j];
		
	}
}



void  matrix :: InsertElement(matrix m , double el , int row, int col)
{
	if (row > m.rows) 
		printf("in InsertElement. Row too large");
	else if (col > m.cols) 
		printf("in InsertElement. Col too large");
	else
		m.data[row][col] = el;
}



void  matrix ::  InsertRowInMatrix(matrix& m, longvector& row, int rowNumber)
{
	int k;
	
	if( row.length > m.cols)
		printf("ERROR in InsertRowInMatrix. Row too large");
	
	else
		for (k=0 ; k<m.cols; k++) m.data[rowNumber][k] = row.data[k];
}



/* we invert the matrix m by solving m * X = identity, with the lapack function dgesv.
 the inverse overwrites the output  */
void  matrix :: InvertMatrix(matrix m)
{
	/*
	 int    i, j;
	 matrix id;
	 int    n , *p , info;
	 double idet, **M, tmp;
	 double old[3][3];
	 extern void dgesv_();
	 
	 
	 n = m.rows;
	 M = m.data;
	 
	 switch(n)
	 {
	 case 1:
	 M[0][0] = 1.0/M[0][0];
	 break;
	 
	 case 2:
	 idet     = 1.0/(M[0][0]*M[1][1] - M[0][1]*M[1][0]);
	 tmp      = M[0][0];
	 M[0][0]  = idet * M[1][1];
	 M[1][1]  = idet * tmp;
	 M[0][1] *= -idet;
	 M[1][0] *= -idet;
	 break;
	 
	 case 3:
	 idet = M[0][0]*M[1][1]*M[2][2] + M[1][0]*M[2][1]*M[0][2] 
	 +    M[2][0]*M[0][1]*M[1][2] - M[0][2]*M[1][1]*M[2][0] 
	 -    M[0][1]*M[1][0]*M[2][2] - M[0][0]*M[2][1]*M[1][2];
	 
	 if (idet == 0.0) 
	 Message("\n Error in invert matrix. Singular matrix");
	 
	 else
	 {
	 idet = 1.0/idet;
	 
	 for (i=0;i<3;i++) for(j=0;j<3;j++) old[i][j] = M[i][j];
	 
	 M[0][0] = (old[1][1]*old[2][2]-old[1][2]*old[2][1])*idet;
	 M[0][1] = (old[0][2]*old[2][1]-old[0][1]*old[2][2])*idet;
	 M[0][2] = (old[0][1]*old[1][2]-old[0][2]*old[1][1])*idet;
	 
	 M[1][0] = (old[2][0]*old[1][2]-old[1][0]*old[2][2])*idet;
	 M[1][1] = (old[0][0]*old[2][2]-old[0][2]*old[2][0])*idet;
	 M[1][2] = (old[0][2]*old[1][0]-old[0][0]*old[1][2])*idet;
	 
	 M[2][0] = (old[1][0]*old[2][1]-old[1][1]*old[2][0])*idet;
	 M[2][1] = (old[0][1]*old[2][0]-old[0][0]*old[2][1])*idet;
	 M[2][2] = (old[0][0]*old[1][1]-old[0][1]*old[1][0])*idet;
	 }
	 break;
	 
	 default:
	 id = IdentityMatrix(n);
	 
	 p  = (int *) Allocate (n *sizeof(int));
	 dgesv_(&n, &n, m.start, &n, p, id.start , &n , &info);
	 
	 if (info > 0) printf("\n Error in InvertMatrix: singular matrix.");
	 
	 Deallocate(p);
	 Deallocate(m.start);
	 Deallocate(m.data);
	 m.start = id.start;
	 m.data  = id.data;
	 }
	 */
}



/* LU decomposition for full matrices with no pivoting.
 * Replaces matrix m with LU,
 * unit Lower triangular and U, upper triangular. Of course, m is lost.
 */
bool matrix :: factorLU()
{
	int r,c;
	int i,j,k;
	double aii;
	double **M;
	bool   ret(true);
	matrix& m=*this;
	
	c = cols;
	r = rows;
	
	if (c != r ) 
	{
		printf("ERROR in FullLU_NoPiv. Only square matrices!");
		ret = false;
	}
	
	else
    {
		M = m.data;
		for (i=0; i<= r-2 ;i++)
		{
			aii = M[i][i];
			if (aii == 0.0) 
			{
				printf("FullLU_NoPiv. Zero pivot number %d.", i+1);
				ret = false;
				break;
			}
			aii = 1.0/aii;
			for (j=i+1; j<r; j++) M[j][i] *= aii;
			for (j=i+1; j<r; j++)
				for (k=i+1; k<r; k++) M[j][k] -= M[j][i]*M[i][k];
			ret = true;
		}
    }
	return ret;
}




double matrix :: determinant() const
{
	/*
	 double det=0.0, sign;
	 int i;
	 
	 if ( rows == 3 || cols == 3)
	 {
	 det = start[0*3 + 0]*start[1*3 + 1]*start[2*3 + 2]
	 +     start[1*3 + 0]*start[2*3 + 1]*start[0*3 + 2]
	 +     start[2*3 + 0]*start[0*3 + 1]*start[1*3 + 2]
	 -     start[0*3 + 2]*start[1*3 + 1]*start[2*3 + 0]
	 -     start[0*3 + 1]*start[1*3 + 0]*start[2*3 + 2]
	 -     start[0*3 + 0]*start[2*3 + 1]*start[1*3 + 2];
	 }
	 else
	 {
	 sign = +1.0;
	 for (i=0; i<rows; i++)
	 {
	 matrix minor = MatrixMinor(*this, 0, i);
	 if (data[0][i] != 0.0) det  += data[0][i]*sign*minor.determinant();
	 sign *= -1.0;
	 }
	 }
	 */
	double det = 0.0;
	return det;
}








/*---------------------------------------------------------------------*/
/*                             MatrixColumn                            */
/*---------------------------------------------------------------------*/
longvector matrix ::  MatrixColumn(matrix m, int c)
{
	longvector v;
	int i;
	
	if (c>m.cols-1)
		printf("MatrixColumn. This matrix does not have this column");

	else
    {
		v.resize(m.rows);
		for (i=0; i<m.rows; i++) v.data[i] = m.data[i][c];
    }
	
	return v;   
}



longvector matrix :: diagonal() const
{
	longvector v(cols);
	for (int i=0; i<cols; i++) v.data[i] = data[i][i];
	return v;
}


/*---------------------------------------------------------------------*/
/*                             MatrixElement                           */
/*---------------------------------------------------------------------*/
double  matrix :: MatrixElement(matrix m,int r, int c)
{
	double ret;
	
	if ( r>=m.rows || c>=m.cols || r<0 || c<0)
    {
		printf("Error in MatrixElement. Matrix has dimensions %d,%d",m.rows,m.cols);
		printf("and the element requested is %d,%d",r,c);
		ret = 0.0;
    }
	else      
		ret =  m.data[r][c];
	
	return ret;
}

extern "C" {void dgeev_(char*, char*, int*, double*, int*, double*, double*,
						double*, int*, double*, int*,
						double*, int*,
						int*);}



/* computes the eigenvalues and eigenvectors of a general full, unsymmetric matrix
 with a lapack routine */
void matrix :: eigendata(matrix& evectors, longvector& evaluesR, longvector& evaluesI)
{
	int   r, lwork, info , dummyint=1;
	char  N='N', V='V';
	double *work=NULL, *dummy=NULL;
	 
	 
	r = rows;
	evaluesR.resize(r);
	evaluesI.resize(r);
	
	lwork = 5*r;
	work = (double *) malloc( lwork * sizeof(double));
	dgeev_(&N, &V          , &r ,
		   this->start     , &r , 
		   evaluesR.data  , evaluesI.data ,
		   dummy           , &dummyint ,
		   evectors.start , &r, 
		   work            , &lwork, 
		   &info);
	free(work);
}



matrix  matrix :: MatrixMinor(matrix m, int row, int col)
{
    matrix minor;
	/*
	 int i, j, r, c;
	 
	 
	 if (m.rows == 1 || m.cols == 1) return NULL;
	 
	 minor = NewMatrix(m.rows-1, m.cols-1);
	 
	 r = -1;
	 for (i=0; i<m.rows; i++)
	 {
	 if (i != row)
	 {
	 r++;
	 c = -1;
	 for (j=0; j<m.cols; j++)
	 {
	 if (j != col)
	 {
	 c++;
	 minor.data[r][c] = m.data[i][j];
	 }
	 }
	 }
	 }
	 */
	return minor;
}



double  matrix :: matrix :: norm(const char type)
{
	int i,j;
	matrix& m=*this;
	double n=0, **M=NULL, sum;
	
	M = m.data;
	
	if (type == 'F' || type == 'f')
	{
		// Frobenius norm
		for (i=0; i < m.rows; i++)
			for (j=0; j < m.cols; j++)
				n += M[i][j] * M[i][j];
		n = sqrt(n);
	}
	
	// 1-norm, maximum column sum
	else if (type == '1')
	{
		n = 0.0;
		for (i=0; i < m.cols; i++)
		{
			sum = 0.0;
			for (j=0; j < m.rows; j++) sum += fabs(M[j][i]);
			n = std::max(n, sum);
		}
	}
	
	
	/* infinite-norm, maximum row sum */
	else if (type == 'i' || type == 'I')
	{
		n = 0.0;
		for (i=0; i < m.rows; i++)
		{
			sum = 0.0;
			for (j=0; j < m.cols; j++) sum += fabs(M[i][j]);
			n = std::max(n, sum);
		}
	}
	
	return n;
}






longvector  matrix :: MatrixRow(matrix m, int r)
{
	longvector v;
	int i;
	
	if (r>m.rows-1)
		printf("ERROR in MatrixRow. This matrix does not have this row");
	else
	{
		v.resize(m.cols);
		for (i=0; i<m.cols; i++) v.data[i] = m.data[r][i];
	}
	return v;
}



/* using BLAS m3 = alpha * m1 * m2 + beta * m3     */
/* since BLAS works with fortran conventions, the matrices must be transposed before
 using them, and the result transposed before exiting, 
 even simpler, just swap the order */
void  matrix :: MatrixTimesMatrix(double alpha, matrix m1, matrix m2, double beta , matrix m3)
{
	/*
	 char notr='N';
	 extern void dgemm_();
	 
	 
	 if ( m1.cols != m2.rows) 	  
	 printf("Bad dimensions in MatrixTimesMatrix");
	 
	 else
	 {
	 if (m3 == NULL) m3 = NewMatrix(m1.rows, m2.cols);
	 dgemm_(&notr, &notr, &(m2.cols) , &(m1.rows) , &(m2.rows) , &alpha , 
	 m2.start   , &(m2.cols) , m1.start   , &(m1.cols) , &beta  , 
	 m3.start   , &(m3.cols));
	 }
	 */
    /*  The good old way
	 for (i=0; i<m1.rows; i++)
	 for (j=0; j<m2.cols; j++)
	 for (k=0; k<m1.cols; k++) m3.data[i][j] += m1.data[i][k] * m2.data[k][j]; */
}



/* using BLAS m3 = alpha * m1 * m2 + beta * m3     */
/* since BLAS works with fortran conventions, the matrices must be transposed before
 using them, and the result transposed before exiting, 
 even simpler, just swap the order */
void  matrix :: MatrixTrTimesMatrix(double alpha, matrix m1, matrix m2, double beta , matrix m3)
{
	/*
	 char tr='T' , notr='N';
	 extern void dgemm_();
	 
	 
	 if ( m1.rows != m2.rows) 	 
	 printf("Bad dimensions in MatrixTimesMatrix");
	 
	 else
	 {
	 if (m3 == NULL) m3 = NewMatrix(m1.cols, m2.cols);
	 dgemm_(&notr, &tr  , &(m2.cols) , &(m1.cols) , &(m2.rows) , &alpha , 
	 m2.start   , &(m2.cols) , m1.start   , &(m1.cols) , &beta  , 
	 m3.start   , &(m3.cols));
	 }
	 */
}




/*---------------------------------------------------------------------*/
/*                             MatrixTimesVector                       */
/*---------------------------------------------------------------------*/
/*  r = alpha*M * v + beta * r */
void  matrix :: MatrixTimesVector(double alpha , matrix m, longvector v, double beta , longvector r)
{
	/*
	 char tr='T';
	 int  inc=1;
	 extern void dgemv_(), zgemv_();
	 //int i, j;
	 //complex double *R, *V, **M;
	 
	 
	 if (m.cols != v.length)  
	 {
	 printf("Error in MatrixTimesVector. Bad sizes");
	 return;
	 }
	 
	 else if (r == NULL) 
	 printf("\n ERROR in MatrixTimesVector. NULL longvector.");
	 
	 else if (m.rows != r.length) 
	 printf("\n Error in MatrixTimes Vector. Matrix rows different from longvector length");
	 
	 else if (m.datatype == DATA_DOUBLE && v.datatype == DATA_DOUBLE)
	 dgemv_(&tr, &(m.cols), &(m.rows), &alpha, m.start, &(m.cols), v.data, &inc, &beta, r.data, &inc);
	 
	 else if (m.datatype == DATA_DCOMPLEX && v.datatype == DATA_DCOMPLEX)
	 {
	 //	zgemv_(&tr     , &(m.cols) , &(m.rows) , &alfa   , (complex double *)m.start , &(m.cols) ,
	 //		   (complex double *) v.data , &inc       , &beta      , (complex double *)r.data , &inc);
	 R = (complex double *) r.data;
	 M = (complex double **) m.data;
	 V = (complex double *) v.data;
	 for (i=0; i<v.length; i++)
	 for (j=0; j<v.length; j++)
	 R[i] = beta*R[i] + alpha * M[i][j]*V[j];
	 }
	 
	 
	 else
	 printf("\n Error in MatrixTimesVector: matrix and longvector are of different type");
	 */
}




/*  r = alpha*M' * v + beta * r */
void  matrix :: MatrixTrTimesVector(double alfa , matrix m, longvector v, double beta , longvector r)
{
	/*
	 char notr='N';
	 int  inc=1;
	 extern void dgemv_();
	 
	 if (m.rows != v.length)  
	 {
	 printf("Error in MatrixTrTimesVector. Bad sizes");
	 return;
	 }
	 
	 else if (r == NULL)
	 printf("\n ERROR in MatrixTrTimesVector. NULL longvector.");
	 
	 else if (m.cols != r.length) 
	 printf("\n Error in MatrixTimes Vector. Matrix cols different from longvector length");
	 
	 else
	 dgemv_(&notr   , &(m.cols) , &(m.rows) , &alfa   , m.start , &(m.cols) ,
	 v.data , &inc       , &beta      , r.data , &inc);
	 */
}


matrix& matrix :: operator=(const matrix& m)
{
	this->resize(m.rows, m.cols);
	for (int i=0; i< rows; i++)
		for (int j=0; j<cols; j++)
			data[i][j] = m.data[i][j];
	return *this;
}


matrix& matrix :: operator+=(const matrix &m)
{
	for (int i=0; i< rows; i++)
		for (int j=0; j<cols; j++)
			data[i][j] += m(i,j);
	
	return *this;
}


matrix&	matrix ::	operator-=(const matrix &m)
{
	for (int i=0; i< rows; i++)
		for (int j=0; j<cols; j++)
			data[i][j] -= m(i,j);
	
	return *this;
}



void matrix :: print()
{
	int i,j,pstart,pend, c;
	matrix& m=*this;
	c = m.cols;
	
	if (m.rows > 0 && c > 0)
    {
		printf("\n");
		
		// row loop
		for (i=0; i< m.rows; i++)
		{
			printf("\nRow %3d: ",i+1);
			pstart = 0;
			pend   = std::min<int>(c, NUM_PER_LINE);
			
			// column loop
			if (m.datatype == DATA_DOUBLE) 
				for (j=pstart; j<pend; j++) printf("% 10.5e ", m.data[i][j] );
			else
				for (j=pstart; j<pend; j++) printf("(% 6.3e,% 6.3e) ", m.data[i][2*j], m.data[i][2*j+1]);
			
			pstart += NUM_PER_LINE;
			pend    = std::min<int>( c , pend + NUM_PER_LINE);
			
			while (pstart <= c )
			{
				printf("\n\t ");
				if (m.datatype == DATA_DOUBLE) 
					for (j=pstart; j<pend; j++) printf("% 10.5e ", m.data[i][j] );
				else
					for (j=pstart; j<pend; j++) printf("(% 6.3e,% 6.3e) ", m.data[i][2*j], m.data[i][2*j+1]);
				
				pstart += NUM_PER_LINE;
				pend    = std::min<int>(m.cols , pend + NUM_PER_LINE);
			}
		}
    }
	fflush(stdout);
}





void  matrix :: PrintMatrixEV(matrix m)
{
	int n, a;
	
	n = m.rows;
	matrix vec(n,n);
	longvector valR(n);
	longvector   valI(n);
	
	m.eigendata(vec, valR, valI);
	
	printf("\n\n Matrix eigenvalues:");
	for (a=0; a<n; a++)
    {
		printf("\n     %2d : % 10.4e", a+1, valR.data[a]); 
		if ( sqrt(valI.data[a]*valI.data[a])  > 1e-6 )
			printf(" + %d* % 10.4e", a, valI.data[a]); 
    }
	
}



void   matrix ::   PrintMatrixToFile(matrix m, FILE *fp)
{
    int r, c;
	
    fprintf(fp, "%d\n%d", m.rows, m.cols);
    for (r=0; r<m.rows; r++)
        for (c=0; c<m.cols; c++)
            fprintf(fp, "\n%e", m.data[r][c]);
}



void matrix ::  resize( const size_t newrows, const size_t newcols)
{
	matrix &m=*this;
    
    if (m.rows != newrows || m.cols != newcols)
    {
        if (m.start != 0) free(m.start);
        if (m.data  != 0) free(m.data);
		
        m.rows  = (int) newrows;
        m.cols  = (int) newcols;
        
		m.data  = (double **) malloc(newrows         * sizeof(double *));
        m.start = (double  *) malloc(newrows*newcols * sizeof(double  ));
		
        // set the pointers in data to each row
        for (int i=0; i < newrows; i++)
            m.data[i] =  &(m.start[i*newcols]) ;
    }
	
	// in any case, set to Zero
	setZero();
}






void matrix :: round() const
{
	int i, tot;
	const matrix& m=*this;
	
	tot = m.rows * m.cols;
	for (i=0; i<tot; i++) m.start[i] *= ( fabs(m.start[i]) < 1e-30 ) ? 0.0 : 1.0;
}




int matrix :: solveFull(longvector& b)
{
    factorLU();   /* LU decomposition of full matrix */
	solveLU(b);
	
    return 1;
}


void matrix :: solveLU(longvector& b)
{
    int i,k,c;
	matrix &LU=*this;
	
    /* forward substitution */
    c = LU.cols;
    for (k=1; k<c; k++)
		for (i=0; i<k; i++) 
			b.data[k] -= LU.data[k][i]*b.data[i];
	
    /* back substitution */
    b.data[c-1] /= LU.data[c-1][c-1];
    for (k=c-2; k>=0; k--)
	{
		for (i=k+1; i<c; i++)
			b.data[k] -= LU.data[k][i]*b.data[i];
		b.data[k] /= LU.data[k][k];
	}
}




void    matrix ::  SymmetrizeMatrix(matrix m)
{
    int r, c, i, j;
	
    r = m.rows;
    c = m.cols;
	
    for (i=0; i<r; i++)
        for (j=0; j<i; j++)
        {
            m.data[i][j] = 0.5*(m.data[i][j]+m.data[j][i]);
            m.data[j][i] = m.data[i][j];
        }
}



void  matrix ::  SwapColumns(matrix m, int col1, int col2)
{
	int r;
	double temp;
    
	if ( col1 > m.cols-1  || col2 > m.cols-1 ) 
		printf("\n ERROR in SwapColumns");
	
	else if ( col1 == col2 ) 
		return;
	
	else 
    {
		for (r = 0; r<m.rows; r++)
		{
			temp             = m.data[r][col1];
			m.data[r][col1] = m.data[r][col2];
			m.data[r][col2] = temp;
		}
    }
	return;
}



/*---------------------------------------------------------------------*/
/*                             SwapRows                                */
/*---------------------------------------------------------------------*/
void  matrix ::  SwapRows(matrix m, int row1, int row2)
{
	int c;
	double temp;
	
	if ( row1 == row2 ) 
		return;
	
	else if ( row1 > m.rows-1  || row2 > m.rows-1 ) 
		printf("\n Error in SwapRows");
    
	else
    {
		for (c = 0; c<m.cols; c++)
		{
			temp             = m.data[row1][c];
			m.data[row1][c] = m.data[row2][c];
			m.data[row2][c] = temp;
		}
    }
	return;
}





/* using BLAS m3 = alpha * m1 * m2 + beta * m3     */
/* since BLAS works with fortran conventions, the matrices must be transposed before
 using them, and the result transposed before exiting, 
 even simpler, just swap the order */
void  matrix :: VectorTrTimesMatrix(double alpha, longvector  v, matrix m, double beta , matrix m3)
{
	/*
	 char notr='N'; 
	 int  one=1.0;
	 extern void dgemm_();
	 
	 
	 if ( v.length != m.rows) 	  
	 printf("Bad dimensions in VectorTimesMatrix");
	 
	 else
	 {
	 if (m3 == NULL) m3 = NewMatrix(1, m.cols);
	 dgemm_(&notr, &notr, &(m.cols) , &one       , &(m.rows) , &alpha , 
	 m.start    , &(m.cols) , v.data    , &(v.length) , &beta  , 
	 m3.start   , &(m3.cols));
	 }
	 */
    /*  The good old way
	 for (i=0; i<m1.rows; i++)
	 for (j=0; j<m2.cols; j++)
	 for (k=0; k<m1.cols; k++) m3.data[i][j] += m1.data[i][k] * m2.data[k][j]; */
}







void matrix :: setZero()
{
	Dzero( start, rows*cols);
}
