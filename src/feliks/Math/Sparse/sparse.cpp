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
/* sparse.cpp
 *
 * functions to work with sparse matrices, following Harwell-Boeing convention.
 *
 * i. romero, october 2001
 *
 * right now only works with Harwell-Boeing storage, but in the future other type
 * of storages should be included
 *
 */
#include <iostream>
#include <stdlib.h>
#include "Math/Sparse/sparse.h"
#include "Math/vector.h"

#include "Model/model.h"
#include "Math/matrix.h"
#include "Io/message.h"

// leave memory at the end, to avoid conficts with dmalloc



sparseMatrix :: sparseMatrix() :
	datatype(DATA_DOUBLE) ,
	condition(0.0),
	factorized(false),
	ddata(NULL),
	rows(0),
	cols(0),
	numbering(NUM_C),
	dumpfmt("")
{
}

 sparseMatrix :: sparseMatrix(model &m, char datatype) :
	datatype(DATA_DOUBLE) ,
	condition(0.0),
	factorized(false),
	ddata(NULL),
	rows(0),
	cols(0),
	numbering(NUM_C),
	dumpfmt("")
{
	// common initializations for all sparseMatrix types
	condition = 0.0;
}


sparseMatrix :: ~sparseMatrix()
{
	delete [] ddata;
}




double sparseMatrix :: getEntry(int row, int col) const
{
	WarningMessage("\n getEntry. SparseMatrix type not implemented");
	return 0.0;
}




void sparseMatrix :: info(ostream &of)
{
	static const char *storageTnames[] = {"Dense", "CSC", "CSC", "SSK", "SSS"};
	static const char *numberingTnames[] = {"Fortran", "C"};
	static const char *dataTnames[] = {"Single precision", "Double precision", "Single complex", "Double complex"};

	of  << endl << endl
		<< "\n System matrix:"
		<< "\n\tNumber of rows             : " << rows
		<< "\n\tNumber of columns          : " << cols
		<< "\n\tCondition number           : " << condition
		<< "\n\tFactorized                 : " << factorized
		<< "\n\tNumbering style            : " << numberingTnames[numbering]
		<< "\n\tStorage                    : " << storageTnames[storage]
		<< "\n\tData type                  : " << dataTnames[datatype];

	if (dumpfmt == "matrixmarket" || dumpfmt == "mm")
	{
		of << "\n\tTangent entries in file    : feliks.tangent.mm";
		dump(dumpfmt);
	}

	else if (dumpfmt == "harwellboeing" || dumpfmt == "hb" || dumpfmt == "harwell-boeing")
	{
		of << "\n\tTangent entries in file    : feliks.tangent.hwb";
		dump(dumpfmt);
	}
}




void  sparseMatrix :: print(FILE *fp)
{
	std::cout<<"\n print sparseMatrix. SparseMatrix type not implemented";
}



// creates a full matrix by filling up entries from sparse matrix
matrix sparseMatrix :: toFull()
{
	matrix full(rows,cols);
	if (datatype == DATA_DOUBLE)
	{
		for (int k=0; k< rows; k++)
			for (int j=0; j<cols; j++)
				full.data[k][j] = getEntry(k,j);
	}
	
	else
	{
	}		
	
	return full;
}







/* Ejemplo de un problema estandar de autovalores y autovectores con una matriz simétrica
la matriz de nombre bcsstk01.rsa en el matrix market */ 

#ifdef testSymmetric

#include "Math/Sparse/compressed.h"
#include "Math/matrix.h"
#include <iostream>


int main()
{
	
const	int nval  = 4;		//num de autovalores deseado


//Definimos la matriz sparse

	SCSCMatrix A("/home/borja/Escritorio/programacion/verano_2008/matrix_market/matrices/bcsstk01.rsa");
	matrix Amatrix=NULL;
	const	int sizeA  = A.getNColumns();

	matrix Amat = NULL;
	Amat = A.toFull();
	

//Definimos el problema estándar de los autovalores
	const	int sizevect = sizeA*nval;
	double eval[nval];
	double evali[nval];
	matrix evect = NewMatrix(sizeA,nval);

	A.computeEigenSystem(nval, evect, eval);
	cout << "The eigenvalues calculate with the function computeEigenSystem are:" << endl;
	for(int i=0; i<nval;i++) cout << i <<" : "<<eval[i] <<endl;	

	A.computeEigenValues(nval, eval);
	cout << endl<<endl<<"The eigenvalues calculate with the function computeEigenValues are:" << endl;
	for(int i=0; i<nval;i++) cout << i <<" : "<<eval[i] <<endl;

	return 0;
}
#endif





/* Ejemplo de un problema estandar de autovalores y autovectores con una matriz simétrica
la matriz de nombre bfw398a.rua en el matrix market */ 

#ifdef testUnSymmetric

#include "Math/Sparse/compressed.h"
#include "Math/matrix.h"
#include <iostream>


int main()
{
	
const	int nval  = 4;		//num de autovalores deseado


//Definimos la matriz sparse

	CSCMatrix A("/home/borja/Escritorio/programacion/verano_2008/matrix_market/matrices/bfw398a.rua");
	const	int sizeA  = A.getNColumns();

	matrix Amat = NULL;
	Amat = A.toFull();
	

//Definimos el problema estándar de los autovalores
	const	int sizevect = sizeA*nval;
	double eval[nval];
	double evali[nval];
	matrix evect = NewMatrix(sizeA,nval);

	A.computeEigenSystem(nval, evect, eval);
	cout << "The eigenvalues calculate with the function computeEigenSystem are:" << endl;
	for(int i=0; i<nval;i++) cout << i <<" : "<<eval[i] <<endl;	

	A.computeEigenValues(nval, eval);
	cout << endl<<endl<<"The eigenvalues calculate with the function computeEigenValues are:" << endl;
	for(int i=0; i<nval;i++) cout << i <<" : "<<eval[i] <<endl;

	PrintMatrixEV( Amat);

	return 0;
}
#endif

