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
 * csc.h
 * 
 * december 2002, ignacio romero
 *
 * data type and functions for sparse matrices in Compressed Sparse Column format , 
 * aka Harwell/Boeing format.
 */

#ifndef _compressed_matrix_h
#define _compressed_matrix_h

#include <fstream>
#include <cstring>

#include "Math/Sparse/sparse.h"
#include "Math/vector.h"
#include "Math/matrix.h"

#ifdef WITHPETSC
#include "petscda.h"
#endif



// the compressedMatrix is a virtual class to support
// all the data for the CSC, CSR, SSCS, SCSR formats
class compressedMatrix : public sparseMatrix{
	
protected:	
	int    *prow;		// vector of nnz ints, the row index of the elements Aij in ddata.
	int    *colstart;	// vector of N+1 ints, the start of each column.
	int    *pdiag;		// an array of ndof ints, the positions of the diagonal terms in ddata
	virtual void timesVectorDouble(longvector& b, longvector& Axb)=0;
	void	csc_matrix_sort_rowidx();
	void	copy_pairs2col (const void* col, int max_nnz, int j);
	void	copy_col2pairs (int j, void* col, int max_nnz);



	int		nnz;
	
	
	
public:
	compressedMatrix();
	compressedMatrix (const char *filename, bool symmetric=true);
	compressedMatrix(model &m, char datatype);
	virtual ~compressedMatrix();
	virtual void createProfile(const model &m)=0;


	// Functions to read a HB format file
	static void         convertDouble(char* num);
	void                ReadFormat( ifstream& file, int& n, int& fmt);
	bool                ReadEntry( ifstream& file, int nval, int fval, int& j, double& val);
	void                expand_symmetric_storage();

	void                setZero();
	virtual void        dump(string format="HB");
	virtual void        info(ostream &of=cout)=0;
	
	virtual void        convertToFortran();
	virtual void        convertToC();
	void                getDiagonal(longvector& diag);
	void                timesVector(longvector& b, longvector& c);
	virtual void        timesDoubleArray(double* v, double* Av)=0;
	virtual void        transpose();

	inline	int			getNProfileTerms() const{ return nnz;}
	inline  int*		getPdiag()const { return pdiag;}
	inline  int*		getProw(){ return prow;}
	inline  const int*  getProw() const { return prow;}
	inline  int*		getColstart(){ return colstart;}
	inline  const int*  getColstart() const{ return colstart;}
	int					getNTermsInColumn(const int col);
};





// the symmetric CSC matrix holds only the data for the diagonal and the lower??
// part of the matrix
class SCSCMatrix : public compressedMatrix{
	
private:

	void assembleDouble (const matrix &block, const std::vector<int>& position);
	void timesVectorDouble(longvector& b, longvector& Axb);
	
			
public:

	SCSCMatrix(model &m, const char datatype='d');
	SCSCMatrix(const SCSCMatrix &m);
	SCSCMatrix(const char *fname);	//constructor to build a scscmatrix according to a HB format file
	SCSCMatrix();
	~SCSCMatrix();
	SCSCMatrix& operator=(const SCSCMatrix &m);
	void createProfile(const model &m);

	
	// functions to add a double to the matrix or a block
	void  addTo(double x, int row, int col);
	void  assemble(const matrix& block, const std::vector<int>& list);
	void  insert(double x, int row, int col);
	
	// other functions
	void	check();
	void	dump(string format);
	double	getEntry(int row, int col);
	void	info(ostream &of=cout);
	void	print(FILE *fp=stdout);
	void	timesDoubleArray(double* v, double* Av);

	
	matrix	toFull();
};






class CSCMatrix : public compressedMatrix{
	
private:
	void	assembleDouble (const matrix& block, const std::vector<int>&  position);
	void	assembleComplex(const matrix& block, const std::vector<int>& position);		
	void	timesVectorDouble(longvector& b, longvector& Axb);
	friend class direct;
	
public:
	CSCMatrix();
	CSCMatrix(int size, int nnz, double *nzelements, int* rowvet, int *colvect);
	CSCMatrix(const char *filename);	//constructor to build a cscmatrix from HB file
	CSCMatrix(model &m, const char datatype='d');
	CSCMatrix(const CSCMatrix &m);
	CSCMatrix(const SCSCMatrix &m);

	~CSCMatrix();
	void createProfile(const model &m);


	CSCMatrix& operator=(const CSCMatrix &m);

	void	timesDoubleArray(double* v, double* Av);


	// functions to add a double to the matrix or a block
	void	addTo(double x, int row, int col);
	void	assemble(const matrix& block, const std::vector<int>& list);			
	void	insert(double x, int row, int col);


	// other functions
	void	check();
	double	getEntry(int row, int col);
	void	info(ostream &of=cout);
	void	print(FILE *fp=stdout);
	matrix	toFull();
	
	
	//Functions to solve the eigenvalue problem
	/*
	void	computeEigenSystem(int nval, matrix& evect, double eval[], double evali[]=NULL, const char *which="LM");
	void	computeEigenValues(int nval, double eval[], double evali[]=NULL, const char *which="LM");

	void	computeGeneraliseEigenSystem(compressedMatrix& M, int nval, matrix& evect, double eval[], double evali[]=NULL, char *which="LM");  
	void	computeGeneraliseEigenValues(compressedMatrix& M,int nval, matrix& evect, double eval[],  double evali[]=NULL, char *which="LM");
	 */
};




#endif
