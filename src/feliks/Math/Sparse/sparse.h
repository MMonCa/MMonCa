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
 * sparse.h
 * 
 * i. romero, october 2001
 *
 * and functions for sparse matrices. 
 *
 * Every sparse matrix is stored into a structure that includes, among other data, 
 * an array of integer arrays and an array of doubles.
 * This is flexible enough to accomodate different storage stategies.
 * Moreover, it has a void point to which anything can be linked.
 */

#ifndef _sparse_h
#define _sparse_h


#include <cstdio>
#include <iostream>
#include <string>

#include "Math/matrix.h"
#include "General/idata.h"

#include "Math/vector.h"

class model;


using namespace std;


enum storageT{  
	DNS ,			// Dense
	CSR ,			// Compressed Sparse Row Format
	CSC ,			// Compressed Sparse Column Format (DEFAULT)
	SCSC,			// symmetric CSC
	SSK ,			// Symmetric Skyline Format
	SSS	};			// Symmetric Sparse Skyline format



class sparseMatrix{
	
	
protected:

	
	// Numbering strategy for matrix pointers
	enum numberingT{
		NUM_FORTRAN,	// fortran type, starting from 1 
		NUM_C};			// C type,       starting from 0
	
	
	enum scalarT{
		SCALAR_SINGLE,
		SCALAR_DOUBLE,
		SCALAR_COMPLEX,
		SCALAR_DCOMPLEX};
	
	storageT			storage;		// the type of storage employed for the data
	dataT				datatype;		// data type real, single, double, complex
	numberingT			numbering;		// either Fortran or C style
	int					rows;			// number of rows
	int					cols;			// number of cols
	std::string			dumpfmt;		// the matrix format that is to be used if dumped to file
	
	double				*ddata;			// where the matrix entries go	
	bool				factorized;		// if LU are computed,
	double				condition;		// condition number of matrix
	
public:
	sparseMatrix(model &am, const char datatype='d');
	sparseMatrix();
	virtual ~sparseMatrix();
	
	//Functions to solve the eigenvalue problem
	virtual void computeEigenSystem(int nval, matrix& evect, double eval[], double evali[]=NULL, const char *which="LM"){};
	virtual void computeEigenValues(int nval, double eval[], double evali[]=NULL, const char *which="LM"){};
	
	
	// functions to manipulate a sparse matrix
	virtual void		addTo(double x, int row, int col)=0;
	virtual void		assemble(const matrix& block, const std::vector<int>& list)=0;
	virtual void		check()=0;
	virtual void		insert(double x, int row, int col)=0;
	inline  void		setConditionNumber(double cond) {condition = cond;}
	virtual void		setZero()=0;
	virtual void		transpose()=0;
	inline  void		setDumpFormat(std::string f) {dumpfmt = f;}
	virtual void		finalizeAssembly(){}
	
	// functions to get data or info from the matrix
	inline const int	getNColumns() const	{return cols;}
	inline int			getNRows()	const {return rows;}
	virtual int			getNTermsInColumn(const int col)=0;
	inline double		getConditionNumber() const {return condition;}
	inline storageT		getStorage() const {return storage;}
	inline double*		getData(){return ddata;}
	inline const double*	getData() const {return ddata;}
	
	// functions to check the status
	inline bool			isFactorized() const	{return factorized;}
	inline void			setNotFactorized(){factorized = false;}
	inline void			setFactorized(){factorized = true;}
	
	// functions to get information on the sparse matrix
	inline dataT		getDatatype() const	{return datatype;}
	virtual void		getDiagonal(longvector& diag)=0;
	virtual double		getEntry(int row, int col) const;
	virtual int			getNProfileTerms() const = 0;
	
	// functions to display information 
	virtual void		dump(std::string format="HB")=0;
	virtual void		info(ostream &of=cout);
	//virtual void		plot()=0;
	virtual void		print(FILE *fp=stdout);
	
	
	// operations
	virtual void		timesVector(longvector& b, longvector& Axb)=0;
	virtual void		timesDoubleArray(double b[], double Axb[])=0;
	virtual matrix		toFull();
	virtual void		convertToFortran()=0;
	virtual void		convertToC()=0;
};




#endif

