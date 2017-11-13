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
 *  petscmatrix.h
 *  feliks
 *
 *  Created by Ignacio Romero on 11/9/09.
 *  Copyright 2009 Universidad Politécnica de Madrid. All rights reserved.
 *
 */


#ifndef _petscmatrix_h
#define _petscmatrix_h


#ifdef WITHPETSC

#include "Math/Sparse/sparse.h"
#include "petscda.h"
#include "petscsnes.h"


class compressedMatrix;
class CSCMatrix;
class model;


class petscmatrix : public sparseMatrix
	{
		
	public:
		petscmatrix();
		petscmatrix(const model &m, const char* datatype = MATSEQAIJ);
		petscmatrix(const compressedMatrix& spm);
		petscmatrix(const petscmatrix& pm);
		~petscmatrix();
		bool createProfile(const model &m, const char* datatype = MATSEQAIJ);

		
		// functions to manipulate a sparse matrix
		void				addTo(double x, int row, int col);
		void				assemble(const matrix& block, const std::vector<int>& list);
		void				check(){}
		virtual void		finalizeAssembly();
		void				insert(double x, int row, int col);
		void				setZero();
		void				transpose();
		
		// functions to get information
		virtual int			getNTermsInColumn(const int col);
		virtual void		getDiagonal(longvector& diag);
		virtual double		getEntry(int row, int col) const;
		virtual int			getNProfileTerms() const;
		
		// functions to display information 
		virtual void		dump(std::string format="HB"){};
		virtual void		info(std::ostream &of=cout);
		virtual void		plot(){}
		virtual void		print(FILE *fp=stdout);
		
		
		// operations
		virtual void		timesVector(longvector& b, longvector& Axb);
		virtual void		timesDoubleArray(double b[], double Axb[]);
		virtual void		convertToFortran(){};
		virtual void		convertToC(){};
		
		Mat A;
	
	private:
		int *colstart;
		int *prow;
		int nnz;
		int rank;
		
		friend class petscsolver;
		friend class petscss;
		friend class eigenproblem;
		
	};


#endif

#endif


