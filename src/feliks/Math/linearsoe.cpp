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
/* linearsoe.cpp
 *
 * linear system of equations
 *
 * i. romero
 * july 2000
 *
 * functions for systems of the form Ax=b
 */

#include "Math/linearsoe.h"
#include "Math/Sparse/compressed.h"

#include <string>

//#include </usr/include/complex.h>
#include <cstdio>
#include "Io/logger.h"
#include "Math/vector.h"
#include "Math/Sparse/sparse.h"
#include "Io/message.h"
#include "Model/model.h"



linearSOE :: linearSOE() :
	Cprec(NULL),
	Rprec(NULL),
	A(0),
	b(0),
	deleteMembersWhenFinished(false)
{}


// allocates the memory space for a new system of equations used to solve
// the problem. By default, the matrix is sparse.
linearSOE :: linearSOE(model &m, const storageT st, const string& dumpfmt) :
	Cprec(NULL),
	Rprec(NULL),
	A(0),
	b(0),
	deleteMembersWhenFinished(true)
{
	dimension = (int) m.getNDofs();
	if     (st == SCSC) 	A = new SCSCMatrix(m);
	else					A = new CSCMatrix(m);
	A->setDumpFormat(dumpfmt);
	
	b = new longvector(dimension);
}



linearSOE :: linearSOE(sparseMatrix &Ac, longvector& bc) :
	deleteMembersWhenFinished(false)
{
    dimension = Ac.getNColumns();
    A         = &Ac;
    b         = &bc;
    Cprec     = NULL;
    Rprec     = NULL;
}




linearSOE :: ~linearSOE()
{		
	if (Cprec != NULL) delete [] Cprec;
	if (Rprec != NULL) delete [] Rprec;
	if (deleteMembersWhenFinished)
	{
		delete A;
		delete b;
	}
}





void linearSOE ::  info(ostream &of)
{
	of << endl << endl << endl;
	printCentered(of, "L i n e a r   s y s t e m   o f   e q u a t i o n s");
	extern int global_complexdata;
	
	of << "\n\n System of equations and tangent:";
	of << "\n\tNumber of equations        : " <<  dimension;
	if (global_complexdata == 1) of << " (complex)";
	A->info(of);
}



void linearSOE :: setLinks(sparseMatrix &Ac, longvector& bc)
{
    dimension = Ac.getNColumns();
    A         = &Ac;
    b         = &bc;
    Cprec     = 0;
    Rprec     = 0;
}



void linearSOE :: zeroA()
{
	A->setZero();
}


void linearSOE :: zeroB()
{
	b->setZero();
}
