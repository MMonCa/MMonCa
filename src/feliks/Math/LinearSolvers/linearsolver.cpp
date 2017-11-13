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
 *  linearsolver.cpp
 *
 *  Created by Ignacio Romero on Wed May 05 2004.
 *
 */

#include <cstring>


#include "Math/LinearSolvers/linearsolver.h"
#include "Math/linearsoe.h"
#include "Io/logger.h"
#include "Main/feliks.h"
#include "General/timer.h"
#include "Io/message.h"

#include "Math/LinearSolvers/cg.h"
#include "Math/LinearSolvers/direct.h"


/* this sets the default solver for linear equations, in accordance
   with the options set in feliks.h
*/
linearSolver :: linearSolver() :
	nprocessors(1),
	rtol(1e-8)
{
}



linearSolver :: linearSolver(const commandLine &cl) :
	nprocessors(1),
	rtol(1e-8)
{
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if       (uc.keyword() == "nprocessors") nprocessors = (int) uc.value();
	}
}



linearSolver :: linearSolver(const linearSolver& ls) :
	name(ls.name),
	rtol(ls.rtol),
	nprocessors(ls.nprocessors),
	storage(ls.storage)
{
}




linearSolver :: ~linearSolver()
{
}




void  linearSolver :: info(ostream &of)
{
	of  << "\n\n\n\n";
	printCentered(of, "S o l v e r   f o r   l i n e a r   e q u a t i o n s");
	of  << endl;
}





/*  scan the command line
*  "linsolver, type = [direct|cholesky|pcg]
*/
void linearSolver :: scan(const commandLine &cl, linearSolver **ls)
{
	linearSolver *ns=0;	
	
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if       ( uc.keyword() == "type" && uc.option() == "pcg")		ns = new cg(cl);
		else if  ( uc.keyword() == "type" && uc.option() == "direct")	ns = new direct(cl);
		else if  ( uc.keyword() == "type" && uc.option() == "cholesky")	ns = new direct(cl);
#ifdef WITHPETSC
		else if  ( uc.keyword() == "type" && uc.option() == "petsc")	ns = new petscsolver(cl);
#endif
	}
	
	if (ns != 0)
	{
		if (*ls != 0) delete *ls;
		*ls = ns;
	}
}




/* this should, in the future, choose among available solvers.
Returns 1 if ok, -1 if cannot solve system.
*/
bool linearSolver :: solve(linearSOE &theSystem, double &cputime, double &mflops)
{
    timer  t;
    double mlsec;
    bool    ret=true;

    if (theSystem.getDimension() == 0) return 1;

    t.start();
	this->solve(theSystem, mflops);
	t.stop();
	
    mlsec   = t.miliseconds();
	cputime = 0.001*mlsec;

    DebugMessage("LinearSOE solved");
    return ret;
}


/* this should, in the future, choose among available solvers.
 Returns 1 if ok, -1 if cannot solve system.
 */
bool linearSolver :: solve(double &cputime, double &mflops)
{
    timer  t;
    double mlsec;
    bool    ret=true;
	
    if (theLinearSOE.getDimension() == 0) return 1;
	
    t.start();
	this->solve(mflops);
	t.stop();
	
    mlsec   = t.miliseconds();
	cputime = 0.001*mlsec;
	
    DebugMessage("LinearSOE solved");
    return ret;
}


/* this should, in the future, choose among available solvers.
 Returns 1 if ok, -1 if cannot solve system.
 */
bool linearSolver :: solve(double &mflops)
{
	return this->solve(theLinearSOE, mflops);
}



