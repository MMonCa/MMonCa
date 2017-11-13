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
* ignacio romero
 * november 1999, translated to c++ in march 2006
 *
 * defines the properties of the analysis. This is the main data type in FELIKS, holding
 * all the information of the finite element analysis.  *
 */


#ifndef _analysis_h
#define _analysis_h

#include <string>
#include <memory>
#include <vector>
#include <list>

#include "Analysis/Stepsolvers/newtonraphson.h"
#include "Analysis/Propfactors/propfactor.h"
#include "Analysis/energies.h"
#include "Analysis/Integrators/integrator.h"
#include "Analysis/Stepsolvers/stepsolver.h"
#include "Analysis/TScontrol/tscontrol.h"
#include "Elements/element.h"
#include "Io/usercommand.h"
#include "Math/linearsoe.h"
#include "Math/Sparse/sparse.h"
#include "Math/LinearSolvers/linearsolver.h"
#include "Model/model.h"


class loading;
class constraint;

namespace feliks
{
    class mmoncaInterface;
}


enum analysisTypeT{
	FEANALYSIS_UNDEFINED,
	FEANALYSIS_STATIC,
	FEANALYSIS_TRANSIENT,
	FEANALYSIS_SPECTRAL,
	FEANALYSIS_HARMONIC,
	FEANALYSIS_QMODAL,
	FEANALYSIS_ASYNCHRONOUS,
	FEANALYSIS_STAGGERED};




class FEanalysis 
{
		
public:
	FEanalysis();
	virtual					~FEanalysis();
	bool					check(ostream &of=cout);
	virtual  void			initialize(){}
	virtual  bool			specificCheck(ostream &of=cout) {return true;}
	void					setInteractive();
	
	// functions to retrieve members of the analysis
	inline  const model&	getMesh()	const			{return *theModel;}
	inline  linearSOE&		getTheLinearSOE() const	{return *theLinearSOE;}
	inline  linearSolver&	getTheLinearSolver() const{return *theLinearSolver;}
	inline  tscontrol&		getTheTSControl() const	{return *theTSControl;}
	inline  stepsolver&		getTheStepSolver()		{return *theStepsolver;}
	
	
	// analysis information
	void					finalMessage(ostream &of=cout, bool ok=true);
	virtual void			info(ostream &os=cout);
	virtual void			infoSummary(ostream &os=cout);
	
	
	/* function to operate on the analysis, obtaining solutions for the linearized
	 set of equations. SolveGlobalSystem solves the system K du= R, where K and R 
	 must be already formed. SolveNewtonRaphsonIteration does the same thing, but
	 first constructs K and R. Both functions update the nodal solution and zero
	 the residual afterwards. */
	virtual bool            solve();
    virtual bool            specificSolve()=0;

	
protected:

	model*                          theModel;
	linearSOE*                      theLinearSOE;			// a pointer, not allocated	
	stepsolver*                     theStepsolver;          // the (nonlinear) solver of every step
	linearSolver*                   theLinearSolver;		// the solver of the system Ax=b
	tscontrol*                      theTSControl;			// time step controller
    
	std::auto_ptr<sparseMatrix>     mass;                   // the mass matrix
	std::auto_ptr<postprocessor>    thePost;
	energies                        theEnergies;			// energy, momenta, errors & more

	size_t                          dofs_per_node;          // maximum number of dofs/node
	bool                            interactive;
	
	virtual	void                    createFEdata(const feliks::mmoncaInterface& d);
	friend class                    newtonraphson;
	friend class                    explicitss;
	friend class                    interactiveFEanalysis;
	
    
    
private:
		
	// prevent analysis copy by declaring it private, and then leaving the function empty
	FEanalysis(const FEanalysis& rhs);
	FEanalysis& operator=(const FEanalysis& rhs);
};



#endif

