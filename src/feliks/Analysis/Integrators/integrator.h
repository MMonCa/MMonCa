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
 * integrator.h
 *
 * ignacio romero
 * june 2000, converted into C++ in dec 2005
 *
 *
 * the integrator class orchestrates the time advancing solution strategy.
 * There are different types of integrators: static, newmark, hht, ...
 * but all have to do the same things: give initial guesses to the dofs, 
 * advance them in the iterations and after the end of each time step.
 * They must be also able to go back and restart a time step solution.
 *
 * Many of the methods are pure virtual because they are specific of the
 * integrators and nothing can be said that applies to all.
 */

#ifndef _integrator_h
#define _integrator_h

#include <string>
#include <ostream>

#include "Model/model.h"
#include "Elements/element.h"
#include "Model/Node/node.h"
#include "Math/linearsoe.h"
#include "Math/LinearSolvers/linearsolver.h"
#include "Math/Sparse/sparse.h"


class assembler;

class integrator{
	
protected:
	std::string  name;
    bool         implicit;
	unsigned	 integrationOrder;
    double       eqtime;					// in (0,1], time at which force is evaluated
	bool		 nodallyBased;				// the updates are node by node
	bool		 variableMass;				// if true, recompute mass every time step
	bool		 spatialUpdates;			// if ture, rotational updates are spatial
	
	
public:		

	double       theTime;
	
	integrator();
	integrator(const commandLine& cl);
	virtual ~integrator();
	
	bool				check();
    static  void		scan(commandLine &cl, integrator **i);

	
	virtual void		advanceSolutionInTime(model &m, const double dt);		// advance variables from tn to tn+1
	void				advanceTimeCounters(const double dt);
	
	virtual void		advanceDofset(vectorDofset &nd)	const;
	virtual void		advanceDofset(scalarDofset &nd)	const;
	virtual void		advanceDofset(directorDofset &nd) const;
	virtual void		advanceDofset(rotationDofset &nd) const;
	
    virtual void		incrementDofset(vectorDofset &ds) const=0;
	virtual void		incrementDofset(scalarDofset &ds) const=0;
	virtual void		incrementDofset(directorDofset &ds) const{ cout << "not implemented incrementDofset" << endl;}
	virtual void		incrementDofset(rotationDofset &ds) const{ cout << "not implemented incrementDofset" << endl;}

	virtual void		startStep(vectorDofset &nd, const double dt){}
	virtual void		startStep(scalarDofset &nd, const double dt){}
	virtual void		startStep(rotationDofset &nd, const double dt){}
	virtual void		startStep(directorDofset &nd, const double dt){}
	
	
    virtual	void		incrementElementSolution(element &e) const;
    virtual void		incrementSolution(vector<node*> &nodes) const;	// increment solution after a step
	virtual void		incrementSolution(vector<element*> &e) const;
	
	
	virtual void		initialize(model &m, linearSolver &linsolver, linearSOE &theSOE, assembler &as);
    virtual void		info(ostream &of=cout);
    
	virtual void		startStep(model &m, const double dt)=0;			// called at the beginning of every t-step
    virtual void		setNodeDof(node &nd, int ndof, double x, double dt)=0;	// enforce a dof, and corresponding rates

	
    inline  std::string  getName() const	{return name;}
	inline  double       getEqTime()  const	{return eqtime;}
	inline	unsigned	 getIntegrationOrder() const {return integrationOrder;}
	inline  double       getTheTime() const	{return theTime;}
	inline  int          isImplicit() const	{return implicit;}
	inline	bool		 isNodallyBased() const		{return nodallyBased;}
};
  

#endif

