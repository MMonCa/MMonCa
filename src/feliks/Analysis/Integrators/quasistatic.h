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
 * static.h
 *
 * i. romero july 2000
 *
 * quasiStatic integrator
 */



#ifndef _quasistatic_h
#define _quasistatic_h


#include "Analysis/Integrators/integrator.h"
#include "Io/usercommand.h"

class assembler;
class node;
class Unode;
class commandLine;


class quasiStatic : public integrator
{

public:
	quasiStatic();
	quasiStatic(const commandLine &cl);
	virtual ~quasiStatic();
	
	virtual void		advanceDofset(vectorDofset &nd) const;
	virtual void		advanceDofset(scalarDofset &nd)	const;
	virtual void		advanceDofset(rotationDofset &nd) const;
	virtual void		advanceDofset(directorDofset &nd) const;
	
	virtual void		incrementDofset(vectorDofset &ds) const;
	virtual void		incrementDofset(scalarDofset &ds) const;
	virtual void		incrementDofset(rotationDofset &ds) const;
	virtual void		incrementDofset(directorDofset &DS) const;
	
	virtual void		startStep(vectorDofset &nd, const double dt){};
	virtual void		startStep(scalarDofset &nd, const double dt){};
	virtual void		startStep(rotationDofset &nd, const double dt){};
	virtual void		startStep(directorDofset &nd, const double dt){};


	

	void  info(ostream &of=cout);
	void  initialize(model &m, linearSolver &ls, linearSOE &theSOE, assembler& as);
	void  incrementNodeSolution(node &nd) const;
	void  setNodeDof(node &nd, int ndof, double x, double dt);
	
	virtual void startStep(model &m, const double dt);
};

#endif


