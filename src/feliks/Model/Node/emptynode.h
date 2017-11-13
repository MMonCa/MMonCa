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
 *  emptynode.h
 *  feliks
 *
 *  Created by Ignacio Romero on 8/30/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */


#ifndef _emptynode_h
#define _emptynode_h


#include "node.h"
#include "Math/tensor.h"


class integrator;
class longvector;

class emptyNode : public node
{
	
public:
	emptyNode(const int label_, const blue::ivector& coor) : node(label_, coor) {}
	emptyNode(const int label_) : node(label_) {}
	virtual ~emptyNode(){}
	
	
	virtual	void	advanceIntegration(const integrator& integ, const double dt){}
	virtual	void	incrementSolution(const integrator& integ){}
	virtual void	computeExplicitAcceleration(){}
	virtual void	computeExplicitVelocity(){}
	virtual void	computeExplicitEulerianAcceleration(){}
	
	// functions to get information about the node
	virtual	int		getNDofs() const				{return 0;}
	virtual nodetypeT	getNodetype() const			{return _EMPTYnode;}
	virtual bool	isOrphan() const				{return false;}
		
	
	// functions to manipulate dofs
	virtual void		recoverFromBackupDofs(){}
	virtual void		setDofIncrement(const int local, const double v){};
	virtual void		setInitialDof(const int local, const double v){};
	virtual void		setDeltaDofsToZero(){};
	
	
	// functions to assign and retrieve dof/equation mapping ID. Local id start from 0
	virtual void		constrainDof(const int local){};
	virtual int			getID (const int local) const{cout << "Error in getID for emptyNode" << endl; return 0;};
	
	
	// functions to retrieve information about free/constrained dofs
	virtual int			getNConstrainedDofs() const	{return 0;};
	virtual bool		isDofConstrained(const int local) const{return true;}
	
	
	// functions to manipulate the forces
	virtual void		assembleForceDof(const int dof, const double f){}
	virtual void		forceReset(){};
	
	virtual void		localizeSolutionIncrement(const longvector& delta,  partitionDofs pt=PARTITION_ALL){};
	
	
	// functions to display or print information about the node
	// and its dofs (1,2,3,...)
	virtual void		printID(ostream &of=cout) const{};
	virtual void		printDOFs(ostream &of=cout) const{};
	
};

#endif


