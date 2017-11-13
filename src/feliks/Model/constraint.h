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
/* constraint.h
 *
 * single point degree-of-freedom constraints
 *
 * ignacio romero
 *
 * the data is identical to that of pointLoads, so we reuse all the functions. 
 * Other options would be
 *  i) to rewrite everything identical to pointload.h and pointload.c
 * ii) to eliminate the constraint data type and just consider pointloads over
 * constrained dofs differently.
 *
 * I have chosen to have a different datatype for clarity but identical implementation.
 * Data types and functions are then just alias to the ones in pointLoad.
 */


#ifndef _constraint_h
#define _constraint_h
	
#include "Analysis/Loading/pointload.h"
#include "Model/Sets/nodeset.h"

#include <iostream>
#include <string>
#include <list>
#include <set>


class commandline;
class faceset;
class integrator;
class model;
class nodeset;
class DOFconstraint;



namespace feliks
{
    namespace topology
    {
        class cell;
    }
}



class constraint
{
	public:
                                constraint(const commandLine& uc);
		virtual                 ~constraint(){};
	
		static	void            scan(const commandLine &cl, ifstream &meshfile, model &m);
        factorcombo&            getScalingFactor() const;
	
	
		virtual list<DOFconstraint> convertToDOFconstraints()=0;	
		virtual void            imposeConstrainedDOFs(const double t, const double dt, const bool isImpl, integrator& i)=0;
		virtual void            initialize(model& m)=0;
		virtual void            markConstrainedDofs()=0;
		virtual void            print(std::ostream& of=std::cout) const;
	
    
    protected:
        int                     factorlabel;
        factorcombo             *fc;
                                constraint();
};


inline factorcombo&             constraint::getScalingFactor() const {return *fc;}





class DOFconstraint : public constraint
	{
		
	public:
                                DOFconstraint(const int ndlabel, const int dofn, const double F, const int flabel);

        int                     getDof() const;
        node&                   getNode() const;
        int                     getNodeLabel() const;
        double                  getValue()	const;
		
        std::list<DOFconstraint>     convertToDOFconstraints();
		virtual void            initialize(model& m);
		void                    imposeConstrainedDOFs(const double t, const double dt, const bool isImplicit, integrator& i);
		void                    markConstrainedDofs();
		void                    print(std::ostream& of=std::cout) const;

		
	private:
		int                     _dof;				// local dof of the node. Starts from 0;
		double                  _value;				// the reference value of the load.
		int                     _nodelabel;			// label of the node where the load is applied
		node                    *_nd;
		
		friend class model;
	};


inline int      DOFconstraint::getDof() const       {return _dof;}
inline node&    DOFconstraint::getNode() const      {return *_nd;}
inline int      DOFconstraint::getNodeLabel() const	{return _nodelabel;}
inline double   DOFconstraint::getValue() const		{return _value;}




class nodesetConstraint : public constraint
	{
	
	public:
                                nodesetConstraint(const commandLine &cl);
                                nodesetConstraint(nodeset* Ndt, int uconst, double uvalue);                        
                                virtual ~nodesetConstraint(){}
		
		
		virtual list<DOFconstraint> convertToDOFconstraints();
		virtual void            imposeConstrainedDOFs(const double t, const double dt, const bool isImplicit, integrator& i);
		virtual void            initialize(model& m);
		virtual void            markConstrainedDofs();
		virtual void            print(std::ostream& of=std::cout) const;
		
	
	protected:
		nodeset*                theNodeset;
		string                  theNodesetName;
        bool                    uconstrained[3], pconstrained, rconstrained[3], dconstrained[2];
        double                  u[3], p, r[3], d[2];
	};





class manifoldConstraint : public nodesetConstraint
{
    
public:
                                manifoldConstraint(const commandLine& cl);
    virtual                     ~manifoldConstraint();
    
    virtual void                initialize(model& m);
    virtual void                print(std::ostream& of=std::cout) const;


private:
    
    feliks::topology::cell*       theBRepConstrainedCell;
    string                      bodyname;
    int                         dimension, mnfldLabel;
    nodeset                     privateNodeset;
};



#endif
