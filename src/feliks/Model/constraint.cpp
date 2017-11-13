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
 *  constraint.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 26/6/09.
 *  Copyright 2009 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include "Model/constraint.h"

#include <exception>
#include <iostream>
#include <typeinfo>

#include "Analysis/Loading/pointload.h"
#include "Analysis/Integrators/integrator.h"
#include "Analysis/Propfactors/propfactor.h"

#include "Geometry/Solidmodel/BRep.h"
#include "Math/Topology/topology.h"
#include "Model/model.h"
#include "Model/Parts/body.h"

#include "Io/io.h"


#include "boost/foreach.hpp"

#define DEBUG_SCANNER 0

using namespace feliks;


constraint :: constraint()
:
    factorlabel(0),
    fc(0)
{
    
}


constraint :: constraint(const commandLine& cl) :
factorlabel(0),
fc(0)
{
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if      (uc.keyword() == "scaling") factorlabel  = uc.integerValue();
	}
}




void constraint :: print(std::ostream& of) const
{
	of                      << "\n       scaling      : " << fc->getLabel();
	if (fc->getLabel() == 0) of << "  (default)";
}





/* the boundary conditions has two formats. 
 
 bc, nodeset = name [, ux = value] [, uy = value] [, uz = value]  [, p = value] 
 [,rotx = value] [, roty = value] [,rotz = value]
 [,dirx = value] [, diry = value]
 [, scaling = value]
 
 bc, body = name, dimension = <int>, label = <int>, [, ux = value] [, uy = value] [, uz = value]  [, p = value] 
 [,rotx = value] [, roty = value] [,rotz = value]
 [,dirx = value] [, diry = value]
 [, scaling = value]
 
 
 */
void constraint :: scan(const commandLine &cl, ifstream &meshfile, model &m)
{
	if (cl.size() < 2)
    {
        logger::mainlog  << "\n Obsolete call to BC block";
        logger::warnings << "\n Obsolete call to BC block";
		//scanBCblock(cl, meshfile, m);
    }
	
	else
	{
		int k;
		for (k=1; k< cl.size(); k++)
		{
			const usercommand  &uc = cl[k];
			
			if (uc.keyword() == "nodeset")
			{
				constraint *sc = new nodesetConstraint(cl);
				m.add(*sc);
				break;
			}
			
			if (uc.keyword() == "body")
            {
                constraint *sc = new manifoldConstraint(cl);
                m.add(*sc);
                break;
            }
		}
		if (k == cl.size()) logger::mainlog << "\n Error in scan constraint" ;
	}
}






/* the format for the boundary conditions is
 bc
 nodelabel1  dofindex value scalingfactor
 nodelabel2  dofindex value scalingfactor
 ...
 [blank line]
 *
 void constraint :: scanBCblock(const commandLine &cl, ifstream &meshfile, model &m)
 {
 string oneLine;
 
 // scan boundary conditions definitions  until blank line or EOF
 while ( commandLine::scanCommandLine(meshfile, oneLine) && oneLine.size() > 0)
 {
 // set the default values, so it does not take the last one if missing
 int vertexlabel  =  0;
 int dofnumber    = -1;
 double bcvalue   =  0.0;
 int factorlabel  =  0;
 
 stringstream strm(oneLine);
 strm >> vertexlabel >> dofnumber >> bcvalue >> factorlabel;
 DOFconstraint *bc = new DOFconstraint(vertexlabel, dofnumber-1, bcvalue, factorlabel);
 if (DEBUG_SCANNER == 1) cout << "adding constraint on node: " << vertexlabel << ", dof " << dofnumber << endl;
 if (bc != NULL) m.add(*bc);
 }
 }
 */





DOFconstraint :: DOFconstraint(const int ndlabel, const int dofn, const double F, const int flabel):
_nodelabel(ndlabel), 
_dof(dofn), 
_value(F)
{
    factorlabel = flabel;
}





list<DOFconstraint> DOFconstraint :: convertToDOFconstraints()
{
	list<DOFconstraint> l;
	l.push_back( DOFconstraint( _nodelabel , _dof, _value, factorlabel) );	
	
	return l;
}



void DOFconstraint :: imposeConstrainedDOFs(const double time, const double dt, const bool isImplicit, integrator& theIntegrator)
{	
	// only need to update dofs which are not zero
	if (_value != 0.0)
	{		
		// here's the different strategy for implicit/explicit methods
		double f   = fc->eval(time, _nd->coordinates() );
		if ( isImplicit )   
		{
			double df  = f - fc->eval(time - dt, _nd->coordinates() );
			_nd->setDofIncrement(_dof, _value*df);
		}
		else
			theIntegrator.setNodeDof(*_nd, _dof, _value*f, dt);
	}	
}


void DOFconstraint :: initialize(model& m)
{
	fc  = &(factorcombo::getPropfactorCombinationFromLabel(factorlabel));
	_nd = &(m.getNode(_nodelabel) );
}



void DOFconstraint ::	markConstrainedDofs()
{
	_nd->constrainDof(_dof);	
}


/* prints the information of a point load.
 * The local degree of freedom is incremented by one, since it is more natural for
 * the user
 */
void DOFconstraint :: print(ostream &of) const
{
	of  << "\n" 
	<< setw(4) << _nodelabel << "    "
	<< setw(2) << _dof+1     << "   "
	<< setw(9) << _value     << "   "
	<< setw(3) << factorlabel << flush;
}




manifoldConstraint :: manifoldConstraint(const commandLine& cl)
:
nodesetConstraint(cl),
theBRepConstrainedCell(0)
{        
    for (int k=1; k< cl.size(); k++)
    {
        const usercommand  &uc = cl[k];
        
        if      (uc.keyword() == "body")         bodyname   = uc.option();
        else if (uc.keyword() == "dimension")    dimension  = uc.integerValue();
        else if (uc.keyword() == "label")        mnfldLabel = uc.integerValue();
    }
}


manifoldConstraint :: ~manifoldConstraint()
{
    
}




void manifoldConstraint :: initialize(model& m)
{
}




void manifoldConstraint :: print(std::ostream& of) const
{
    of << "\n\n Constraint on a body boundary";
    of << "\n Name of the body   : " << bodyname;
    of << "\n Boundary dimension : " << dimension;
    of << "\n Boundary label     : " << mnfldLabel;
    
    if (uconstrained[0]) of << "\n       ux imposed   : " << u[0];
	if (uconstrained[1]) of << "\n       uy imposed   : " << u[1];
	if (uconstrained[2]) of << "\n       uz imposed   : " << u[2];
	if (pconstrained)    of << "\n       p  imposed   : " << p;
	if (rconstrained[0]) of << "\n     rotx imposed   : " << r[0];
	if (rconstrained[1]) of << "\n     roty imposed   : " << r[1];
	if (rconstrained[2]) of << "\n     rotz imposed   : " << r[2];
	if (dconstrained[0]) of << "\n     dirx imposed   : " << d[0];
	if (dconstrained[1]) of << "\n     diry imposed   : " << d[1];
	
    constraint::print(of);
}


nodesetConstraint :: nodesetConstraint(nodeset* Ndst, int uconst, double uvalue) :
constraint(),
theNodeset(Ndst),
theNodesetName(Ndst->getName())
{	
    uconstrained[0] = uconstrained[1] = uconstrained[2] = false;
	pconstrained    = false;
	rconstrained[0] = rconstrained[1] = rconstrained[2] = false;
	dconstrained[0] = dconstrained[1] = false;
    
    
    if (uconst == 0) uconstrained[0] = true;
    else if (uconst == 1) 
    {
        uconstrained[1] = true;
        u[1] = 0.0;
    }    
    else if (uconst == 2) 
    {
        uconstrained[2] = true;
        u[2] = 0.0;   
    }
 }


nodesetConstraint :: nodesetConstraint(const commandLine &cl) :
constraint(cl),
theNodeset(0),
theNodesetName()
{	
    uconstrained[0] = uconstrained[1] = uconstrained[2] = false;
	pconstrained    = false;
	rconstrained[0] = rconstrained[1] = rconstrained[2] = false;
	dconstrained[0] = dconstrained[1] = false;
    
    
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		
		if      (uc.keyword() == "nodeset")	theNodesetName = uc.option();
        else if	(uc.keyword() ==  "ux")		u[0]   = uc.value(), uconstrained[0] = true;
		else if	(uc.keyword() ==  "uy")		u[1]   = uc.value(), uconstrained[1] = true;
		else if	(uc.keyword() ==  "uz")		u[2]   = uc.value(), uconstrained[2] = true;
		
		else if	(uc.keyword() ==  "p")		p      = uc.value(), pconstrained    = true;
		
		else if	(uc.keyword() ==  "rotx")	r[0]   = uc.value(), rconstrained[0] = true;
		else if	(uc.keyword() ==  "roty")	r[1]   = uc.value(), rconstrained[1] = true;
		else if	(uc.keyword() ==  "rotz")	r[2]   = uc.value(), rconstrained[2] = true;
		
		else if	(uc.keyword() ==  "dirx")	d[0]   = uc.value(), dconstrained[0] = true;
		else if	(uc.keyword() ==  "diry")	d[1]   = uc.value(), dconstrained[1] = true;
	}
}



list<DOFconstraint> nodesetConstraint :: convertToDOFconstraints()
{
	list<DOFconstraint> l;
	
	set<node*>::iterator ni = theNodeset->begin();
	while ( ni != theNodeset->end() )
	{
		if (uconstrained[0]) l.push_back( DOFconstraint( (**ni).getLabel() , 0, u[0], fc->getLabel() ) );
		if (uconstrained[1]) l.push_back( DOFconstraint( (**ni).getLabel() , 1, u[1], fc->getLabel() ) );
		if (uconstrained[2]) l.push_back( DOFconstraint( (**ni).getLabel() , 2, u[2], fc->getLabel() ) );
		
		++ni;
	}
	
	return l;	
}



void nodesetConstraint :: initialize(model& m)
{
	fc         = &(factorcombo::getPropfactorCombinationFromLabel(factorlabel));
	theNodeset = &(m.getNodeset(theNodesetName));
}




void nodesetConstraint :: imposeConstrainedDOFs(const double time, const double dt, const bool isImplicit,integrator& theIntegrator)
{
	// only need to update dofs which are not zero
	for (int i=0; i<3; i++)
	{
        if (uconstrained[i] && u[i] == 0.0) // MMM Created just for controlvolumes, no matter the value, we have to update every time step the boundary conditions, as new nodes are generated.
		{
            set<node*>::iterator iter = theNodeset->begin();
            while (iter != theNodeset->end() )
			{
				node& nd = **iter;
				// here's the different strategy for implicit/explicit methods
				if ( !isImplicit )   
					theIntegrator.setNodeDof(nd, i, 0, dt);
                
                ++iter;
			}
   
        }
		if (uconstrained[i] && u[i] != 0.0)
		{		
            set<node*>::iterator iter = theNodeset->begin();
            while (iter != theNodeset->end() )
			{
				node& nd = **iter;
				
				// here's the different strategy for implicit/explicit methods
				double f   = fc->eval(time, nd.coordinates() );
				if ( isImplicit )   
				{
					double df  = f - fc->eval(time - dt, nd.coordinates() );
					nd.setDofIncrement( i , u[i]*df);
				}
				else
					theIntegrator.setNodeDof(nd, i, u[i]*f, dt);
                
                ++iter;
			}
		}
	}	

    if (pconstrained && p != 0.0)
    {		
        set<node*>::iterator iter = theNodeset->begin();
        while (iter != theNodeset->end() )
        {
            node& nd = **iter;
            
            // here's the different strategy for implicit/explicit methods
            double f   = fc->eval(time, nd.coordinates() );
            if ( isImplicit )   
            {
                double df  = f - fc->eval(time - dt, nd.coordinates() );
                nd.setDofIncrement( 3 , p*df);
            }
            else
                theIntegrator.setNodeDof(nd, 3, p*f, dt);
            
            ++iter;
        }
    }

	// only need to update dofs which are not zero
	for (int i=0; i<3; i++)
	{
		if (rconstrained[i] && r[i] != 0.0)
		{		
            set<node*>::iterator iter = theNodeset->begin();
            while (iter != theNodeset->end() )
			{
				node& nd = **iter;
				
				// here's the different strategy for implicit/explicit methods
				double f   = fc->eval(time, nd.coordinates() );
				if ( isImplicit )   
				{
					double df  = f - fc->eval(time - dt, nd.coordinates() );
					nd.setDofIncrement( i+3 , r[i]*df);
				}
				else
					theIntegrator.setNodeDof(nd, i+3, r[i]*f, dt);
                
                ++iter;
			}
		}
	}	


	for (int i=0; i<2; i++)
	{
		if (dconstrained[i] && d[i] != 0.0)
		{		
            set<node*>::iterator iter = theNodeset->begin();
            while (iter != theNodeset->end() )
			{
				node& nd = **iter;
				
				// here's the different strategy for implicit/explicit methods
				double f   = fc->eval(time, nd.coordinates() );
				if ( isImplicit )   
				{
					double df  = f - fc->eval(time - dt, nd.coordinates() );
					nd.setDofIncrement( i+3 , d[i]*df);
				}
				else
					theIntegrator.setNodeDof(nd, i+3, d[i]*f, dt);
                
                ++iter;
			}
		}
	}	
}



void nodesetConstraint :: markConstrainedDofs()
{
	for (int i=0; i<3; i++)
    {
		if (uconstrained[i])
        {
			BOOST_FOREACH(node* n, *theNodeset) n->constrainDof(i);
        }
    }
    
    if (pconstrained)
    {
        BOOST_FOREACH(node* n, *theNodeset)
        {
            n->constrainDof(3);
        }
    }
    
	for (int i=0; i<3; i++)
	{
        if (rconstrained[i])
        {
			BOOST_FOREACH(node* n, *theNodeset)
                n->constrainDof(3+i);
        }
    }
    
	for (int i=0; i<2; i++)
    {
		if (dconstrained[i])
        {
			BOOST_FOREACH(node* n, *theNodeset)
                n->constrainDof(3+i);
        }
    }
}



void nodesetConstraint :: print(std::ostream& of) const
{
	of << "\n\n Nodeset constraint on nodeset " << theNodesetName;
    if (uconstrained[0]) of << "\n       ux imposed   : " << u[0];
	if (uconstrained[1]) of << "\n       uy imposed   : " << u[1];
	if (uconstrained[2]) of << "\n       uz imposed   : " << u[2];
	if (pconstrained)    of << "\n       p  imposed   : " << p;
	if (rconstrained[0]) of << "\n     rotx imposed   : " << r[0];
	if (rconstrained[1]) of << "\n     roty imposed   : " << r[1];
	if (rconstrained[2]) of << "\n     rotz imposed   : " << r[2];
	if (dconstrained[0]) of << "\n     dirx imposed   : " << d[0];
	if (dconstrained[1]) of << "\n     diry imposed   : " << d[1];
	
    constraint::print(of);
}



