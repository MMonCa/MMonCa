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

//
//  poorbody.cpp
//  feliks
//
//  Created by Romero Ignacio on 4/2/12.
//  Copyright (c) 2012 Universidad Politécnica de Madrid. All rights reserved.
//

#include <iostream>
#include <utility>
#include <sstream>
#include <boost/foreach.hpp>
#include "poorbody.h"

#include "Elements/element.h"
#include "Elements/evalspot.h"

#include "General/feliksutil.h"
#include "Io/io.h"
#include "Model/model.h"
#include "Model/Parts/modelpart.h"
#include "Model/Node/node.h"
#include "Model/Node/emptynode.h"

#ifdef WITHTBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_for_each.h"
#endif

#define DEBUG_SCANNER 0


using namespace blue;






poorbody :: poorbody(const std::string& name, const std::string& geo, model& m, const size_t dim)
:
    modelpart(name, m),
    allNodes(),
    externalNodes(),
    remeshafterthesesteps(100000),
    _geometricDescription(geo),
    _isParsed(false),
    _hasVariableMass(false),
    _massComputed(false),
    _externalNodesFound(false),
    _dimension(dim),
    mass(0.0),
    coordinationNumber(0)
{
 	// an empty nodeset with the external nodes
	externalNodes.setName(getName() + "_external_nodes" );
	externalNodes.declareInternal();	
	m.add(externalNodes);
	
	// an empty nodeset with all the nodes of the body
	allNodes.setName( getName() + "_allnodes" );
	allNodes.declareInternal();
	m.add(allNodes);
}




poorbody :: poorbody(const std::string& name, model& m)
:
    modelpart(name, m),
    allNodes(),
    externalNodes(),
    remeshafterthesesteps(100000),
    _geometricDescription("undefined"),
    _isParsed(false),
    _hasVariableMass(false),
    _massComputed(false),
    _externalNodesFound(false),
    _dimension(3),
    mass(0.0),
    coordinationNumber(0)
{
 	// an empty nodeset with the external nodes
	externalNodes.setName(getName() + "_external_nodes" );
	externalNodes.declareInternal();	
	m.add(externalNodes);
	
	// an empty nodeset with all the nodes of the body
	allNodes.setName( getName() + "_allnodes" );
	allNodes.declareInternal();
	m.add(allNodes);
}





poorbody :: ~poorbody()
{
    // the externalNodes and allNodes are internal nodesets. This means that
	// the user has not defined them but rather are created simultaneously with the
	// body. They belong to the body and have not meaning without it, so we eliminate them
	linkedModel.remove(externalNodes);
	linkedModel.remove(allNodes);
    
    // remove the nodes the vectors, and later, the vectors themselves.
    BOOST_FOREACH(element*  e, elements)     delete e;
	elements.clear();
    
    BOOST_FOREACH(node*     n, nodes)  delete n;
    nodes.clear();
    
    //BOOST_FOREACH(evalspot* v, evalspots) delete v;
	evalspots.clear();
    
    externalNodes.clear();
    allNodes.clear();
}






#ifdef WITHTBB
// data structure for tbb threaded version
class tbbAdvanceNodesInTime{
public:
    tbbAdvanceNodesInTime( vector<node*>& nds, const integrator &integ, const double dt) 
    : my_nodes(&nds), my_integ(&integ), my_dt(dt){}
    
    void operator()( const tbb::blocked_range<size_t>& r) const
    {
        vector<node*> &nds( *my_nodes );
        for (size_t a=r.begin(); a!=r.end(); a++)
            nds[a]->advanceIntegration(*my_integ, my_dt);
        
    }
    
private:
    vector<node*> *const my_nodes;
    integrator     const *my_integ;
    const double   my_dt;
};
#endif




class advanceInTimeWorker
{
private:
    vector<node*> *const my_nodes;
    integrator     const *my_integ;
    const double   my_dt;
    
public:
    advanceInTimeWorker(vector<node*>& nds, const integrator &integ, const double dt) 
    : my_nodes(&nds), my_integ(&integ), my_dt(dt){}
    
    void operator()(node* nd) const { nd->advanceIntegration(*my_integ, my_dt);}
};




// when a solution has been obtained, all the body variables and the ones
// in its elements must be updated accordingly. 
void poorbody :: advanceInTime(const integrator& integ, const double dt)
{
	if ( !isActive() ) return;
    
    updateCurrentState(dofset::tn1);
    feliks_for_each(elements.begin(),  elements.end(),  mem_fun(&element::commitCurrentState) );
    feliks_for_each(evalspots.begin(), evalspots.end(), mem_fun(&evalspot::commitCurrentState) );
    feliks_for_each(nodes.begin(),     nodes.end(),     advanceInTimeWorker(nodes, integ, dt));
    
    
    /*
     // advance the nodal variables III
     #ifdef WITHTBB
     
     tbb::parallel_for(tbb::blocked_range<size_t>(0,nodes.size()), 
     tbbAdvanceNodesInTime(nodes, integ, dt));
     
     #else
     
     vector<node*>::iterator iter = nodes.begin();
     while ( iter != nodes.end() )
     {
     (*iter)->advanceIntegration(integ, dt);
     ++iter;
     }   
     
     #endif
     */
	
	if (_hasVariableMass)	computeMass();	
}




bool poorbody :: check() const
{
	for_each(elements.begin(), elements.end(), mem_fun(&element::check) );
    
 	return true;
}



#ifdef WITHTBB
// data structure for tbb threaded version
class tbbComputeExplicitAcceleration
{
public:
	tbbComputeExplicitAcceleration(vector<node*>& nodes) : my_nodes(&nodes){}
	
	void operator()( const tbb::blocked_range<size_t>& r) const
	{
		vector<node*> &nn( *my_nodes );		
		for (size_t a=r.begin(); a != r.end(); ++a)	
			nn[a]->computeExplicitAcceleration();
	}
	
private:
	vector<node*> *const my_nodes;
};
#endif




void poorbody :: computeExplicitAcceleration()
{
#ifdef WITHTBB
	tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
					  tbbComputeExplicitAcceleration(nodes));
#else
	for_each(nodes.begin(), nodes.end(), mem_fun(&node::computeExplicitAcceleration) );	
	
#endif
	
}



#ifdef WITHTBB
// data structure for tbb threaded version
class tbbcomputeExplicitEulerianAcceleration{
public:
	tbbcomputeExplicitEulerianAcceleration( vector<node*>& nodes) : my_nodes(&nodes){}
	
	void operator()( const tbb::blocked_range<size_t>& r) const
	{
		vector<node*> &nn( *my_nodes );		
		for (size_t a=r.begin(); a != r.end(); a++)	
			nn[a]->computeExplicitEulerianAcceleration();
	}
	
private:
	vector<node*> *const my_nodes;
};
#endif



void poorbody :: computeExplicitEulerianAcceleration()
{
#ifdef WITHTBB
	tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
					  tbbcomputeExplicitEulerianAcceleration(nodes),
					  tbb::auto_partitioner() );
#else
	for_each(nodes.begin(), nodes.end(), mem_fun(&node::computeExplicitEulerianAcceleration) );	
	
#endif
	
}

#ifdef WITHTBB
// data structure for tbb threaded version
class tbbComputeExplicitVelocity{
public:
	tbbComputeExplicitVelocity( vector<node*>& nodes) : my_nodes(&nodes){}
	
	void operator()( const tbb::blocked_range<size_t>& r) const
	{
		vector<node*> &nn( *my_nodes );		
		for (size_t a=r.begin(); a != r.end(); a++)	
			nn[a]->computeExplicitVelocity();
	}
	
private:
	vector<node*> *const my_nodes;
};
#endif


void poorbody :: computeExplicitVelocity()
{
#ifdef WITHTBB
	tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
					  tbbComputeExplicitVelocity(nodes),
					  tbb::auto_partitioner() );
#else
	for_each(nodes.begin(), nodes.end(), mem_fun(&node::computeExplicitVelocity) );	
	
#endif
	
}


void poorbody :: computeExternalNodesAndFaces()
{
	if (_externalNodesFound) return;
	
	elset tmp(linkedModel);
	
	// create a temporary elset with all the elements of the body
    BOOST_FOREACH(element* e, elements)
    tmp.add(*e);
    
	// compute the external nodes of the elset and copy them to externalNodes
	tmp.findExternalNodesAndFaces();
	externalNodes.add(tmp.getExternalNodes());
	externalFaces.add(tmp.getExternalFaces());
	
	//the tmp elset gets automatically erased
	_externalNodesFound = true;
}


void poorbody :: computeMass()
{
    if ( !isActive() ) return;
    
    feliks_for_each(nodes.begin(),    nodes.end(),    mem_fun(&node::resetMass) );
    feliks_for_each(elements.begin(), elements.end(), mem_fun(&element::lumpMassToNodes) );
	
	mass = 0.0;
    BOOST_FOREACH(node* nd, nodes) mass += nd->getMass();
    _massComputed = true;
}



size_t  poorbody :: dimension() const
{
    return _dimension;
}


std::set<element*>& poorbody :: getElementsTouchingNode(node& nd)
{
    return elementsTouchingNode[&nd];
}


const faceset& poorbody :: getExternalFaces() const
{	
	return externalFaces;
}



double poorbody :: getMass() const
{
	if (!_massComputed) const_cast<poorbody*>(this)->computeMass();
	return mass;
}





feliks::meshing::solidmesh& poorbody ::  getTheSolidmesh()
{
    throw std::runtime_error("This is a fake function for compatibility it should never be called.");
    
    feliks::meshing::solidmesh* m(0);
    return *m;
}




const feliks::meshing::solidmesh& poorbody ::  getTheSolidmesh() const
{
    throw std::runtime_error("This is a fake function for compatibility it should never be called.");
    
    feliks::meshing::solidmesh* m(0);
    return *m;
}




std::set<node*> poorbody :: getNodesFrom0Cells(const std::set<feliks::topology::cell*>& cells) const
{
    throw std::runtime_error("This is a fake function for compatibility it should never be called.");
    
    std::set<node*> m;
    return m;
}




void poorbody :: incrementSolution(const integrator& i)
{
	if (!active) return;
	
	if ( i.isNodallyBased() )
		i.incrementSolution( nodes );
	else
		i.incrementSolution( elements );
}





void poorbody :: initialize()
{
    if (isInitialized()) return;
    
    BOOST_FOREACH(element* el, elements)
    {
        el->setParentSubmesh(*this);
        el->initialize();
    }
    DebugMessage("Elements initialized");
    
	
    // Initialize nodes
	feliks_for_each( nodes.begin(), nodes.end(), mem_fun(&node::initializeDofs) );
	computeMass();
	DebugMessage("Nodes initialized");
    
    // mark boundary nodes
    computeExternalNodesAndFaces();
    feliks_for_each( externalNodes.begin(), externalNodes.end(), mem_fun(&node::setOnBoundary) );
    
    
    // create elementsInNode
    BOOST_FOREACH(element* el, elements)
    {
        for (size_t a=0; a< el->getNNodes(); a++)
        {
            node& nd = el->getNode(a);
            elementsTouchingNode[&nd].insert(el);
        }
    }
    
    
	setInitialized();
	computeMass();
}


#ifdef WITHTBB
// data structure for tbb threaded version
class tbbMaxElementEigenvalue{
public:
	double norm;
	
	tbbMaxElementEigenvalue( vector<element*>& elements) : my_elements(&elements), maxev(0.0), inElement(0){}
	tbbMaxElementEigenvalue( tbbMaxElementEigenvalue& nn, tbb::split) : my_elements(nn.my_elements), maxev(0.0), inElement(0){}
	void join(const tbbMaxElementEigenvalue& y)
    {
        if (maxev < y.maxev)
        {
            maxev = y.maxev;
            inElement = y.inElement;
        }
    }
	
	void operator()( const tbb::blocked_range<size_t>& r)
	{
		vector<element*> &nn( *my_elements );
		for (size_t a=r.begin(); a != r.end(); a++)
		{
			double elev = nn[a]->maxEigenvalue();
			if (elev > maxev)
			{
				maxev = elev;
				inElement = nn[a]->getLabel();
			}
		}
	}
	
public:
	vector<element*> *const my_elements;
	double	maxev;
	int		inElement;
};
#endif




/* finds the maximum eigenvalue for all the elements of the model. In fact, not quite
 the element eigenvalues as much as an upper bound.
 If an element returns 0 as maximum eigenvalue, it means that it does not compute any
 eigenvalue at all
 */
double poorbody :: maxEigenvalueEstimate() const
{
    double  maxev(0.0);
	int maxLabel = 0;
	poorbody* bb = const_cast<poorbody*>(this);

#ifdef WITHTBB
    
	tbbMaxElementEigenvalue tbbmax( bb->elements );
	tbb::parallel_reduce(tbb::blocked_range<size_t>(0, elements.size()),
						 tbbmax,
						 tbb::auto_partitioner() );
	maxev    = tbbmax.maxev;
	maxLabel = tbbmax.inElement;
    
#else
	
	double eigv(0.0);
    BOOST_FOREACH(element* el, bb->elements)
	{
		eigv = el->maxEigenvalue();
		if (maxev < eigv)
		{
			maxev    = eigv;
			maxLabel = el->getLabel();
		}
	}
    
#endif
	
	if (maxev > 1e10)
	{
        std::cout << "enormous eigenvalue in element ";
		exit(0);
	}
    return maxev;
}




// preprocess geometric information of the finite element model of the body
// several tasks are included here
void poorbody :: parse()
{
	if ( _isParsed ) return;
	
	// Task 1: add, to each node, the elements that connect to it
	std::vector<element*>::iterator iter = elements.begin();
    /*
     while (iter != elements.end())
     {
     for (int a=0; a< (*iter)->getNNodes(); a++)
     (*iter)->getNode(a).addElementToConnectivity(**iter);
     ++iter;
     }
     */
    
	// task 3: we count the number of elements containing a node
	vector<node*>::iterator niter = nodes.begin();
	while (niter != nodes.end())
	{
		//(*niter)->elementCounterReset();
		++niter;
	}
	
	iter = elements.begin();
	while (iter != elements.end())
	{
		for (size_t a=0; a< (*iter)->getNNodes(); a++)
		{
            //	(*iter)->getNode(a).incrementElementCount();
		}
		++iter;
	}	
	
	
	colorDisjointElements();
    
    tmpNodes.clear();    
	
	_isParsed = true;
	logger::mainlog << "\n[ Poorbody " << getName() << " parsed]";
}



void poorbody :: parseInformationReset()
{
}



void poorbody :: print(std::ostream& os)
{
    modelpart::print(os);
    os  << "\n Geometry           : " << _geometricDescription
    << "\n Num external nodes : " << externalNodes.size()
    << "\n Exterior nodes name: " << externalNodes.getName();
}





#ifdef WITHTBB
// data structure for tbb threaded version
class tbbResetNodalForces{	
	
public:
	tbbResetNodalForces( vector<node*>& nds) : my_nodes(&nds){}
	
	void operator()( const tbb::blocked_range<size_t>& r) const
	{
		vector<node*> &nn( *my_nodes );
		for (size_t a=r.begin(); a != r.end(); a++)	nn[a]->forceReset();
	}
	
private:
	vector<node*> *const my_nodes;		
};
#endif


// set all the nodal forces to zero
void  poorbody :: resetNodalMasses()
{
	for_each( nodes.begin(), nodes.end(), mem_fun(&node::resetMass));	
    
}

// set all the nodal forces to zero
void  poorbody :: resetNodalForces()
{
	
#ifdef WITHTBB	
	tbb::parallel_for(tbb::blocked_range<size_t>(0,nodes.size()), tbbResetNodalForces(nodes), tbb::auto_partitioner() );	
#else
	for_each( nodes.begin(), nodes.end(), mem_fun(&node::forceReset));	
#endif
}


#ifdef WITHTBB
// data structure for tbb threaded version
class tbbResetSolutionIncrements{	
	
public:
	tbbResetSolutionIncrements( vector<node*>& nds) : my_nodes(&nds){}
	
	void operator()( const tbb::blocked_range<size_t>& r) const
	{
		vector<node*> &nn( *my_nodes );
		for (size_t a=r.begin(); a != r.end(); a++)	nn[a]->setDeltaDofsToZero();
	}
	
private:
	vector<node*> *const my_nodes;		
};
#endif



// set all the nodal solution increments dU to zero
void  poorbody :: resetSolutionIncrements()
{
#ifdef WITHTBB	
	tbb::parallel_for(tbb::blocked_range<size_t>(0,nodes.size()), tbbResetSolutionIncrements(nodes), tbb::auto_partitioner() );	
#else
	for_each( nodes.begin() , nodes.end(), mem_fun(&node::setDeltaDofsToZero) );
#endif
}




void poorbody :: rewindSolution()
{
	// Recover nodal dofs from tn
	for_each( nodes.begin() , nodes.end(), mem_fun(&node::recoverFromBackupDofs) );
}



// this function scans a body from the commandline. Only regularbodies are defined in
// this way, but we put it here so that in the future maybe more elaborated things
// are added
poorbody& poorbody :: scan(const commandLine &cl, model &m)
{
	poorbody *bd(0);
	return *bd;
}





void poorbody :: scanElements(const commandLine &cl, ifstream &modelfile)
{
}



// when scanning nodes, put then in a temporary map and only create the nodes when the
// elements are created, because the type of node is not known until it is associated
// with an element
void poorbody :: scanNodes(const commandLine &cl, ifstream &modelfile)
{    
    ivector coor;
	int     label;
	
    // scan vertex coordinates until blank line or EOF
    string oneLine;
    while( getline(modelfile, oneLine) && oneLine > "                                     ")
    {
		coor.setZero();
		stringstream str(oneLine);
		if (DEBUG_SCANNER) cout << "|" << str.str() << "|" << "\n";
		str >> label >> coor[0] >> coor[1] >> coor[2];
		
		tmpNodes[label] = std::make_pair<ivector,node*>(coor,0);
	}
}




void poorbody :: updateCurrentState(const dofset::evaluation_time when)
{
    feliks_for_each(elements.begin(), elements.end(),
                  std::bind2nd(std::mem_fun(&element::updateCurrentState), when));
    
    feliks_for_each(evalspots.begin(), evalspots.end(),
                  std::bind2nd(std::mem_fun(&evalspot::updateCurrentState), when));   
}


