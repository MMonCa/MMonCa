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
 *  modelpart.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 4/11/08.
 *  Copyright 2008 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "Elements/evalspot.h"

#include "Analysis/Dofs/dofset.h"
#include "General/feliksutil.h"
#include "Io/logger.h"
#include "Io/message.h"
#include "Model/Parts/modelpart.h"
#include "Model/Sets/elset.h"


#include <algorithm>
#include <vector>
#include <numeric>
#include <sstream>
#include "boost/foreach.hpp"

#ifdef WITHTBB
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_for_each.h"
#endif

using namespace blue;



modelpart :: modelpart(const std::string& name, model& m) :
    _name(name),
	linkedModel(m),
	active(true),
    _elementsSorted(false),
    _nodesSorted(false), 
	_scheduledForDeletion(false),
	_isInitialized(false)
{
}




modelpart :: ~modelpart()
{
    linkedModel.remove( _externalNodes );
    _externalNodes.clear();
	
    // remove the nodes the vectors, and later, the vectors themselves.
    BOOST_FOREACH(element*  e, elements)    delete e;
	elements.clear();

    BOOST_FOREACH(node*     n, nodes)       delete n;
    nodes.clear();
    
    BOOST_FOREACH(evalspot* v, evalspots)   delete v;
	evalspots.clear();
}




#ifdef WITHTBB
class accumulatedAcceleration
{

    
public:
    double                          value;
    const std::vector<node*>&       my_nodes;

    
    accumulatedAcceleration(const std::vector<node*>& nodes) :
        value(0.0), my_nodes(nodes) {}
    
    accumulatedAcceleration(accumulatedAcceleration& s, tbb::split ) :
        my_nodes(s.my_nodes), value(0.0)    {}
    
    void operator()(const tbb::blocked_range<size_t>& r)
    {
        const size_t re = r.end();
        for (size_t e = r.begin(); e != re; e++)
            value += my_nodes[e]->acceleration().squaredNorm();
    }
    
    void join( accumulatedAcceleration& rhs )
    {
        value += rhs.value;
    }
};

#endif




double  modelpart :: accelerationNormSquared() const
{
    double acnorm2=0.0;
    
#ifdef WITHTBB
    
    accumulatedAcceleration acc(nodes);
	tbb::parallel_reduce(tbb::blocked_range<size_t>(0, nodes.size()), acc);
    acnorm2 = acc.value;
    
#else
    
    std::vector<node*>::const_iterator       iter  = nodes.begin();
    const std::vector<node*>::const_iterator eiter = nodes.end();
    while ( iter != eiter )
    {
        acnorm2 += (*iter)->acceleration().squaredNorm();
        ++iter;
    }
    
    
#endif
    if (isnan(acnorm2))
        cout << " ERROR:: NaN in Acceleration\n";    
    return acnorm2;
}




void modelpart :: add(const element& e)
{
	elements.push_back(const_cast<element*>(&e));
	linkedModel.updateElementLabels(e.getLabel() );
    if (e.getCellP() != 0)
        cell2element[ e.getCellP() ] = &e;
}




void modelpart :: add(const evalspot& ev)
{
    evalspots.push_back( const_cast<evalspot*>(&ev) );
}




void modelpart :: add(const node& nd)
{
    nodes.push_back(const_cast<node*>(&nd));
    _nodesSorted = false;
    linkedModel.updateNodeLabels(nd.getLabel());
    if (nd.getCellP() != 0)
        cell2node[ nd.getCellP() ] = &nd;
}




class changeSignOfNodalForcesWorker
{
public:
	changeSignOfNodalForcesWorker(){}
	void operator()(node* nd) const { nd->force().changeSign();}
};




void modelpart :: changeSignOfNodalForces()
{
    changeSignOfNodalForcesWorker w;
    feliks_for_each(nodes.begin(), nodes.end(), w);
}




bool modelpart :: check() const
{
    int no = std::count_if( nodes.begin(), nodes.end(), mem_fun( &node::isOrphan) );

	if (no >0)
        logger::mainlog << "\n Orphan nodes detected in part " << getName();
    
    return true;
}




void modelpart :: colorDisjointElements()
{
}




elset& modelpart :: getCompleteElset()
{
	if (_theElements == 0)
		_theElements = new elset(linkedModel, elements);
	
	return *_theElements;
}



// This retrieves an element from a body, given its label. We must be careful
// because the element must be in the body. To know this, we should use the
// function hasElement() before calling getElement()
element& modelpart :: getElement(const int label) const
{		
	// using the sorted vector compare until element has the same label as a dummy element
	vector<element*>::const_iterator iter = elements.begin();
    while (iter != elements.end() )
    {
        if ( (*iter)->getLabel() == label) break;
        ++iter;
    }
    
    
	// it might happen that there is no element with the given label
	if (iter == elements.end() )
		throw runtime_error("The element sought is not in the body provided");
	
	return **iter;
}




nodeset& modelpart :: getExternalNodes()
{
    return _externalNodes;
}




const nodeset& modelpart :: getExternalNodes() const
{
    return _externalNodes;
}




// warning: the body must have the node!!
node& modelpart :: getNode(const int label) const
{
	// we do this cast because, despite what it looks like, getting a node
	// from a body modifies the body itself. This is because the body is reordered
	// but we want to keep that detail hidded from the user, so the function is const
	modelpart *bd = const_cast<modelpart*>(this);
	if (!_nodesSorted) bd->sortNodes();
    
	// using the sorted vector compare until element has the same label as a dummy element
	Unode dummy(label);
	vector<node*>::const_iterator iter;
	iter = lower_bound( nodes.begin(), nodes.end(), &dummy, node::nodepLess );
	
	// it might happen that there is no element with the given label
	if (iter >= nodes.end() )
	{
		stringstream st("Node with label ");
		st << "Node with label " << label << "is not in the modelpart";
		throw runtime_error(st.str());
	}
	return **iter;
}




int modelpart :: getNVisibleElements() const
{	
	return std::count_if( elements.begin(), elements.end(), mem_fun(&element::isVisible) );
}




bool modelpart :: hasNode(const int label) const
{
	modelpart *bd = const_cast<modelpart*>(this);
	bd->sortNodes();
    
	Unode dummy(label);	
	return binary_search ( nodes.begin(), nodes.end(), &dummy, node::nodepLess );
}




bool modelpart :: integrateDampingMatrix(assembler& as) const
{
    bool ret = false;
    if ( active )
    {
        ret = true;
        BOOST_FOREACH( element* ev, elements)
        ev->integrateDampingMatrix(as);			
    }
    return ret;
}




bool modelpart :: integrateDEnergy() const
{
    if ( !active ) return false;
    
    feliks_for_each(evalspots.begin(), evalspots.end(), mem_fun(&evalspot::integrateDEnergy));
    feliks_for_each(elements.begin(), elements.end(), mem_fun(&element::integrateDEnergy));
    
    return true;
}




#ifdef WITHTBB

// data structure for tbb threaded version
class tbbComputeEnergySpots
{
public:
	energies tenergy;
	
	tbbComputeEnergySpots(vector<evalspot*>& ev, const dofset::evaluation_time& when) :  my_spots(ev), evtime(when)
    { 
        tenergy.setZero(); 
    }
    
	tbbComputeEnergySpots( tbbComputeEnergySpots& nn, tbb::split) : my_spots(nn.my_spots), evtime(nn.evtime)
    { 
        tenergy.setZero(); 
    }
    
	void join(const tbbComputeEnergySpots& y) 
    {
        tenergy += y.tenergy;
    }
	
	void operator()( const tbb::blocked_range<size_t>& r)
	{
        const size_t re = r.end();
        
        for (int e = r.begin(); e != re; e++)
        {
			energies elenergy;
			my_spots[e]->integrateEnergy(elenergy, evtime);
			tenergy += elenergy;
		}
	}
	
private:
    std::vector<evalspot*>&         my_spots;
    const dofset::evaluation_time   evtime;
};




class tbbComputeEnergyElements
{
public:
	energies tenergy;
	
	tbbComputeEnergyElements( vector<element*>& ev, const dofset::evaluation_time& when) : 
        my_elements(ev), evtime(when)
    { 
        tenergy.setZero(); 
    }
    
	tbbComputeEnergyElements( tbbComputeEnergyElements& nn, tbb::split) : my_elements(nn.my_elements), evtime(nn.evtime)
    { 
        tenergy.setZero(); 
    }
    
	void join(const tbbComputeEnergyElements& y) 
    {
        tenergy += y.tenergy;
    }
	
	void operator()( const tbb::blocked_range<size_t>& r)
	{
        const size_t re = r.end();

        for (size_t e= r.begin(); e != re; e++)
        {
			energies elenergy;
			my_elements[e]->integrateEnergy(elenergy, evtime);
			tenergy += elenergy;
		}
	}
	
private:
    std::vector<element*>&          my_elements;
    const dofset::evaluation_time   evtime;
};

#endif






bool modelpart :: integrateEnergy(energies& energy, const dofset::evaluation_time& when) const
{
	energy.setZero();
    
#ifdef WITHTBB
    
	modelpart* bb = const_cast<modelpart*>(this);
	tbbComputeEnergySpots tbbenergy(bb->evalspots, when);
	tbb::parallel_reduce(tbb::blocked_range<size_t>(0, evalspots.size()), tbbenergy);
    energy  = tbbenergy.tenergy;
	
    tbbComputeEnergyElements tbbenergyE(bb->elements, when);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, elements.size()), tbbenergyE);    
    energy += tbbenergyE.tenergy;

    
#else

    BOOST_FOREACH(evalspot* ev, evalspots)
    {
        energies elenergy;
        ev->integrateEnergy(elenergy, when);
        energy += elenergy;
    }

	BOOST_FOREACH(element* e, elements)
	{
		energies elenergy;
		e->integrateEnergy(elenergy, when);
		energy += elenergy;
	}
	
#endif
    
	return true;
}



bool modelpart :: integrateMassMatrix(assembler& as) const
{
    bool ret = false;
    if (active)
    {
        ret = true;
        
        BOOST_FOREACH(evalspot* e, evalspots)
                e->integrateMassMatrix(as);
        
        BOOST_FOREACH(element* e, elements)
                e->integrateMassMatrix(as);        
    }
    return ret;
}




void modelpart :: localizeSolutionIncrement(const longvector& delta,  partitionDofs pt)
{
    BOOST_FOREACH(node* n, nodes)
		n->localizeSolutionIncrement(delta, pt);
} 




/* Go through all the elements in the body and mark those that have
 at least a node with a constrained dof. This makes the assembly of
 the residual and tangent faster
 */
void modelpart :: markConstrainedElements()
{
	// Initialize elements and set dofs/node
	for_each( elements.begin(), elements.end(), mem_fun(&element::markConstrained) );
}



size_t modelpart :: maxDofsPerNode() const
{
	size_t m = 0;
    
	BOOST_FOREACH(node* n, nodes) 
        m = std::max<int>( m, n->getNDofs() );
    
	return m;
}




void modelpart :: print(std::ostream&of)
{
    of  << "\n"
        << "\n Part name          : " << getName()
        << "\n Manifold dimension : " << dimension()
        << "\n Center             : ("<< center[0] << ", " << center[1] << ", " << center[2] << ")"
        << "\n Orientation        : ("<< orientation[0] << ", "<<orientation[1]<<", "<<orientation[2] << ")"
        << "\n Number of elements : " << elements.size()
        << "\n Number of nodes    : " << nodes.size()
        << "\n Number of spots    : " << evalspots.size();
}




bool modelpart :: remove(element& e)
{
    std::vector<element*>::iterator i = find( elements.begin(), elements.end(), &e);
    
    if ( i == elements.end() ) return false;
    
	elements.erase(i);
    return true;
}



bool modelpart :: remove(node& n)
{
    std::vector<node*>::iterator i = find( nodes.begin(), nodes.end(), &n);
    
    if ( i == nodes.end() ) return false;
    
	nodes.erase(i);
    return true;
}




void modelpart :: replaceNode(node& deleted, node& survivor)
{
	vector<element*>::iterator iter=elements.begin();
	while (iter != elements.end() )
	{
		(*iter)->replaceNode(deleted, survivor);
		++iter;
	}
}



// modify the coordinates of all the body nodes by rotating them and
// shifting them.
void modelpart :: rotateAndShift()
{
#ifndef WITHEIGEN
	iquaternion qorientation( orientation );	
#else
	Vector3d vv( orientation );
	iquaternion qorientation( AngleAxisd(vv.norm(), vv.normalized() ) );
#endif
	
	irotation rotation(qorientation);
    BOOST_FOREACH(node* n, nodes)
	{
		ivector coord( n->getReferenceCoordinates() );
		
		// rotate and shift reference coordinates
		coord  = rotation*coord;
		coord += center;
		
		// store modified coordinates
		n->setReferenceCoordinates(coord);
	}
	
}



void modelpart :: sortElements()
{
	if (_elementsSorted) return;
	
	// sort by label
	std::sort(elements.begin(), elements.end(), element::epLess);

	// elimiate repeated values
	vector<element*>::iterator end_unique = std::unique(elements.begin(), elements.end());
	elements.erase(end_unique, elements.end());
	
	_elementsSorted = true;
}




void modelpart :: sortNodes()
{
	if (_nodesSorted) return;
	
	// sort by label
	std::sort(nodes.begin(), nodes.end(), node::nodepLess);
	_nodesSorted = true;
}



void  printNodes(const modelpart& b, ostream &of)
{
	// check to avoid printing empty bodies
	if (b.nodes.empty()) return;
	
	// special header for the first node
	of << "\n\n  Body : " << b.getName() << "\n\n";
    of << "\n label             coord1               coord2               coord3";
	
	BOOST_FOREACH( node* nd, b.nodes )
	{
		of << "\n" << setw(6) << noshowpos << nd->getLabel();
		of << setprecision(6) << scientific << right << showpos;
		of << nd->getReferenceCoordinates();
    }
    of << "\n" << noshowpos;
}

