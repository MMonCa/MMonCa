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
 *  smoother.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 11/11/08.
 *  Copyright 2008 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include "smoother.h"

#include "Geometry/quadpoint.h"
#include "Math/matrix.h"
#include "Model/model.h"
#include "Model/Parts/body.h"
#include "Model/Parts/poorbody.h"
#include "Model/Node/node.h"
#include "Math/tensor.h"
#include "Math/vector.h"
#include "Math/Topology/topology.h"
#include "Io/Postprocessor/resultdata.h"

#include <set>
#include <vector>
#include <iostream>
#include <utility>
#include "boost/foreach.hpp"

#define  CDEBUG 0

using namespace blue;


/*private functions */
/* these are the iterpolation polynomials for several element types, each of them
 computed in a different function  */
static void P_1_x		(const ivector& coor, longvector &p);
static void P_1_x_y		(const ivector& coor, longvector &p);
static void P_1_x_y_xy	(const ivector& coor, longvector &p);
static void P_206		(const ivector& coor, longvector &p);
static void P_303		(const ivector& coor, longvector &p);
static void P_304		(const ivector& coor, longvector &p);
static void P_308		(const ivector& coor, longvector &p);
static void P_310		(const ivector& coor, longvector &p);




void SPRsmoother :: printNodalGradients(const model& m, const std::string& gradname, std::ostream& of)
{
	if ( resultdata::getResultCode(gradname) == RESULT_NONE )
		of << std::endl << " There is no gradient called " << gradname << std::endl;
	
	else
	{
		of	<< std::endl << std::endl << " Projected gradient called " << gradname << std::endl;
        
		resultdata grad(gradname);
        BOOST_FOREACH(body* bd, m.theBodies )
        {
            BOOST_FOREACH( node* nd, bd->nodes)
            {
				this->smoothToNode(grad, *nd, *bd);
				of << endl << setw(8) << nd->getLabel() << "  " << grad;
			}
		}
	}
	of << endl;
}

/*
 // functor to compare the distance of two nodes to a reference one
 class distanceFunctor{
 private:
 const node* ref;
 
 public:
 distanceFunctor(const node& nd) : ref(&nd){}
 bool operator()(const node* pt1, const node* pt2) const
 {
 return squaredNorm( *pt1 - *ref) < squaredNorm( *pt2 - *ref);
 }
 };
 */


static node* nearestNode(const node& nd, std::set<node*> neighbors)
{
	ivector ndpos( nd.coordinates() );
	
	set<node*>::iterator iter=neighbors.begin();
	node* nearest = *iter;
	++iter;
	while (iter != neighbors.end() )
	{
		ivector diter ((*iter)->coordinates() - ndpos);
		ivector dnear (nearest->coordinates() - ndpos);
		
		if ( diter.squaredNorm() < dnear.squaredNorm() ) nearest = *iter;
		++iter;
	}
	return nearest;
}




/* This function computes the "smoothed" stress at a node using the superconvergence
 patch recovery method of Zienkiewicz & Zhu
 
 */
bool SPRsmoother :: smoothToNode(resultdata& r, node& nd, body& md)
{
    
	return true;
}


bool SPRsmoother :: smoothToNode(resultdata& r, node& nd, poorbody& md)
{	
	// quick return for isolated nodes or nodes in invisible elements
    //const size_t dim = md.dimension();
    std::set<element* >& touching = md.getElementsTouchingNode(nd);
    
    
    // quick return for isolated nodes or nodes in invisible elements
	if (touching.empty()) return false;
	
	
	bool cornerNode(false), edgeNode(false);
	r.setZero();
	
	// first, we need to find out the center node for the patch that we are going to
	// use for the smoothing. If the node nd belongs to the interior, its itself. Otherwise
	// we must find a node in the neighborhood of nd that is in the interior, and use the latter.
	// in the neighborhood means that it should be connected by an edge to the former
	node *pnode(0);
	if ( !nd.isOnBoundary() ) 
		pnode      = &nd;
    
	else
	{
		// find the nearest interior node
		std::set<node*> neighbors;
        BOOST_FOREACH(element* el, touching)
		{
			for (int b=0; b< el->getNNodes(); b++)
				if ( !el->getNode(b).isOnBoundary() ) 
					neighbors.insert( &(el->getNode(b)) );
		}
		
		// if the node is on the boundary and there are no interior neighbors, use the node itself, but flag it
		// this happens, for example, in a 1 element thick slab
		if ( neighbors.empty() ) 
		{
			cornerNode = true;
			pnode = &nd;
		}
		
		// if the node is on the boundary, but is it not a corner node, use the closest interior node as center of
		// the spr patch
		else
		{		
			edgeNode = true;
			pnode = nearestNode(nd, neighbors);
		}
	}
	node& patchNode(*pnode);
    std::set<element* >& touchingPatchnode = md.getElementsTouchingNode(patchNode);
    
    
	
	// we need to find out the number of nodes in the SPR patch
	std::set<node*> patchnodes;
    BOOST_FOREACH(element* el, touchingPatchnode)
    {
		for (int b=0; b< el->getNNodes(); b++)
			patchnodes.insert( &(el->getNode(b) ) );
    }
	
    if (CDEBUG)
    {
        cout << endl << "Patch for node " << nd.getLabel() 
        << " has " << patchnodes.size() << " nodes " << flush;
    }
    
	
	// find out the form of the interpolation function P(x,y,z) in the neighbour elements
	// this only works for nodes that connect elements of the same type ... it should be improved
	int nnode = (*(touchingPatchnode.begin()))->getNNodes();
	
	/* select the correct interpolation polynomial */
	const int ndm(3);
	void    (*Pxyz)(const ivector& coor , longvector &p)=NULL;
	switch (100*ndm + nnode)
    {
		case 102: 
		case 302: Pxyz = P_1_x;                      break; //linear element
            //case 103: Pxyz = P_1_x_x2;                   break;
            //case 104: Pxyz = P_1_x_x2_x3;                break;
		case 203: Pxyz = P_1_x_y;                    break;
		case 204: Pxyz = P_1_x_y_xy;                 break;
		case 206: Pxyz = P_206;		break;
		case 303: Pxyz = P_303;		break;
		case 304: Pxyz = P_304;		break;
		case 308: Pxyz = P_308;     break;
		case 310: Pxyz = P_310;     break;
		default: 
			throw runtime_error("\n Error in smoothToNode. Smoothing not programmed for this element ");
			break;
    }
	
	
	// To find the SPR values a small system of equations has to be solved. Find its dimension 
    // and allocate the space for the matrix and RHS
	// the dimension of the polynomial basis is equal to the number of nodes (isoparametric)
	int			npoly = nnode;
	matrix    A(npoly, npoly);
	
	// this an array of projected gradients, one per component of the required gradient
	std::vector<longvector*> intgrad;
	for (int a=0; a < r.getNComponents(); a++) 
		intgrad.push_back(  new longvector(npoly) );
	
	// a longvector that collects the interpolation polynomial base
	longvector p(npoly);
	
    if (CDEBUG)
    {
        cout << endl << "Node " << nd.getLabel() << " is contained in  " 
        <<  touchingPatchnode.size() << " elements " << flush;
        BOOST_FOREACH(element* el, touchingPatchnode)
        cout << " " << el->getLabel() << flush;
    }
    
	// Construct the least-squares system of equations by looping over all the Gauss points
    BOOST_FOREACH(element* el, touchingPatchnode)
    {
		size_t		ngp;
		quadpoint   gauss[10];
		//cout << endl << "Container element has label " << el.getLabel();
		
		
		if ( el->isVisible() )
		{
			// What makes SPR so fast is that the next quadrature rule is reduced. One can select a full quadrature rule
			// but it becomes (approx 8 times) slower (for 3d) without any gain, except in material models with history
			// variables, where the SPR quadpoints do not contain updated history variables, unless taken care of
			
			SPRQuadraturePoints(el->getNNodes(), el->getGeometryType(), cornerNode || edgeNode, &ngp , gauss);
			resultdata  elgrad(r);
			ivector		gpcoor;
			double		dvol;
			
            
			for (int ip=0; ip<ngp; ip++)
			{
				// get from the element the stress at the Gauss point and its coordintes
				el->gradient(ip, elgrad, gpcoor, dvol);
                
				// evaluate the interpolation polynomial
				Pxyz(gpcoor, p);
				
				// evaluate the interpolation matrix (symmetric)
				for (int j=0; j<npoly; j++)
					for (int k=0; k<npoly;  k++) 
						A.data[j][k] += p.data[j] * p.data[k] * dvol;
				
				
				// evaluate the projected gradient
				for (int b=0; b< intgrad.size(); b++)
					for (int j=0; j<npoly; j++)
						intgrad[b]->data[j] += p.data[j] * elgrad[b] * dvol;
			}
		}
	}
	
	// Solve the systems of equations. intgrad[] now holds the best coefficients for the interpolation
	if ( !A.factorLU() )
	{
		cout << "Error in SPR smoother. Singular smoothing operator for node " << nd.getLabel() <<  endl;
	}
	else
		for (int b=0; b < intgrad.size(); b++)
			A.solveLU(*intgrad[b]);
	
	// evaluate the polynomial longvector P at the original node (not at the center of the patch), 
	// multiply by Coeff, to get the recovered grad
	Pxyz(nd.coordinates(), p);
	for (int b=0; b < intgrad.size(); b++)
	{
		r[b] = p.dot( *intgrad[b]);
		delete intgrad[b];
	}
    
	return true;
    
}



/* interpolation polynomial [1 x ] for linear 1d element */
static void P_1_x(const ivector& coor, longvector &p)
{
	p.data[0] = 1.0;
	p.data[1] = coor[0];
}



/* interpolation polynomial [1 x y] for CST */
static void P_1_x_y(const ivector& coor, longvector &p)
{
	p.data[0] = 1.0;
	p.data[1] = coor[0];
	p.data[2] = coor[1];
}



/* interpolation polynomial [1 x y xy] for QUAD */
static void P_1_x_y_xy(const ivector& coor, longvector &p)
{
	p.data[0] = 1.0;
	p.data[1] = coor[0];
	p.data[2] = coor[1];
	p.data[3] = coor[0]*coor[1];
}


/* interpolation polynomial [1 x y z x^2 y^2 xy ] for QTRI */
static void P_206(const ivector& coor, longvector &p)
{
	double x(coor[0]), y(coor[1]);
	
	p.data[0] = 1.0;
	p.data[1] = x;
	p.data[2] = y;
	p.data[3] = x*x;
	p.data[4] = x*y;	
	p.data[5] = y*y;
}




/* interpolation polynomial [1 x y z] for TRIANGULAR shell */
static void P_303(const ivector& coor, longvector &p)
{
	p.data[0] = 1.0;
	p.data[1] = coor[0];
	p.data[2] = coor[1];
}



/* interpolation polynomial [1 x y z] for TET */
static void P_304(const ivector& coor, longvector &p)
{
	p.data[0] = 1.0;
	p.data[1] = coor[0];
	p.data[2] = coor[1];
	p.data[3] = coor[2];
}



/* interpolation polynomial [1 x y z xy yz zx xyz] for BRICK */
static void P_308(const ivector& coor, longvector &p)
{
	double x(coor[0]), y(coor[1]), z(coor[2]);
    
	p.data[0] = 1.0;
	p.data[1] = x;
	p.data[2] = y;
	p.data[3] = z;
	p.data[4] = x*y;
	p.data[5] = y*z;
	p.data[6] = z*x;
	p.data[7] = x*y*z;
}




/* interpolation polynomial [1 x y z x^2 y^2 z^2 xy yz zx] for QTET */
static void P_310(const ivector& coor, longvector &p)
{
	double x(coor[0]), y(coor[1]), z(coor[2]);
	
	p.data[0] = 1.0;
	p.data[1] = x;
	p.data[2] = y;
	p.data[3] = z;
	p.data[4] = x*x;
	p.data[5] = y*y;
	p.data[6] = z*z;
	p.data[7] = x*y;
	p.data[8] = y*z;
	p.data[9] = z*x;
}



bool lumpsmoother :: smoothToNode(resultdata& r, node& nd, body& md)
{	
	bool cornerNode = false;
    bool edgeNode   = false;
	const size_t  dim = md.dimension();
    std::set<feliks::topology::cell* > touching = nd.getCellP()->touchingCells(dim);
    
	// first, we need to find out the center node for the patch that we are going to
	// use for the smoothing. If the node nd belongs to the interior, its itself. Otherwise
	// we must find a node in the neighborhood of nd that is in the interior, and use the latter
	node *pnode(0);
	if ( !nd.isOnBoundary() ) 
	{
		pnode      = &nd;
		cornerNode = false;
	}
	else
	{
        std::set<node*> neighbors;
        BOOST_FOREACH(feliks::topology::cell* c, touching)
		{
			const element* el = md.cell2element[c];
			for (int b=0; b< el->getNNodes(); b++)
				if ( !el->getNode(b).isOnBoundary() ) 
					neighbors.insert( &(el->getNode(b)) );
			
		}
		
		// if the node is on the boundary and there are no interior neighbors, use the node itself, but flag it
		// this happens, for example, in a 1 element thick slab
		if ( neighbors.empty() ) 
		{
			cornerNode = true;
			pnode = &nd;
		}
		
		// if the node is on the boundary, but is it not a corner node, use the closest interior node as center of
		// the spr patch
		else
		{		
			edgeNode = true;
			pnode = nearestNode(nd, neighbors);
		}
	}
	
	
	r.setZero();
	double totalVol(0.0);
	
	// Construct the lumped gradient by looping over all the Gauss points
    BOOST_FOREACH(feliks::topology::cell* c, touching)
    {
        const element& el = *(md.cell2element[c]);
        
        int ngp = numberOfQuadraturePoints(el.getNNodes(), el.getGeometryType());
        resultdata  elgrad(r);
        ivector		gpcoor;
        double		dvol;
		
		for (int ip=0; ip<ngp; ip++)
		{
			// get from the element the stress at the Gauss point and its coordintes
			el.gradient(ip, elgrad, gpcoor, dvol);
			
			r += elgrad*dvol;
			totalVol += dvol;
		}
    }
	
	r *= (1.0/totalVol);
    
	return true;
}


bool lumpsmoother :: smoothToNode(resultdata& r, node& nd, poorbody& md)
{
    return false;
}

