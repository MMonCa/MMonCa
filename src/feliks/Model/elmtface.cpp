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
 *  elmtface.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 24/1/05.
 *
 *  the face nodes are stored in an integer array of fixed length for performance issues
 *  This makes the allocation much faster ... but less flexible. 
 */

#include "elmtface.h"

#include "General/idata.h"
#include "Elements/element.h"
#include "assert.h"
#include <vector>

using namespace std;
using namespace blue;


elmtface :: elmtface(const element &e_, const int facenum) :
	e(&e_)
{
	nnodes = e->getNNodesPerFace();
	assert( nnodes <= FELIKS_MAX_FACENODES);
	for (int i=0; i<nnodes; ++i) 
	{
		nodes[i]  = &(e->getNodeInFace(facenum, i));
		labels[i] = nodes[i]->getLabel(); 
	}
}




// copy constructor
elmtface :: elmtface(const elmtface &ef)
{
	e		    = ef.e;
	nnodes		= ef.nnodes;
	for (int i=0; i<nnodes; ++i) 
	{
		nodes[i]  = ef.nodes[i];
		labels[i] = ef.labels[i];
	}
}


elmtface :: elmtface(element *ec,node *n1,node *n2,node *n3,node  *n4) 
{
	e = ec ;
	nnodes = 4 ;
	nodes[0] = n1 ; 
	nodes[1] = n2 ;
	nodes[2] = n3 ;
	nodes[3] = n4 ;
	for (int i=0; i< nnodes ; i++) labels[i] = nodes[i]->getLabel();
}

// this function is used to compare two elmtface, and see if they are identical, i.e., if
// they contain the same nodes AND they belong to different elements
// it is used sometimes used by sorting algorithms
// the only thing that must be tested is if both element faces contain the same node (labels)
// and this we must do without touching the order of the labels, which indicate the node
// order within the element
// if (*this) is smaller than ef2 return -1
// if (*this) is equal to ef2, return -1
// if (*this) is larger than ef2, return +1
//
// the faces should not be modified after the comparison.

int elmtface ::  compare(const elmtface &ef2) const 
{
	int cmp;
	
	switch (getNNodes())
	{
		// make a faster comparison for triangles, that does not require memory allocation/deallocation
		case 3:
			
			// create tmp A and B to hold the labels, and orden them
			int A[3], B[3];
			A[0] = labels[0];
			A[1] = labels[1];
			A[2] = labels[2];
			if ( A[0] > A[1] ) swap(A[0], A[1]);
			if ( A[1] > A[2] ) swap(A[1], A[2]);
			if ( A[0] > A[1] ) swap(A[0], A[1]);
			
			B[0] = ef2.labels[0];
			B[1] = ef2.labels[1];
			B[2] = ef2.labels[2];
			if ( B[0] > B[1] ) swap(B[0], B[1]);
			if ( B[1] > B[2] ) swap(B[1], B[2]);
			if ( B[0] > B[1] ) swap(B[0], B[1]);
			
			for (int i=0; i<3; i++)
			{
				if ( A[i] < B[i] ) return -1;
				if ( A[i] > B[i] ) return +1;
			}
			cmp = 0;
			
			break;
			
			
		default:
			idata id1, id2;
			id1 = NewIdata( getNNodes() );
			id2 = NewIdata( getNNodes() );
			FillIdata(id1, getNNodes(), const_cast<int*>(labels));
			FillIdata(id2, getNNodes(), const_cast<int*>(ef2.labels));
	
			IdataOrder(id1);
			IdataOrder(id2);
			cmp = IdataCompare(id1, id2);
	
			FreeIdata(id1);
			FreeIdata(id2);
	}
	return cmp;
}



ivector elmtface :: getCenter() const
{
	ivector c;
	c.setZero();
	for (int i=0; i<nnodes; ++i) c += nodes[i]->coordinates();
	return c*(1.0/static_cast<double>(nnodes));
}


int elmtface ::	getElementLabel() const 
{
	return e->getLabel();
}



nodeset elmtface :: getNodes() const
{
	nodeset nds;
	for (int i=0; i<nnodes; ++i) nds.add(getNode(i));
	
	return nds;
}



bool operator<(const elmtface &ef1, const elmtface &ef2) 
{
	return ef1.compare(ef2)<0;
}



bool operator==(const elmtface &ef1, const elmtface &ef2)
{
	return ef1.compare(ef2) == 0;
}


void elmtface :: print(std::ostream& of) const
{
	of  << "\n Face of element  " << setw(5) << e->getLabel() << ": [";
	for (int i=0; i<nnodes; ++i) of << setw(5) << labels[i];;
	of << "]";
}
