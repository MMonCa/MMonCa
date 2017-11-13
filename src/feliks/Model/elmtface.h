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
 *  elmtface.h
 *  MacFELIKS
 *
 *  Created by Ignacio Romero on 24/1/05.
 *
  a face is given by a (ordered) list of the labels of the nodes that belong to it.
 */

#ifndef _elmtface_h
#define _elmtface_h

#include <vector>
#include <iostream>

#include "Main/feliks.h"
#include "Math/tensor.h"
#include "Model/Sets/nodeset.h"

class element;
class node;



class elmtface {
	
	
public:

	elmtface(const element &el, const int facenum_);
	elmtface(const elmtface &ef);		// copy constructor
	elmtface(element *el,node *n1,node *n2,node *n3,node  *n4);

	
	// retrieve information
	blue::ivector			getCenter() const;
	int						getElementLabel() const;
	inline const element&	getElement(){return *e;}
	int						getNNodes() const {return nnodes;}
	
	// retrieve individual nodes or all nodes
	const node&				getNode(int local) const;
	nodeset					getNodes() const;

	
	// display information
	void					print(std::ostream& of=std::cout) const;
	
	
	
	// this functions defines an order relation by comparing labels
	friend bool operator<(const elmtface &ef1, const elmtface &ef2);
	friend bool operator==(const elmtface &ef1, const elmtface &ef2);


private:	
    int			nnodes;
	int			labels[FELIKS_MAX_FACENODES];
	const node* nodes[FELIKS_MAX_FACENODES];
	const element	  *e;					// a pointer to the element to which the face belongs
	
	int   compare(const elmtface &ef2) const;
	friend    class elset;
};

inline const node& elmtface :: getNode(int local) const{return *(nodes[local]);}




#endif
