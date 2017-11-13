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
 *  faceset.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 31/3/10.
 *  Copyright 2010 Universidad Politécnica de Madrid. All rights reserved.
 *
 *  There is an important difference between nodeset, elsets and facesets. The 
 *  latter is not a collection of pointers, but faces. This is because it seems
 *  reasonable that faces themselves have no need to be reused in other facesets,
 *  nor in any other computation.
 *
 */

#include "faceset.h"

#include <algorithm>
#include <vector>
#include <iomanip>
#include "boost/foreach.hpp"
#include "Model/Sets/nodeset.h"


faceset :: faceset() :
	name("unnamed faceset")
{
	nodesInFaceset.declareInternal();
}


faceset :: faceset(const std::string& n) :
	sorted(true)
{
	setName(n);
	nodesInFaceset.declareInternal();
}


faceset :: ~faceset()
{
 	faces.clear();
} 


void faceset :: add(const elmtface &f)
{
	faces.push_back(f);
	sorted = false;
}


void faceset :: add(const faceset &fs)
{
	faces.insert(faces.end(), fs.begin(), fs.end() );
	sorted = false;
}


void faceset :: eliminateDuplicatedFaces()
{
	// eliminate every entry in the vector that appears twice. Note that a face in a body
	// can not appear more than twice.
	sort();
	//int ntama = faces.size() ; 
	std::vector<elmtface>::iterator laux = faces.begin()+1 ;
	while(laux < faces.end() ) {
		if (*(laux)==*(laux-1)) {
			laux = faces.erase(laux-1,laux+1) ;
			++laux ; 
		}
		else 
			++laux ; 
	}
	// faces.erase( std::unique( faces.begin(), faces.end() ), faces.end() );
	//int ntama2 = faces.size();
}



nodeset& faceset :: getNodes()
{
	const_cast<faceset&>(*this).sort();
	
	return nodesInFaceset;
}


void faceset ::	info(std::ostream& of) const
{
	of << "\n Faceset            : " << std::setw(25) << name << " (n faces: " << size() << ")";
}


void faceset ::	print(std::ostream& of) const
{
	info(of);
	std::vector<elmtface>::const_iterator iter = faces.begin();
	while ( iter != faces.end())
	{
		iter->print(of);
		++iter;
	}
}


void faceset :: setName(const std::string& n)
{
	name = n;
	nodesInFaceset.setName(name + "_nodes");
}


void faceset :: sort()
{
	if (sorted || empty() ) return;
	std::sort(faces.begin(), faces.end());
	
	// create a nodeset with the unique nodes in the faceset
	nodesInFaceset.empty();
	BOOST_FOREACH( elmtface f, faces) nodesInFaceset.add( f.getNodes() );
	
	sorted = true;
}

