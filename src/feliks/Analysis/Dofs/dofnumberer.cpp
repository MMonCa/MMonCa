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
 *  dofnumberer.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 25/6/09.
 *  Copyright 2009 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include "Analysis/Dofs/dofnumberer.h"

#include "Elements/element.h"
#include "Model/Parts/modelpart.h"
#include "Model/Node/node.h"
#include "idmap.h"

#include <vector>

using namespace std;

DOFnumberer :: DOFnumberer() :
	_nextDOF(0)
{
	
}



void DOFnumberer ::	resetDOFNumbers(const std::list<modelpart*>::iterator start, const std::list<modelpart*>::iterator end)
{
	std::list<modelpart*>::iterator i=start;
	while (i != end)
	{
		resetDOFNumbers(**i);
		++i;
	}
}



void DOFnumberer :: resetDOFNumbers(modelpart& b)
{
	std::vector<node*>::iterator iter=b.nodes.begin();	
	while (iter != b.nodes.end() )
	{	
		node &nd = **iter;
		
		if (nd.getNodetype() == _Unode)
			dynamic_cast<Unode&>(nd).getUDS().getIDs().reset() ;
		
		else if (nd.getNodetype() == _UPnode)
		{
			dynamic_cast<UPnode&>(nd).getUDS().getIDs().reset() ;
			dynamic_cast<UPnode&>(nd).getPDS().getIDs().reset() ;
		}		
		
		
		else if (nd.getNodetype() == _UDnode)
		{
			dynamic_cast<UDnode&>(nd).getUDS().getIDs().reset() ;
			dynamic_cast<UDnode&>(nd).getDDS().getIDs().reset() ;
		}		
		
		else if (nd.getNodetype() == _Pnode)
			dynamic_cast<Pnode&>(nd).getPDS().getIDs().reset() ;

		++iter;
	}
}



void plainDOFnumberer :: assignDOFNumbers(const std::list<modelpart*>::iterator start, 
                                          const std::list<modelpart*>::iterator end)
{
	std::list<modelpart*>::iterator i=start;
	while (i != end)
	{
		assignDOFNumbers(**i);
		++i;
	}
}



void plainDOFnumberer :: assignDOFNumbers(modelpart& b)
{
	std::vector<node*>:: iterator iter=b.nodes.begin();	
	while (iter != b.nodes.end() )
	{	
		node &nd = **iter;
		
		if (nd.getNodetype() == _Unode)
			fillIDs( dynamic_cast<Unode&>(nd).getUDS().getIDs() );
		
		else if (nd.getNodetype() == _UPnode)
		{
			fillIDs( dynamic_cast<UPnode&>(nd).getUDS().getIDs() );
			fillIDs( dynamic_cast<UPnode&>(nd).getPDS().getIDs() );
		}		
				
		else if (nd.getNodetype() == _UDnode)
		{
			fillIDs( dynamic_cast<UDnode&>(nd).getUDS().getIDs() );
			fillIDs( dynamic_cast<UDnode&>(nd).getDDS().getIDs() );
		}		
		
		else if (nd.getNodetype() == _Pnode)
			fillIDs( dynamic_cast<Pnode&>(nd).getPDS().getIDs() );
		
		
		++iter;
	}	
}



void plainDOFnumberer :: assignDOFNumbers(element& b)
{
    for (int a=0; a< b.getNNodes(); a++)
	{	
		node &nd = b.getNode(a);
		
		if (nd.getNodetype() == _Unode)
			fillIDs( dynamic_cast<Unode&>(nd).getUDS().getIDs() );
		
		else if (nd.getNodetype() == _UPnode)
		{
			fillIDs( dynamic_cast<UPnode&>(nd).getUDS().getIDs() );
			fillIDs( dynamic_cast<UPnode&>(nd).getPDS().getIDs() );
		}		
		
		
		else if (nd.getNodetype() == _UDnode)
		{
			fillIDs( dynamic_cast<UDnode&>(nd).getUDS().getIDs() );
			fillIDs( dynamic_cast<UDnode&>(nd).getDDS().getIDs() );
		}		
		
		else if (nd.getNodetype() == _Pnode)
			fillIDs( dynamic_cast<Pnode&>(nd).getPDS().getIDs() );		
	}	
}




void DOFnumberer ::	fillIDs(IDmap& m)
{
	for (int a=0; a< m._ndof; a++)
		if ( m.isUndefined(a) ) 	m._map[a] = _nextDOF++;
}



splitDOFnumberer :: splitDOFnumberer() :
	_UnextDOF(0),
	_PnextDOF(0)
{}



void splitDOFnumberer :: assignDOFNumbers(const std::list<modelpart*>::iterator start, 
                                          const std::list<modelpart*>::iterator end)
{
	std::list<modelpart*>::iterator i=start;
	while (i != end)
	{
		assignDOFNumbers(**i);
		++i;
	}
}


void splitDOFnumberer :: assignDOFNumbers(modelpart& b)
{
	//assignUDOFNumbers(b); // MMM This could be a total breakout under some circ... works properly so far with the split chorin projection method.
	assignPDOFNumbers(b);
}



void splitDOFnumberer :: assignDOFNumbers(element& e)
{
    cout << "Do not use splitDOFnumberer for single element";
}



void splitDOFnumberer :: assignUDOFNumbers(const std::list<modelpart*>::iterator start, 
                                           const std::list<modelpart*>::iterator end)
{
	std::list<modelpart*>::iterator i=start;
	while (i != end)
	{
		assignUDOFNumbers(**i);
		++i;
	}
}



void splitDOFnumberer :: assignUDOFNumbers(modelpart& b)
{
	// loop through all the nodes numbering the U dofs
	std::vector<node*>:: iterator iter=b.nodes.begin();
	setNextDOFtoAssign(_UnextDOF);
	while (iter != b.nodes.end() )
	{	
		node &nd = **iter;
		
		if (nd.getNodetype() == _Unode || nd.getNodetype() == _URnode ||
			nd.getNodetype() == _UDnode || nd.getNodetype() == _Pnode)
			throw runtime_error("Error in splitDOFnumberer. This is only for coupled problems");
		
		else if (nd.getNodetype() == _UPnode)
			fillIDs( dynamic_cast<UPnode&>(nd).getUDS().getIDs() );
		++iter;
	}	
	_UnextDOF = getLastDOFAssigned() + 1;
}



void splitDOFnumberer :: assignPDOFNumbers(const std::list<modelpart*>::iterator start, const std::list<modelpart*>::iterator end)
{
	std::list<modelpart*>::iterator i=start;
	while (i != end)
	{
		assignPDOFNumbers(**i);
		++i;
	}
}



void splitDOFnumberer :: assignPDOFNumbers(modelpart& b)
{
	// loop through all the nodes numbering the P dofs
	std::vector<node*>:: iterator iter=b.nodes.begin();
	_PnextDOF = _UnextDOF; // MMM. Careful on this.
    setNextDOFtoAssign(_PnextDOF);
	while (iter != b.nodes.end() )
	{	
		node &nd = **iter;
		
		if (nd.getNodetype() == _Unode || nd.getNodetype() == _URnode ||
			nd.getNodetype() == _UDnode || nd.getNodetype() == _Pnode)
			throw runtime_error("Error in splitDOFnumberer. This is only for coupled problems");
		
		else if (nd.getNodetype() == _UPnode)
			fillIDs( dynamic_cast<UPnode&>(nd).getPDS().getIDs() );
		
		++iter;
	}	
	_PnextDOF = getLastDOFAssigned() + 1;
}


