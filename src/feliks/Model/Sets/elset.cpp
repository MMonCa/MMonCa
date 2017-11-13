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
 *  elset.c
 *  feliks
 *
 *  Created by Ignacio Romero on Sept 2 2003.
 *  moved to C++ in feb 2006
 *
 *  An elset is a collection of element pointers (plus some extra data structure)
 *  Note that the elements must be defined and stored elsewhere because erasing
 *  the elset does not eliminate the elements themselves (similarly nodesets)
 */


#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <ostream>
#include <iomanip>
#include <vector>
#include <boost/foreach.hpp>

#include "Model/Sets/elset.h"
#include "Elements/element.h"
#include "Io/logger.h"
#include "Model/model.h"
#include "Io/message.h"
#include "Io/usercommand.h"


using namespace std;




elset :: elset(model &m) 
:	theName("unnamed elset"),
    externalFacesFound(false), externalNodesFound(false),
	linkedModel(m)
{
}


elset :: elset(model &m, vector<element*> &elmts_)
:	theName("unnamed elset"), 
    externalFacesFound(false), externalNodesFound(false), 
	linkedModel(m)
{
    BOOST_FOREACH( element* e, elmts_)
    {
        this->insert(e);
    }
}




elset :: elset(const string& name_, model &m) :
	theName(name_),
    externalFacesFound(false), externalNodesFound(false),
	linkedModel(m)
{
}



elset :: elset(const elset& c)
:	theName(c.theName), 
	linkedModel(c.linkedModel),
	externalFacesFound(c.externalFacesFound),
	externalNodes(c.externalNodes), externalNodesFound(c.externalNodesFound)
{
	set<element*>::const_iterator iter = c.begin();
	while( iter != c.end() )
	{
		this->insert(*iter);
		++iter;
	}
}




elset :: ~elset()
{
    this->clear();
}



void elset :: add(const element &e)
{
	this->insert(const_cast<element*>(&e));
	externalFacesFound = false;
}


void elset :: add(elset &es)
{
	set<element*>::iterator iter = es.begin();
	while ( iter != es.end() )
	{
		this->insert( *iter );
		++iter;
	}	
	externalFacesFound = false;
}



// this function should be called once a solution is given as good
// it allows to update certain internal variables of the elset that
// depend on the solution, for example the bounding box. Notice
// that these variables remain constant even if the solution changes.
// They are only modified when the solution is taken as good
void elset :: advanceInTime()
{

}



bool elset :: check()
{
	bool ok(true);

	if (!ok) logger::warnings << "Error checking elset " << theName << endl;
	return true;
}




/*  this function parses the elset and finds the external faces of its elements.
*	This is done by looking at the faces that are not shared by 2 elmts of the set.
*
*  The goal of this function is to fill up the array es->extElmtface. This array, which is
*	if size nExtFaces, must be first allocated. Each one of its elements is an (pointer)
*  elmtface data type, which contains all the information about the face. See elmtface.h
*/
void elset :: findExternalFaces()
{
    // quick return
    if (externalFacesFound) return;
	
    // first e create a faceset with all faces
    set<element*>::iterator iter = this->begin();
    while (iter != this->end() )
    {
	    for (int i=0; i< (*iter)->nFaces(); i++)
			externalFaces.add( elmtface(**iter,i) );
        ++iter;
    }
	
    //  We sort the tables in terms of the faces
	externalFaces.sort();
	
//	cout << "Elset of size " << size() << endl;
//	cout << "externalFaces sorted, but not unique(" << externalFaces.size() << ")";
//	for_each(externalFaces.begin(), externalFaces.end(), mem_fun_ref(&elmtface::print));
//	cout << endl;
	
	externalFaces.eliminateDuplicatedFaces();
		
//	cout << endl <<"externalFaces sorted, and unique (" << externalFaces.size() << ")";
//	for_each(externalFaces.begin(), externalFaces.end(), mem_fun_ref(&elmtface::print)); 
//	cout << endl;
		
    externalFacesFound = true;
}




void elset :: findExternalNodesAndFaces()
{
	if (externalNodesFound) return;
	
    // the external faces are needed to compute the external nodes
    if (!externalFacesFound) findExternalFaces();
	if (externalFaces.empty()) return;

	
	// get all the nodes in the external faces
	vector<elmtface>::iterator iter=externalFaces.begin();
	while (iter != externalFaces.end())
	{		
		elmtface &ef(*iter);
		for (int i=0; i< ef.getNNodes(); i++)
		{
			const node  &nd(ef.getNode(i));
			externalNodes.add(nd);
		}
		++iter;
	}
	
	
	// set flag
    externalNodesFound = true;
}



const nodeset& elset :: getExternalNodes() const
{
	elset *es=const_cast<elset*>(this);
	es->findExternalNodesAndFaces();
	return externalNodes;
}


string& elset::name()
{
    return theName;
}



const string& elset::name() const
{
    return theName;
}



void  elset::print(std::ostream &os) const
{
    os  << "\n\n Elset                 : " << this->theName
		<< "\n\t total number of elmts : " << this->size() << "\n";
}





void elset :: printGradients(ostream& of)
{
    set<element*>::iterator iter = this->begin();
    while (iter != this->end() )
	{
		(*iter)->printGradients(of);
		++iter;
	}
}


void elset :: replaceNode(node& deleted, node& survivor)
{
	externalNodes.replaceNode(deleted, survivor);
}






// format:
// elset, label = ..., name = ...
//   1 , 2, 3,
//   4 , 5, 6, 7:20
//   [empty line]
//

elset* elset :: scan(commandLine &cl, ifstream &meshfile, model &m)
{
	// empty elset
	elset *eset = new elset(m);
	
	// scan the set label and name
	eset->theName  = "";
	int n = cl.size();
	if (n > 1)
	{
		for (int i=1; i<n; i++)
		{
			const usercommand &uc = cl[i];
			if (uc.keyword() == "name")   eset->theName  = uc.option();
		}
	}		
	
	// scan the vertices in the set reading a comma separated list
	string line;
	getline(meshfile, line);
	
	// scan to see if at least there if one element
	int elabel, read;
	char *cline;
	cline = new char[line.size()+1];
	strcpy(cline, line.c_str());
	char *str = strtok(cline, ", ");
	char *dash(NULL);
	int from, to;

	if (str != NULL)
	{
		// read a range
		dash = strchr(str, ':');
		if ( dash != NULL)
		{
			read = sscanf(str, "%d:%d", &from, &to);
			if (read == 2) for (int kk=from; kk <= to; kk++) eset->add( m.getElement(kk) );
		}
		
		else
		{
			read = sscanf(str, "%d", &elabel);
			if (read == 1)		eset->add( m.getElement(elabel));
		}

		// read more nodes if there are any
		while ( (str = strtok(NULL, ", "))  != NULL)
		{
			// read a range
			dash = strchr(str, ':');
			if ( dash != NULL)
			{
				read = sscanf(str, "%d:%d", &from, &to);
				if (read == 2) for (int kk=from; kk <= to; kk++) eset->add( m.getElement(kk) );
			}
			
			else
			{
				read = sscanf(str, "%d", &elabel);
				if (read == 1)  eset->add( m.getElement(elabel));
			}
		}
	}
	delete []cline;
			
	
	
	// read if there are more lines
	while ( getline(meshfile,line) && line.size() > 0)
	{
		// scan to see if at least there if one node
		cline = new char[line.size()+1];
		strcpy(cline, line.c_str());
		char *str = strtok(cline, ", ");
		if (str != NULL)
		{
			read = sscanf(str, "%d", &elabel);
			if (read == 1) 
			{
				element &e(m.getElement(elabel));
				eset->add(e);
			}
			while ( (str = strtok(NULL, ", "))  != NULL)
			{
				read = sscanf(str, "%d", &elabel);
				if (read == 1)  
				{
					element &e(m.getElement(elabel));
					eset->add(e);
				}
			}
		}
		delete [] cline;
	}
	
	return  eset;
}





