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
 *  elset.h
 *  feliks
 *
 *  Created by Ignacio Romero on Thu Sep 2 2003.
 *
 */

#ifndef _elset_h
#define _elset_h

#include <string>
#include <cstdio>
#include <set>

#include "Model/Sets/faceset.h"
#include "Elements/element.h"
#include "Model/Sets/nodeset.h"



class model;
class node;


class elset : public set<element*>
{
	
public:
                            elset(model &m);
                            elset(model &m, vector<element*> &elmts);
                            elset(const string& name, model &m);
                            elset(const elset& c);
    virtual                 ~elset();

	static elset*           scan(commandLine &cl, ifstream &meshfile, model &m);
	
	void                    findExternalFaces();
	void                    findExternalNodesAndFaces();
	void                    replaceNode(node& deleted, node& survivor);
	
	
	void					add(const element &e);
	void					add(elset &es);
	void					advanceInTime();
	bool					check();	
	faceset&				getExternalFaces() {return externalFaces;};
	const nodeset&			getExternalNodes() const;
	string&                 name();
    const string&           name() const;
	
	
	// print information
	void                    print(std::ostream &os=std::cout) const;
	void					printGradients(ostream& of=cout);
	
	
private:
    std::string             theName;
	model&                  linkedModel;	
	
	// a table of element faces in the exterior of the elset
	bool                    externalFacesFound;
	faceset                 externalFaces;
	
	// a table of nodes in the exterior of the elset
	bool                    externalNodesFound;
	nodeset                 externalNodes;
	
                            elset();
	friend class            model;
};

#endif

