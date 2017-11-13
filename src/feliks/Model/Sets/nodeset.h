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
 *  nodeset.h
 *  feliks
 *
 *  Created by Ignacio Romero on Thu Jul 24 2003.
 *
 * nodesets are lists of nodes used throughout the code for collecting information
 * or treating a set of nodes in a similar way. The main feature of a nodeset is that
 * it is not the owner of its nodes. 
 *
 * Nodesets are initialized from the
 * model file at the beginning of the execution of feliks, but can also be created
 * afterwards, and even modificated. The main function to add a node to a nodeset
 * is the function 'add', which is the only one that should be used after the creation
 * of the nodeset.
 * 
 * To make the process faster, and avoid resizing of the memory, it is advisable (not
 * mandatory) to reserve the required memory prior to filling up the node data. To that
 * purpose, the function 'reserve' is provided. If during execution, size exceeds the
 * reserved memory, the class will take care of it, but it will consume some time.
 *
 * nodesets might contain nodes belonging to different bodies and hence the global
 * collection of nodesets is stored in the model
 *
 * 
 */

#ifndef _nodeset_h
#define _nodeset_h
	
#include <iostream>

#include <set>
#include <string>

#include "Model/Node/node.h"

using namespace std;

class model;
class commandLine;

class nodeset : public std::set<node*>
{	
	
public:
                            nodeset();
                            nodeset(const string& name);
                            nodeset(const vector<int>& nodelabels, model &m);
                            nodeset(const nodeset& copy);
                            nodeset& operator=(const nodeset &rhs);
    virtual                 ~nodeset();
		
	static nodeset&			scan(const commandLine &cl, ifstream &meshfile, model &m);
	
	void					add(const node &nd);
	void					add(const nodeset &nds);
	void					advanceInTime();
	
	// get information from the nodeset
	const blue::ivector			getMaxVelocity() const;
	const blue::ivector           getMeanVelocity() const ; 
	const string			getName() const;
	bool					hasNode(node &nd);
	bool					isInternal() const;
	
	
	// operations with nodesets
	void					accumulateInertialForces(const blue::ivector& acceleration);
	inline void				declareInternal();
	void					initializeRotations(nodetypeT ntype, int generate, blue::ivector& dir, blue::ivector& rpoint);
	void					replaceNode(node& deleted, node& survivor);
	void					setName(const std::string& name_);
	
	// outputs
	void					info(ostream&of=cout) const;
	void					print(ostream &of=cout) const;
	void					printDofs(ostream &of=cout) const;
	void					printRates(int raten, ostream &of=cout) const;

    
private:
    string                  name;
	bool                    internal;
	blue::ivector                 maxVelocity;
	double                  maxVelocityTime;	
	blue::ivector                 meanVelocity;
	double                  meanVelocityTime;	
	
    
	friend class            model;
};



// inlines
inline void				nodeset ::	declareInternal()		{internal=true;}
inline const string		nodeset ::	getName()	const		{return name;}
inline bool				nodeset ::  isInternal() const		{return internal;}
inline void				nodeset	::	setName(const std::string& name_) {name = name_;}


#endif

