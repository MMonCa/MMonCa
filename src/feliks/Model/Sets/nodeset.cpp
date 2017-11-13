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
 *  nodeset.cpp
 *  ris
 *
 *  Created by Ignacio Romero on Thu Jul 24 2003.
 *  converted to C++ jan 2006
 *
 */

#include <iostream>
#include <iomanip>
#include <iterator>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <boost/foreach.hpp>

#include "Model/Sets/nodeset.h"
#include "Model/model.h"

#include "Math/tensor.h"
#include "Io/logger.h"
#include "Io/message.h"
#include "Io/usercommand.h"

static int cleanLineSize(const std::string& line);

using namespace std;
using namespace blue;



nodeset :: nodeset() :
name("unnamed_nodeset"),
internal(true),
maxVelocity(0.0,0.0,0.0),
maxVelocityTime(-1.0)
{
}


nodeset :: nodeset(const string& name_) :
name(name_), 
internal(false)
{
}



nodeset :: nodeset(const vector<int>& nodelabels, model &m) :
    internal(false), 
    name("unnamed_nodeset")
{
	vector<int>::const_iterator iter=nodelabels.begin();	
	while (iter != nodelabels.end())
	{
		node &nd(m.getNode(*iter));
		this->insert(&nd);
		++iter;
	}
}


// copies the content of one nodeset into another. Note that the variable 'internal'
// should not be copied, because it signals the way memory for the nodeset was originally
// allocated, and affects how the nodeset is deleted
nodeset :: nodeset(const nodeset& original)
:	
    name(original.name)
{
    BOOST_FOREACH(node* n, original)
        this->insert(n);
}


// copies the content of one nodeset into another. Note that the variable 'internal'
// should not be copied, because it signals the way memory for the nodeset was originally
// allocated, and affects how the nodeset is deleted
nodeset&  nodeset :: operator=(const nodeset &rhs)
{
	name = rhs.name;
    
    BOOST_FOREACH(node* n, rhs)
        this->insert(n);
        
	return *this;
}


nodeset:: ~nodeset()
{
	this->clear();	
}



void nodeset ::  accumulateInertialForces(const ivector& acceleration)
{
	set<node*>::iterator iter = this->begin();
	while ( iter != this->end() )
	{
		ivector force( (*iter)->getMass() * acceleration );
		(*iter)->assembleForce( force );
		++iter;
	}
}



void nodeset :: add(const node &nd)
{
	this->insert(const_cast<node*>(&nd));
}




void nodeset :: add(const nodeset &nds)
{
	set<node*>::const_iterator iter = nds.begin();
	while ( iter != nds.end() )
	{
		this->insert( *iter );
		++iter;
	}
}


// actions to be taken when a solution is obtained which can be considered
// to be converged at the current time step, and the time counter is
// advanced
void nodeset :: advanceInTime()
{
}





// we call this method to obtain the maximum velocity of all the nodes of the nodeset.
// however, this is often called at the beginning of a time step, after the velocity
// has been set to zero as a predictor.
// To bypass this, we compute the last maximum velocity, the max velocity at time tn
const ivector nodeset :: getMaxVelocity() const
{
	extern double global_tn1;
	
	if (maxVelocityTime > 0.0 && (global_tn1 - maxVelocityTime) < 1e-3*maxVelocityTime)
		return maxVelocity;
	
	ivector vmax;
	double  nvmax(0.0);
	
	// search for the maximum velocity in the main node list
	set<node*>::const_iterator iter = this->begin();
	while( iter != this->end() )
	{
		// there is a problem with getVelocity() in explicit methods, where the predictor becomes zero
		ivector v( (*iter)->velocity(dofset::tn) );
		double vnorm = v.norm();
		
		if (vnorm > nvmax)
		{
			nvmax = vnorm;
			vmax  = v;
		}
		++iter;
	}
	
	
	// this function is supposed to be const, but actually it is not because
	// we store the velocity in case it needs to be reused
	const_cast<nodeset*>(this)->maxVelocityTime = global_tn1;
	const_cast<nodeset*>(this)->maxVelocity     = vmax;
	
	return vmax;
}


// we call this method to obtain the Mean  velocity of all the nodes of the nodeset.
// however, this is often called at the beginning of a time step, after the velocity
// has been set to zero as a predictor.
// To bypass this, we compute the last maximum velocity, the max velocity at time tn
const ivector nodeset :: getMeanVelocity() const
{
	extern double global_tn1;
	
	if (meanVelocityTime > 0.0 && (global_tn1 - meanVelocityTime) < 1e-3*meanVelocityTime)
		return meanVelocity;
	
	ivector vmean(0.,0.,0.) ;
	int     nv(0) ; 
	
	// search for the maximum velocity in the main node list
	set<node*>::const_iterator iter = this->begin();
	while( iter != this->end() )
	{
		// there is a problem with getVelocity() in explicit methods, where the predictor becomes zero
		ivector v( (*iter)->velocity(dofset::tn) );
		nv++ ; 
		vmean = vmean + v ; 
		++iter;
	}
	vmean  *= 1.0/ nv ; 
	
	
	// this function is supposed to be const, but actually it is not because
	// we store the velocity in case it needs to be reused
	const_cast<nodeset*>(this)->meanVelocityTime = global_tn1;
	const_cast<nodeset*>(this)->meanVelocity     = vmean ;
	
	return vmean;
}


// we call this method to obtain the Scale  velocity of all the nodes of the nodeset.
// however, this is often called at the beginning of a time step, after the velocity
// has been set to zero as a predictor.
// To bypass this, we compute the last maximum velocity, the max velocity at time tn
/*
 const double nodeset :: getScaleVelocity() const
 {
 extern double global_tn1;
 
 if (scaleVelocityTime > 0.0 && (global_tn1 - scaleVelocityTime) < 1e-3*scaleVelocityTime)
 return scaleVelocity;
 int nv(0.);
 double vscale(0.0) ;
 ivector vmean(0.,0.,0.) ;
 vmean = getMeanVelocity() ; 
 
 // search for the maximum velocity in the main node list
 vector<node*>::const_iterator iter = this->begin();
 while( iter != this->end() )
 {
 // there is a problem with getVelocity() in explicit methods, where the predictor becomes zero
 ivector v( (*iter)->velocity(dofset::tn) );
 nv++ ; 
 vscale += (v-vmean).norm(); 
 ++iter;
 }
 vscale   *= 1.0/ nv ; 
 
 
 // this function is supposed to be const, but actually it is not because
 // we store the velocity in case it needs to be reused
 const_cast<nodeset*>(this)->scaleVelocityTime = global_tn1;
 const_cast<nodeset*>(this)->scaleVelocity     = vscale  ;
 
 return vscale;
 }
 */


bool nodeset :: hasNode(node &nd)
{
    return (this->find(&nd) != this->end());
}



void nodeset :: info(ostream&of) const
{
	of << "\n Nodeset            : " << std::setw(25) << name << " (n nodes: " << size() << ")";
}



void  nodeset :: initializeRotations(nodetypeT ntype, int generate, ivector& dir, ivector& rpoint)
{
	// walk through the array of nodes initialiting the rotation/director
    BOOST_FOREACH(node* nd, (*this))
	{
		switch (ntype)
		{
                // generate directors
			case(_UDnode):
				
				// fixed directors
				if (generate == 0)
				{
					dir.normalize();
					nd->initializeRotationalDofs(dir);
				}
				
				// spherical coordinates, "dir" is the center
				else if (generate == 1)
				{
					ivector &center(dir);
					ivector radius ( nd->getReferenceCoordinates() - center);
					radius.normalize();
					nd->initializeRotationalDofs(radius);
				}
				
				// cylindrical coordintes, "dir" is the direction, "rpoint" is the center
				else if (generate == 2)
				{
					ivector axis(dir), center(rpoint);
					ivector radius ( nd->getReferenceCoordinates() - center);
					ivector axb = radius.cross(axis);
					radius = axis.cross(axb);
					radius.normalize();
					nd->initializeRotationalDofs(radius);
				}
				
				else
					logger::warnings << "Rotation initialization not programmed yet" << "\n";
				
				break;
                
                // generate quaternions
			case (_URnode):
				logger::warnings << "Unable to generate quaternions at the moment" << "\n";
				break;
				
                // generate rotation matrices
			default:
				logger::warnings << "Unable to generate rotation matrices at the moment" << "\n";
				break;
		}
	}
}



void nodeset :: print(ostream &of) const
{
    of  << "\n Nodeset : ";
	
	of.setf(ios::left);
	of	<< setw(35) << name 
    << " [nnode: " << setw(5) << size() << "]" << flush;
	if ( internal ) of << " (internal)" << flush;
	
	if ( Debug_Level() > 0)
	{
		of  << "\t Nodes:" ;
		int c = 0;
        set<node*>::const_iterator iter=this->begin();
        while ( iter != this->end() )
        {
            if (c%6) of << "\n    ";
            of << " " << setw(7) << (*iter)->getLabel() << ",";
            ++iter;
        }
	}
}




void  nodeset :: printDofs(ostream &of) const 
{
	set<node*>::const_iterator iter = this->begin();	
	while (iter != this->end() ) 
	{
		(*iter)->printDOFs(of);
		++iter;
	}
}



// raten = 1: velocity, 2= acceleration, 0 = all
void nodeset :: printRates(int raten, ostream &of) const 
{
	set<node*>::const_iterator iter = this->begin();	
	while (iter != this->end() )  
	{
		(*iter)->printRates(raten, of);
		++iter;
	}
}



void nodeset ::	replaceNode(node& deleted, node& survivor)
{
    this->erase( this->find(&deleted) );
    this->insert( &survivor);
}




// format:
// nodeset, label = ..., name = ...
//   1 , 2, 3,
//   4 , 5, 6, 7:20
//   [empty line]
//

nodeset& nodeset :: scan(const commandLine &cl, ifstream &meshfile, model &m)
{
	// empty nodeset
	nodeset *ndset = new nodeset();
	
	ndset->name  = "";
	int n = cl.size();
	if (n > 1)
	{
		for (int i=1; i<n; i++)
		{
			const usercommand &uc = cl[i];
			if (uc.keyword() == "name")   ndset->name   = uc.option();            
		}
	}		
    
	
	// scan the vertices in the set reading a comma separated list
	string line;
	getline(meshfile, line);
	
	// scan to see if at least there if one node
	int ndlabel, read;
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
			if (read == 2) 	for (int kk=from; kk <= to; kk++) ndset->add( m.getNode(kk) );
		}
		
		else
		{
			read = sscanf(str, "%d", &ndlabel);
			if (read == 1)		ndset->add(m.getNode(ndlabel));
		}
		
		
		// read more nodes if there are any
		while ( (str = strtok(NULL, ", "))  != NULL)
		{
			// read a range of nodes
			dash = strchr(str, ':');
			if ( dash != NULL)
			{
				read = sscanf(str, "%d:%d", &from, &to);
				if (read == 2) for (int kk=from; kk <= to; kk++) ndset->add( m.getNode(kk) );
			}
			
			else
			{
				read = sscanf(str, "%d", &ndlabel);
				if (read == 1)  ndset->add( m.getNode(ndlabel));
			}
		}
	}
	delete []cline;
	
	
	// read if there are more lines
	while ( getline(meshfile,line) && cleanLineSize(line) > 0)
	{
		// scan to see if at least there if one node
		cline = new char[line.size()+1];
		strcpy(cline, line.c_str());
		char *str = strtok(cline, ", ");
		if (str != NULL)
		{
			read = sscanf(str, "%d", &ndlabel);
			if (read == 1) 
			{
				node &nd(m.getNode(ndlabel));
				ndset->add(nd);
			}
			while ( (str = strtok(NULL, ", "))  != NULL)
			{
				read = sscanf(str, "%d", &ndlabel);
				if (read == 1)  
				{
					node &nd(m.getNode(ndlabel));
					ndset->add(nd);
				}
			}
		}
		delete [] cline;
	}
	
	return  *ndset;
}



// checks the size of a line, excluding spaces and tabs
static int cleanLineSize(const std::string& line)
{
	int count(0);
	
	for (int i=0; i<line.size(); i++)
		if ( line[i] != ' ' && line[i] != '\t') count++;
	
	return count;
}

