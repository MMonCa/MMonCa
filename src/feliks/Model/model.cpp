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
 * model.cpp
 *
 * iro, november 1999. revised july 2003
 *
 note for further development. When loading the model and the model,
 the first thing to store is the element types. This set the
 number of dofs per vertex maximum... then goes the vertices, elements...
 *
 *
 * converted to C++ in jan 2006
 */


#include "Model/model.h"

#include <algorithm>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <list>
#include <numeric>
#include <string>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <boost/foreach.hpp>

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "Analysis/Dofs/dofset.h"
#include "Analysis/Dofs/dofnumberer.h"
#include "Elements/element.h"
#include "Elements/eltype.h"
#include "General/feliksutil.h"
#include "Main/feliks.h"
#include "Main/feliksinterface.h"

#include "Model/constraint.h"
#include "Model/initcondition.h"
#include "Model/initrate.h"

#include "Model/Parts/poorbody.h"
#include "Model/Parts/body.h"
#include "Model/Sets/elset.h"
#include "Model/Sets/nodeset.h"
#include "Model/Node/node.h"
#include "Model/Node/emptynode.h"
#include "Analysis/Loading/pointload.h"

#include "General/idata.h"
#include "Io/logger.h"
#include "Io/message.h"
#include "Materials/material.h"
#include "Io/usercommand.h"

#ifdef WITHMPI
#include "mpi.h"
#endif


#define DEBUG_SCANNER 0

using namespace std;
using namespace blue;



model :: model() :
dofsPerNode(0), nDofs(0), theNumberer(0), mass(0.0),
smallest_elmt(0), smallest_node(0),
_mustPrint(false), _isParsed(false), _isInitialized(false),
modelfilename("")
{
	// create a default body in the model. Every orphan element and node goes here
	poorbody* defb= new poorbody("Default", *this);
	add(*defb);
	
	// set the default DOFnumberer
	theNumberer = new plainDOFnumberer();
}




// syntax is as follows:
void model :: fill(const feliks::mmoncaInterface& mi)
{		
	DebugMessage("Filling the model");
	
    // creating nodes
    {
        size_t nn = mi.vcoordinates_.size()/3;
        
        for (size_t a=0; a<nn; a++)
        {
            ivector coor( mi.vcoordinates_[3*a], mi.vcoordinates_[3*a+1], mi.vcoordinates_[3*a+2] );
            node* nd = new Unode (a, coor);
        
            thePoorbodies.front()->add(*nd);
        }
    }
    
    
    // creating elements
    {
        size_t ne = mi.connectivity_.size()/8;
        node** vlist = new node*[8];
        eltype*  et  = eltype::allEltypes[0];
        
        for (size_t a=0; a<ne; a++)
        {
            for (size_t b=0; b<8; b++)
                vlist[b] = thePoorbodies.front()->nodes[ mi.connectivity_[8*a+b] ];
            
            element *e = et->createElement(a, 0, 8, vlist, mi.theElements_[a]);
            thePoorbodies.front()->add(*e);
        }
        
        delete [] vlist;
    }
}




bool    model:: isMEinterpolation(const commandLine& cline)
{
    for (int k=1; k< cline.size(); k++)
    {
        const usercommand  &tmpuc = cline[k];
        if  (tmpuc.keyword() == "interpolation" && tmpuc.option() == "me")
            return true;
    }
    return false;
}



// the destructor goes through all the data stored inside the model and
// deletes it. There is no need to go later and search for all the nodes,
// all the elements, etc. and delete them.
model :: ~model()
{		
	// free the constraints
	while ( !spconstraints.empty() )
	{
		delete spconstraints.front();
		spconstraints.pop_front();
	}	
	
	// free the pointloads
	while ( !theLoads.empty() ) 
	{
		delete theLoads.front();
		theLoads.pop_front();
	}
	
	// The elsets do not belong to the model. Just delete the container
    elsets.clear();
	
	// free the nodesets
	while ( !nodesets.empty() )
	{
		// only the nodesets created by the user are allocated
		// the model holds a pointer for all, so that they can be accesses from anywhere,
		// but some are "owned" by larger entities, such as bodies
		if ( !nodesets.begin()->second->isInternal() ) 	
			delete nodesets.begin()->second;
		nodesets.erase(nodesets.begin() );
	}	
	
	// free the initial conditions
	while ( !initConditions.empty() )
	{
		delete initConditions.front();
		initConditions.pop_front();
	}	
	
	// free the initial rates
	while ( !initRates.empty() )
	{
		delete initRates.front();
		initRates.pop_front();
	}	
	
	
    // free the poorbodies
    {
        list<poorbody*>::iterator iter = thePoorbodies.begin();
        while (iter != thePoorbodies.end() )
        {
            delete *iter;
            ++iter;
        }
        thePoorbodies.clear();
    }
    
	
	// free the DOFnumberer
	delete theNumberer;
}




void model :: accumulateEnergy(energies &energy, const dofset::evaluation_time& when)
{    
    energy.setZero();
	
	// loop over the parts
    {
        std::list<modelpart*>::iterator iter = theParts.begin();
        while ( iter != theParts.end() )
        {
            // get the energies in one body
            energies partEnergy;
            bool something = (*iter)->integrateEnergy(partEnergy, when);
		
            // and accumulate them into the global energies, in case something has been computed
            if (something) energy += partEnergy;
            ++iter;
        }
    }

	
	// loop over the external forces to include their potencial
    {
        extern double global_tn1;
        list<loading*>:: iterator iter=theLoads.begin();
        while ( iter != theLoads.end() )
        {
            (*iter)->accumulateEnergy(energy, global_tn1);
            ++iter;
        }
    }
    
	DebugMessage("Energies in model accumulated.");
}




void  model :: add(body& b)
{
}


void  model :: add(controlvolume& cv)
{
}

void  model :: add(poorbody& b)
{
	thePoorbodies.push_back(&b);
    theParts.push_back(&b);
}


void  model :: add(meshlessbody& mb)
{
}

// this function sends a constraint reference in the arguments
// the model holds a list of constraint pointers and thus we
// push the memory address of the constraint
void  model :: add(constraint &c)
{
	spconstraints.push_back(&c);
}



void model :: add(contactpair& cp)
{
}


void model :: add(elset &es)
{
    elsets[es.name()] = &es;
}


void model :: add(faceset &fs)
{
    /*
     facesets[ fs.getName() ] = &fs;
     nodesets[ fs.getNodes().getName() ] = &(fs.getNodes());
     */
}



void model :: add(initcondition& c)
{
    cout << "Here "<<endl;
	initConditions.push_back(&c);
}


void model :: add(initrate& c)
{
	initRates.push_back(&c);
}



void model :: add(nodeset &nds)
{
	nodesets[ nds.getName() ] = &nds;
}




void model :: add(loading& load)
{
    theLoads.push_back(&load);
}




// order is very important here
// first we advance in time the theContactpairs. These hold no degrees of freedom, but their
// update might depend on the way the bodies ended up in the last equilibrium instant
// bodies do not get influenced by the way contact pairs are updated so they can be
// updated before or after the theContactpairs.
void model :: advanceInTime(const integrator& i, const double dt)
{
    updateCurrentState(dofset::tn1);
    {
        std::list<modelpart*>::iterator piter = theParts.begin();
        while (piter != theParts.end())
        {
            (*piter)->advanceInTime(i,dt);
            ++piter;
        }
    }
    
    std::map< std::string, elset*>::iterator iter = elsets.begin();
    while (iter != elsets.end())
    {
        (*iter).second->advanceInTime();
        ++iter;
    }

    
    std::map< std::string, nodeset*>::iterator niter = nodesets.begin();
    while (niter != nodesets.end())
    {
        (*niter).second->advanceInTime();
        ++niter;
    }    
}





// returns true if there is a contactpair connecting the two bodies
bool model :: areContactConnected(const modelpart& bd1, const modelpart& bd2) const
{
	return false;
}


/* this function checks the model. We should add more features in the future, for
 * a more comprehensive verification.
 */
bool model :: check(ostream& of)
{
	for_each(thePoorbodies.begin(), thePoorbodies.end(), mem_fun(&poorbody::check) );
	of  << "\n   [ Mesh checked ]";
    return true;
}






// this function is use to monitor if the solution has blowed up
double computeAccelerationNorm(const model& theModel)
{
	double nac2 = 0.0;
    
    std::list<modelpart*>::const_iterator iter = theModel.theParts.begin();
    while (iter != theModel.theParts.end())
    {
        nac2 += (*iter)->accelerationNormSquared();
        ++iter;
    }
    
	return sqrt(nac2);
}




void model :: computeExplicitAcceleration()
{
    for_each(thePoorbodies.begin(), thePoorbodies.end(), mem_fun(&poorbody::computeExplicitAcceleration) );
}




void model :: computeExplicitVelocity()
{
}

void model :: computeExplicitEulerianAcceleration()
{
}



void model :: computeMass()
{
	for_each(theParts.begin(), theParts.end(), mem_fun(&modelpart::computeMass) );    
}




double model :: computeSmallestEigenvalueEstimate()
{
	double numerator   = 0.0;
	double denominator = 0.0;
	
    BOOST_FOREACH(modelpart* p, theParts)
	{
		// compute nodal forces doted with velocity
        BOOST_FOREACH(node* nd, p->nodes)
        {
            numerator   += nd->force().dot( nd->velocity() );
            denominator += nd->velocity().squaredNorm() * nd->getMass();
        }
    }
    
	extern double global_dt;
	double lambda0 = numerator/denominator/global_dt; // dimensions of 1/T^2
	return lambda0;
}




void model :: computeNodalDamping()
{
}




// this function is use to monitor if the solution has blowed up
double computeVelocityNorm(const model& m)
{
	double nac2 = 0.0;
	BOOST_FOREACH(modelpart* p, m.theParts)
	{
        BOOST_FOREACH(node* nd, p->nodes)
        {
            nac2 += nd->velocity().squaredNorm();
        }
	}
	return sqrt(nac2);
}




body& model :: getBody(const string &name) const
{
	list<body*>:: const_iterator iter=theBodies.begin();
    return *(*iter);
}




meshlessbody& model :: getMeBody(const string &name) const
{
	list<meshlessbody*>:: const_iterator iter=theMeshlessbodies.begin();
    return *(*iter);
}




controlvolume& model :: getControlVolume(const string &name) const
{
	list<controlvolume*>:: const_iterator iter=theControlVolumes.begin();
    return *(*iter);
}




modelpart& model :: getPart(const string &name) const
{
	list<modelpart*>:: const_iterator iter=theParts.begin();	
	while (iter != theParts.end() )
	{		
        if ( (*iter)->getName() == name) break;
		++iter;
    }
	
	if ( iter == theParts.end())
	{
		stringstream st;
		st << "No part exists with name " << name << " exists.";
		cout << "No part exists with name " << name << " exists." << endl;
		throw runtime_error(st.str());
	}
    return *(*iter);
}




poorbody& model :: getPoorbody(const string &name) const
{
	list<poorbody*>:: const_iterator iter=thePoorbodies.begin();
	while (iter != thePoorbodies.end() )
	{		
        if ( (*iter)->getName() == name) break;
		++iter;
    }
	
	if ( iter == thePoorbodies.end())
	{
		stringstream st;
		st << "No poorbody exists with name " << name << " exists.";
		throw runtime_error(st.str());
	}
    return *(*iter);
}




// to get an element the program must search into all bodies until it finds it
element& model :: getElement(const int label) const
{	
	element *ret(0);
    
	BOOST_FOREACH(modelpart* p, theParts)
	{
        element&  e = p->getElement(label);
        if ( &e != 0) 
        {
            ret = &e;
            break;   
        }
	}
	
	
	// it might happen that there is no element with the given label
	if (ret == 0) 
	{		
		string msg("No element exists with label ");
		stringstream msg2;
		msg2 << label;
		msg += msg2.str();
		
		throw runtime_error(msg);
	}
	return *ret;
}



elset& model :: getElset(const string &name) const
{
    std::map< std::string, elset* >::const_iterator it = elsets.find(name);
    
    if ( it == elsets.end() )
        throw runtime_error("Elset not found");
    
    return *( it->second );
}



faceset& model :: getFaceset(const string& name) const
{
    std::map< std::string, faceset* >::const_iterator it = facesets.find(name);
    
    if ( it == facesets.end() )
        throw runtime_error("Faceset not found");
    
    return *( it->second );
}




size_t model :: getNElements() const
{
	size_t n = 0;
    
    BOOST_FOREACH( modelpart* b, theParts)
    {
        n += b->elements.size();
    }
    
	return n;
}




size_t model :: getNNodes() const
{
	size_t n = 0;
    
	BOOST_FOREACH( modelpart* b, theParts)
    {
        n += b->nodes.size();
    }
    
	return n;
}




size_t model :: getNEvalspots() const
{
	size_t n = 0;
    
	BOOST_FOREACH( modelpart* b, theParts)
    {
        n += b->evalspots.size();
    }
    
	return n;
}




node& model :: getNode(const int label)
{
	node *ret(0);
	
    BOOST_FOREACH(modelpart* p, theParts)
	{
        if (p->hasNode(label)) ret = &(p->getNode(label));
		if (ret != 0) break;
	}
	
	// it might happen that there is no node with the given label
	if (ret == 0) 
	{
		stringstream  error;
		error << "No node exists with label " << label;
		throw runtime_error(error.str());
	}
	
	return *ret;	
}



const node& model :: getNode(const int label) const
{
    return const_cast<model*>(this)->getNode(label);
}



nodeset& model :: getNodeset(const string& name)
{
    std::map< std::string, nodeset* >::iterator i;
    i = nodesets.find(name);
    
    if ( i == nodesets.end() )
    {
        logger::mainlog << "Error: Nodeset " << name << " not found." << endl;
        throw runtime_error("Nodeset not found");
    }
    
    return *( i->second );
}



const nodeset& model :: getNodeset(const string& name) const
{
    return const_cast<model*>(this)->getNodeset(name);
    
}


int model :: getNParts()	const       
{
    int num = thePoorbodies.size();
    return num;
}


size_t model :: getNVisibleElements() const
{
	size_t n = 0;
	BOOST_FOREACH(modelpart* b, theParts) 
    n += b->getNVisibleElements();
    
	return n;
}




/*  Loop through the linked list of constraints and set the corresponding 
 degrees of freedom */
void model :: imposeInitialConditions()
{
	list<initcondition*>:: iterator iter=initConditions.begin();
	
	while (iter != initConditions.end() )
	{
		node *nd = &getNode( (*iter)->getNodeLabel() );
		nd->setInitialDof( (*iter)->getDof(), (*iter)->getValue() );
		
		++iter;
	}
}





/*  Loop through the linked list of constraints and set the initial rates at the 
 corresponding nodes
 */
void model :: imposeInitialRates()
{
    BOOST_FOREACH(initrate* ir, initRates) ir->impose();
}




void model :: incrementSolution(const integrator& i)
{
	BOOST_FOREACH(modelpart* b , theParts) b->incrementSolution(i);
    updateCurrentState(dofset::tna);
	DebugMessage("Mesh solution incremented");
}




/*
 * The initialization has to be done in the order indicated
 */
void model :: initialize()
{
    // quick return
    if (_isInitialized) return;
    DebugMessage("initializing model");
    
	
	// initialize bodies
	for_each( thePoorbodies.begin(), thePoorbodies.end(), mem_fun(&poorbody::initialize) );
	
	
	// calculate and store maxDofsPerNode
	dofsPerNode = 0;
    {
        list<poorbody*>::iterator biter=thePoorbodies.begin();
        while (biter != thePoorbodies.end() )
        {
            dofsPerNode = std::max<size_t>( dofsPerNode, (*biter)->maxDofsPerNode());
            ++biter;
        }
    }
    
    // initialize loadings
    {
        list<loading*>:: iterator pliter=theLoads.begin();
        while (pliter != theLoads.end() )
        {
            (*pliter)->initialize(*this);
            ++pliter;
        }
        DebugMessage("Loads initialized");
    }
	
    // initialize constraints
	list<constraint*>:: iterator spiter=spconstraints.begin();	
	spiter=spconstraints.begin();	
	while (spiter != spconstraints.end() )
	{
		(**spiter).initialize(*this);
		++spiter;
	}
	DebugMessage("Constraints initialized");
	
    
    // initialize initrates
	list<initrate*>:: iterator iriter=initRates.begin();
	iriter=initRates.begin();
	while (iriter != initRates.end() )
	{
		(**iriter).initialize(*this);
		++iriter;
	}
	DebugMessage("Initial rates initialized");

	
	// initialize contactpairs
	// This initialization is more complex, because, in some cases, when a 
	// contactpair is initialized, more contactpairs are created and thus
	// the vector needs to be resized and the iterators lose meaning
    // initialize the contact pairs, creating the contact elements
	// and adding them to the linked list
	// we cannot use iterators because the initialization might create new contactpairs
	// that invalidate the iterators
    
	
    // set the ID map in free and constrained dofs
    markConstrainedDofs();
	theNumberer->setNextDOFtoAssign(0);
	theNumberer->assignDOFNumbers( theParts.begin(), theParts.end() );
	setNDofs(theNumberer->getLastDOFAssigned()+1);
	for_each( thePoorbodies.begin(), thePoorbodies.end(), mem_fun(&poorbody::markConstrainedElements) );
	
	
#ifdef WITHMPI
	// complete nodal information of boundary nodes with partition info
	map< int, vector<int> >::iterator miter = tmpBoundaryNodes.begin();
	while ( miter != tmpBoundaryNodes.end() )
	{
		node &nn = getNode( miter->first );
		boundaryConnectivity.push_back( make_pair( &nn, miter->second ) );
		++miter;
	}
#endif
	
    
    // parsing operations over the geometry
    parse();
	
    DebugMessage("Mesh initialized");
	_isInitialized = true;	
}






void  model :: info(ostream &of)
{
	of	<< "\n\n\n\n";
	printCentered(of,  "M o d e l    i n f o r m a t i o n");
	of  << "\n"
	<< "\n   Number of bodies           : " << setw(6) << theBodies.size()
    << "\n   Number of poorbodies       : " << setw(6) << thePoorbodies.size()
    << "\n   Number of meshlessbodies   : " << setw(6) << theMeshlessbodies.size()
    << "\n   Number of control volumes  : " << setw(6) << theControlVolumes.size()
	<< "\n   Number of contactpairs     : " << setw(6) << theContactpairs.size()
	<< "\n   Number of elements         : " << setw(6) << getNElements()
    << "\n   Number of evalspots        : " << setw(6) << getNEvalspots()
	<< "\n   Number of nodes            : " << setw(6) << getNNodes()
	<< "\n   Number of elsets           : " << setw(6) << elsets.size()
	<< "\n   Number of nodesets         : " << setw(6) << nodesets.size()
	<< "\n   Number of point loads      : " << setw(6) << theLoads.size()
	<< "\n   Number of point bc         : " << setw(6) << spconstraints.size()
	<< "\n   Number of materials        : " << setw(6) << material::getNMaterials()
	<< "\n   Number of elem. types      : " << setw(6) << eltype::getNEltypes()
	<< "\n   Number of unknowns         : " << setw(6) << nDofs
	<< "\n   Max. dofs/node             : " << setw(6) << dofsPerNode << endl;
}





/* extracts from the global increment vector the data corresponding to
 each element and stores it in the nodal data */
void model :: localizeSolutionIncrement(const longvector& delta,  partitionDofs pt)
{
    BOOST_FOREACH(modelpart* p, theParts)
	{
		p->localizeSolutionIncrement(delta, pt);
	}
}




/* loops through the constraints in the model creating calling their corresponding
 nodes and marking their constrained dofs as CONSTRAINED
 *
 * Glued nodes which depend on other as marked as constrained at this point, so that
 * they are not counted as free dofs, and they are not given ID
 */
void model :: markConstrainedDofs()
{
	list<constraint*>:: iterator iter=spconstraints.begin();	
	while (iter != spconstraints.end() )
	{
		(**iter).markConstrainedDofs();
		++iter;
	}	
}




double model :: maxEigenvalueEstimate() const 
{
    double  maxev(0.0);
	
    BOOST_FOREACH(modelpart* p, theParts)
    {
		double eigv = p->maxEigenvalueEstimate();
		if (maxev < eigv)
		{
			maxev = eigv;
		}
	}	
	
	if (maxev == 0) logger::warnings << "Maximum model eigenvalue is zero." << endl;
	
    return maxev;
}




/* message indicates if there must be a logged message with the file name */
bool model :: openModelFile(const commandLine &cl, char *message, ifstream &modelfile)
{
    /* the first command is modelfile = filename */
    const usercommand &uc = cl[0];
	modelfilename = uc.option();
	
    // open the file with the corresponding access mode */
    modelfile.open(modelfilename.c_str());
	
    /* Log message, error if file not open */
    if (modelfile.fail())  ErrorMessage("no model file named '%s' exists. Exiting F EL I K S", modelfilename.c_str() );
    if (strlen(message)>0) logger::toLogFile2(message, "%s", modelfilename.c_str() );
	
    return true;
}



/*  this function gathers several operations that one does after all the geometric
 *   entities have been defined and collected.
 *   Right now it does:
 *      -- parse the bodies
 *
 */
void model :: parse()
{
    //quick return
	if (_isParsed) return;
	
	logger::PARSEtime.start();
//	for_each(theBodies.begin(), theBodies.end(), mem_fun(&body::parse));
	logger::PARSEtime.stop();
    
    _isParsed = true;
}




void model :: print(const std::string& what, ostream &of) const
{
    if ( (what == "bodies" || what == "all") && !theBodies.empty() )
    {
        stringstream title;
        title << "B o d i e s   (" << theBodies.size() << ")";
        string str(title.str());
        of << "\n\n\n";
        printCentered(of, str);
        
        list<body*>::const_iterator iter=theBodies.begin();
        while (iter != theBodies.end() )
        {
           // (*iter)->print(of);
            ++iter;
        }
        of << "\n";
    }
    
    if ( (what == "bodies" || what == "all") && !thePoorbodies.empty() )
    {        
        stringstream title;
        title << "P o o r b o d i e s   (" << thePoorbodies.size() << ")";
        string str(title.str());
        of << "\n\n\n";
        printCentered(of, str);
        
        list<poorbody*>::const_iterator piter=thePoorbodies.begin();
        while (piter != thePoorbodies.end() )
        {
            (*piter)->print(of);
            ++piter;
        }
        of << "\n";
    }
    
    
	if (_mustPrint || what == "nodes" || (what == "all" && Debug_Level() > 0))
	{
        stringstream title;
        title << "N o d e s  (" << getNNodes() << ")" << "\n";
        string str(title.str());
        of << "\n\n\n";
        printCentered(of, str);
        of << "\n";
        
        list<modelpart*>::const_iterator iter = theParts.begin();
        while( iter != theParts.end() )
        {
            printNodes(**iter, of);
            ++iter;
        }	        
    }
    
    if (_mustPrint || what == "elements" || (what == "all" && Debug_Level() > 0))
    {
        stringstream title;
        title << "E l e m e n t s (" << getNElements() << ")";
        string str(title.str());
        of << "\n\n\n";
        printCentered(of, str);
        of << "\n";
        
        list<body*>::const_iterator iter = theBodies.begin();
        while( iter != theBodies.end() )
        {
          //  printBodyElements(**iter, of);
            ++iter;
        }
        
        of << "\n";
    }
    
    if ((_mustPrint || what == "nodesets" || what == "all") && !nodesets.empty())
    {
     	of << "\n\n\n";
        stringstream title;
        title << "N o d e    S e t s   (" << nodesets.size() << ")";
        string str(title.str());
        printCentered(of, str);
        of << "\n";
        
        map< string, nodeset* >::const_iterator iter = nodesets.begin();
        while ( iter != nodesets.end() )
        {
            iter->second->print(of);
            ++iter;
        }
        of << endl;
    }
    
    if ((_mustPrint || what == "elsets" || what == "all") && !elsets.empty())
    {
     	of << "\n\n\n";
        stringstream title;
        title << "E l e m e n t    S e t s   (" << elsets.size() << ")" << "\n";
        string str(title.str());
        printCentered(of, str);
        of << "\n";
        
        map< string, elset*>::const_iterator iter=elsets.begin();
        while (iter !=  elsets.end())
        {
            iter->second->print(of);
            ++iter;
        }
        of << "\n";
    }
    
    /*
     if ((mustPrint || what == "facesets" || what == "all") && !facesets.empty())
     {
     int ntotalFacesets=0;
     list<body*>::const_iterator iter = theBodies.begin();
     while (iter != theBodies.end())
     {
     ntotalFacesets += (*iter)->getFacesets().size();
     ++iter;
     }
     
     of << "\n\n\n";
     stringstream title;
     title << "F a c e    S e t s   (" << ntotalFacesets << ")";
     string str(title.str());
     printCentered(of, str);	
     if (ntotalFacesets == 0) return;
     
     iter = theBodies.begin();
     while (iter != theBodies.end())
     {
     of << "\n";
     printBodyFacesets(**iter, of);
     ++iter;
     } 
     
     of << "\n";
     
     }
     */
    
    if ((_mustPrint || what == "constraints" || what == "all") && !spconstraints.empty())
    {
        of << "\n\n\n";
        stringstream title;
        title << "C o n s t r a i n t s   (" << spconstraints.size() << ")";
        string str(title.str());
        printCentered(of, str);
        
        list<constraint*>::const_iterator iter=spconstraints.begin();
        while (iter != spconstraints.end() )
        {
            (*iter)->print(of);
            ++iter;
        }
        of << "\n";        
    }
    
    if ((_mustPrint || what == "loads" || what == "all") && !theLoads.empty())
    {
        of << "\n\n\n";
        stringstream title;
        title << "L o a d s  (" << theLoads.size() << ")";	
        string str(title.str());
        printCentered(of, str);
        of << "\n";
        
        BOOST_FOREACH(loading *ld, theLoads)
            ld->print(of);
        
        of << endl;
    }
    
    if ((_mustPrint || what == "initialconditions" || what == "all") && !initConditions.empty())
    {
        of << "\n\n\n";
        stringstream title;
        title << "I n i t i a l   C o n d i t i o n s   (" << initConditions.size() << ")";
        string str(title.str());
        printCentered(of, str);
        of << "\n\n";
        
        of << " node     dof     value   ignore";
        
        list<initcondition*>::const_iterator iter = initConditions.begin();	
        while (iter != initConditions.end() )
        {
            (*iter)->print(of);
            ++iter;
        }
        of << "\n";
        
    }
    
    if ((_mustPrint || what == "initialrates" || what == "all") && !initRates.empty())
    {
        of << "\n\n\n";
        stringstream title;
        title << "I n i t i a l   R a t e s   (" << initRates.size() << ")"; 
        string str(title.str());
        printCentered(of, str);
        
        
        list<initrate*>::const_iterator iter=initRates.begin();
        while (iter != initRates.end() )
        {
            (*iter)->print(of);
            ++iter;
        }
        of << "\n";

    }
    of << flush;
}




void model :: gatherCurrentSolution(std::vector<double>& sol) const
{
    // for mmonca interface
    sol.resize(getNNodes()*3);
    
    size_t a=0;
    BOOST_FOREACH( modelpart* p, theParts)
    {
        BOOST_FOREACH( node* nd, p->nodes)
        {
            ivector& u = nd->getUDS().displacement();
            sol[a++] = u[0];
            sol[a++] = u[1];
            sol[a++] = u[2];
        }
	}
}




// displays in screen and log file the nodes and their dofs
void model :: printCurrentSolution(ostream &of) const 
{
    of  << "\n\n                   Current Nodal solution"
        << "\n                   ----------------------\n\n";
	
	BOOST_FOREACH( modelpart* p, theParts)
    {
        BOOST_FOREACH( node* nd, p->nodes)
        {
            of << right << setw(4) << nd->getLabel();
            nd->printDOFs(of);
            of << "\n";
        }
	}
    of << endl;
}





void model :: printDetailedSolution(ostream &of) const
{
	of << "\n\n\n";
	stringstream title;
	title << "N o d e s  (" << getNNodes() << ")" << "\n";
	string str(title.str());
	printCentered(of, str);
	of << "\n";
	
	BOOST_FOREACH( modelpart* p, theParts)
    {
        BOOST_FOREACH( node* nd, p->nodes)
        {
            of << *nd;
        }
    }
    of << endl;
}




// prints the gradients of an element, given its label
// if label is <=0, print all gradients
void model :: printGradients(int elabel, ostream& os) const
{
    element &elmt(getElement(elabel));
    elmt.printGradients(os);
}




void model :: printNodalForces(ostream &of) const
{
	BOOST_FOREACH( modelpart* p, theParts)
    {
        BOOST_FOREACH( node* n, p->nodes)
        {
            of << "\n Label " << n->getLabel() << "  " << n->force() ;	
        }
    }
}




void  model :: printReactions(ostream& of) const
{
    of << "\n\n           Nodal reactions (R = fint)";
    of << "\n          ----------------------------";
    
	BOOST_FOREACH( modelpart* p, theParts)
    {
        size_t    dofs = p->maxDofsPerNode();
        double *ftot  = new double[dofs]();
        
        BOOST_FOREACH( node* n, p->nodes)
        {
            of << "\n Label " << n->getLabel() << ": ";
            for (int i=0; i<dofs; i++)
            {
                of <<  -(n->force()(i));
                ftot[i] += n->force()(i);
            }
        }
        
        of << "\n          ----------------------------------";
        of << "\n Sum =        ";
        for (int i=0; i<dofs; i++) 
            of << -ftot[i];
        
        delete []ftot;
    }
}    




void  model :: printReactionsInConstraints() const
{
//	BOOST_FOREACH(const body* bd, theBodies)
    {
//        printBodyReactionsInConstraints(*bd);
    }
}




bool model :: recompute(bool& ghostnodesappeared)
{	
    return false;
}




bool  model :: remove(contactpair &cp)
{
	size_t original( theContactpairs.size() );
	
	theContactpairs.erase(std::remove( theContactpairs.begin(), theContactpairs.end(), &cp ), theContactpairs.end() );
	
	// return found if at least an element has been deleted
    return original > theContactpairs.size();
}




/*  This function only removes the elset from the list in the model, it does not
 delete the elset itself.
 */
bool model :: remove(elset& es)
{
    return ( elsets.erase(es.name()) == 1 );
}




/*  This function only removes the faceset from the list in the model, it does not
 delete the nodeset itself.
 */
bool model :: remove(faceset& nds)
{
    /*
     size_t original( facesets.size() );
     
     //remove from the list the nodeset that has the name of nds
     facesets.erase( nds.getName() );
     
     // return found if at least an element has been deleted
     return original > facesets.size();
     */
    return false;
}




/*  This function only removes the nodeset from the list in the model, it does not
 delete the nodeset itself.
 */
bool model :: remove(nodeset& nds)
{
    return ( nodesets.erase(nds.getName()) == 1 );
}




// set all the nodal masses to zero
void  model :: resetNodalMasses()
{
//	for_each( theBodies.begin(), theBodies.end(), mem_fun(&body::resetNodalMasses) );
	for_each( thePoorbodies.begin(), thePoorbodies.end(), mem_fun(&poorbody::resetNodalMasses) );
}




// set all the nodal forces to zero
void  model :: resetNodalForces()
{
//	for_each( theBodies.begin(), theBodies.end(), mem_fun(&body::resetNodalForces) );
	for_each( thePoorbodies.begin(), thePoorbodies.end(), mem_fun(&poorbody::resetNodalForces) );
}




/* this function stores the solution at time tn in the variables at time tn+1.
 It is useful when, after trying to converge a Newton Raphson solution, we arrive
 to a point where the data at time tn+1 is too bad to converge, and we rather start
 again from the beginning.
 This function rewinds the nodal degrees of freedom, the element internal variables,
 and the history varibles
 */
void model :: rewindSolution()
{
//	for_each( theBodies.begin(),     theBodies.end(),     mem_fun(&body::rewindSolution) );
	for_each( thePoorbodies.begin(), thePoorbodies.end(), mem_fun(&poorbody::rewindSolution) );
 }




// set all the nodal solution increments dU to zero
void  model :: resetSolutionIncrements()
{
//	for_each( theBodies.begin(),     theBodies.end(), mem_fun(&body::resetSolutionIncrements) );
	for_each( thePoorbodies.begin(), thePoorbodies.end(), mem_fun(&poorbody::resetSolutionIncrements) );
}






/*  Format 
 *   directors, nodeset = name |all,
 *              input = [spherical | cyclindrical | fixed ]
 *              d1=x, d2=x, d3=x,               (if type is "director")
 *              q1=x, q2=x, q3=x, q4=x          (if type is "quaternion")
 *              m11=x, m12=x, m13=x, ...        (if type is "matrix")
 */
void model :: scanDirectors(const commandLine &cl)
{
    int generate=-1;
	ivector		 dir, pp;
	nodeset      *nds(0);
	bool         allNodes(false);
    
    // The first command is "directors" and is mandatory.
	
	// Scan the rest of command
	for (int i=1; i< cl.size(); i++)
    {
        const usercommand &uc = cl[i];
		if      (uc.keyword() == "nodeset" && uc.option() == "all") allNodes = true;
		else if (uc.keyword() == "nodeset") nds = &getNodeset(uc.option());
		
		else if (uc.keyword() == "input") 
		{
			// command is "input = [director | quaternion | matrix ] and is mandatory
			if      (uc.option() == "fixed" )		generate = 0;
			else if (uc.option() == "spherical")	generate = 1;
			else if (uc.option() == "cylindrical")	generate = 2;
			else logger::warnings << "Unknown input coordinate system for rotations: " << uc.option() ;
		}
		
		else if (uc.keyword() == "d1") dir[0] = uc.value();
		else if (uc.keyword() == "d2") dir[1] = uc.value();
		else if (uc.keyword() == "d3") dir[2] = uc.value();
		else if (uc.keyword() == "p1") pp[0]  = uc.value();
		else if (uc.keyword() == "p2") pp[1]  = uc.value();
		else if (uc.keyword() == "p3") pp[2]  = uc.value();
		
	}
	
	// quick and dirty fix
	if (allNodes)
	{
		nds = new nodeset("allNodes");
		list<modelpart*>::iterator iter = theParts.begin();
		while( iter != theParts.end() )
		{
            BOOST_FOREACH( node* nd, (*iter)->nodes)
            nds->add( *nd);
			++iter;
		}
		
		nds->initializeRotations(_UDnode, generate, dir, pp);
	}
	
	else if (nds != 0) 
		nds->initializeRotations(_UDnode, generate, dir, pp);
}







// when scanning nodes, put then in a temporary map and only create the nodes when the
// elements are created, because the type of node is not known until it is associated
// with an element
void model :: scanBoundaryNodes(const commandLine &cl, ifstream &modelfile)
{
#ifdef WITHMPI
	int     label;
	
    // scan vertex coordinates until blank line or EOF
    string oneLine;
 	while ( commandLine::scanCommandLine(modelfile, oneLine) && oneLine.size() > 0)
	{
		stringstream str(oneLine);
		
		// read label and type
		str >> label;
		
		// read node labels
		int k=0, partitionLabel;
		while (str >> partitionLabel )
		{
			tmpBoundaryNodes[label].push_back( partitionLabel );
			k++;
		}
		
		if (DEBUG_SCANNER) cout << "|" << str.str() << "|" << 
			"npart = " << tmpBoundaryNodes[label].size()   <<
			"first = " << tmpBoundaryNodes[label][0] << "\n";
	}
#endif
}




void model :: scanElements(const commandLine &cl, ifstream &modelfile)
{
    std::string bodyname("Default");
    
    // Scan the rest of command
	for (int i=1; i< cl.size(); i++)
    {
        const usercommand &uc = cl[i];
		if      (uc.keyword() == "body") bodyname = uc.option();
    }
    
    modelpart& thepart( this->getPart(bodyname) );
    dynamic_cast<poorbody&>(thepart).scanElements(cl, modelfile);
}




// when scanning nodes, put then in a temporary map and only create the nodes when the
// elements are created, because the type of node is not known until it is associated
// with an element
void model :: scanNodes(const commandLine &cl, ifstream &modelfile)
{
    std::string bodyname("Default");
    
    // Scan the rest of command
	for (int i=1; i< cl.size(); i++)
    {
        const usercommand &uc = cl[i];
		if      (uc.keyword() == "body") bodyname = uc.option();
    }
    
    modelpart& thepart( this->getPart(bodyname) );
    dynamic_cast<poorbody&>(thepart).scanNodes(cl, modelfile);
}




void model :: setNodalMasses(const longvector& v)
{
	list<body*>::iterator iter = theBodies.begin();
	while( iter != theBodies.end() )
	{
//		(*iter)->setNodalMasses(v);
		++iter;
	}				
}




void model :: updateCurrentState(const dofset::evaluation_time& when)
{    
    feliks_for_each(theParts.begin(), theParts.end(), 
                  std::bind2nd(std::mem_fun(&modelpart::updateCurrentState), when));    
}




void model :: updateElementLabels(const int label)
{
    smallest_elmt = std::min<int>(smallest_elmt, label);
}




void model :: updateNodeLabels(const int label)
{
    smallest_node = std::min<int>(smallest_node, label);
}


