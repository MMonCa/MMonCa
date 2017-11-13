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
 * model.h  
 *
 * the model holds the entities responsible for contributing to the equations
 * that are subsequently solved. 
 *
 * iro, november 1999
 * converted to C++ jan 2006
 */


#pragma once
#ifndef _feliks_model_h
#define _feliks_model_h

#include <cstdio>
#include <fstream>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>


#include "Analysis/Dofs/dofnumberer.h"
#include "Analysis/Dofs/dofset.h"
#include "Math/tensor.h"


class body;
class constraint;
class contactpair;
class commandLine;
class DOFnumberer;
class element;
class elset;
class evaluation_time;
class energies;
class faceset;
class initcondition;
class initrate;
class integrator;
class loading;
class longvector;
class node;
class nodeset;
class poorbody;
class meshlessbody;
class controlvolume;

namespace feliks
{
    class mmoncaInterface;
};



//enum  MEcontactType{ MEtoME , MEtoCV};


class model
{

public:
		
	model();
	~model();
        
	
	// tasks that must be done at the beginning
	bool                            check(std::ostream &of=std::cout);
	void                            fill(const feliks::mmoncaInterface& mi);
	void                            initialize();
	inline void                     setNDofs(const int n);
	inline void                     setStaggeredNumberer();

	
	// functions to add/remove members to the model
	void                            add(body& bd);
	void                            add(constraint& cons);
	void                            add(contactpair& cp);
    void                            add(controlvolume& cv);    
	void                            add(elset& es);
	void                            add(faceset& es);
	void                            add(initcondition& c);
	void                            add(initrate& c);
	void                            add(nodeset &nds);
	void                            add(poorbody& bd);
    void                            add(meshlessbody& mb);
	void                            add(loading& pload);
	
	bool                            remove(contactpair &cp);
	bool                            remove(elset& es);
	bool                            remove(faceset& es);
	bool                            remove(nodeset& nds);

    

	// functions to retrieve individual members of the model
	DOFnumberer&                    theDOFNumberer();
    poorbody&                       getPoorbody(const std::string &name) const;
	body&                           getBody(const std::string &name) const;
	controlvolume&                  getControlVolume(const std::string &name) const;
	meshlessbody&                   getMeBody(const std::string &name) const;
    modelpart&                      getPart(const std::string &name) const;
    element&                        getElement(const int label) const;
	elset&                          getElset(const std::string &name) const;
	faceset&                        getFaceset(const std::string& name) const;
	node&                           getNode(const int label);
	const node&                     getNode(const int label) const;
    nodeset&                        getNodeset(const std::string& name);
	const nodeset&                  getNodeset(const std::string& name) const;


    
    
    // the model is, essentially, the container of parts. Then, bodies
    // and contactpairs are also collected independently for easy access,
    // but are duplicated here. Only the list theParts owns the parts
    std::list<modelpart*>           theParts;
    std::list<body*>                theBodies;
    std::list<poorbody*>            thePoorbodies;
    std::list<meshlessbody*>        theMeshlessbodies;
    std::list<controlvolume*>       theControlVolumes;
    std::list<contactpair*>         theContactpairs;


    

	// functions to get information about the model
    bool                            areContactConnected(const modelpart& bd1, const modelpart& bd2) const;

	double                          computeSmallestEigenvalueEstimate();
	inline size_t                   getDofsPerNode() const		{return dofsPerNode;}
	inline std::string              getFilename() const			{return modelfilename;}
	inline double                   getMass() const;
	inline size_t                   getNDofs()	const;
    size_t                          getNElements() const;
    size_t                          getNEvalspots() const;
	size_t                          getNNodes() const;
	size_t                          getNVisibleElements() const;
	int                             getNParts()	const;
    double                          maxEigenvalueEstimate() const;
	inline int                      internalNodeLabel()	const	{return smallest_node;}
	inline int                      internalElmtLabel()	const	{return smallest_elmt;}
    bool                            isMEinterpolation(const commandLine& cline);
    
	
	// functions to print data of the model
	void                            info(std::ostream &of=std::cout);
	void                            print(const std::string& what, std::ostream &of=std::cout) const;
	void                            printCurrentSolution(std::ostream &of=std::cout) const;
	void                            printDetailedSolution(std::ostream &of=std::cout) const;
	void                            printGradients(int n, std::ostream& of=std::cout) const;   // n=-1 -> all
	void                            printNodalForces(std::ostream &of=std::cout) const;
	void                            printReactions(std::ostream& of=std::cout) const;
	void                            printReactionsInConstraints() const;
	
    
	
	// compute over the whole model
    void                            computeMass();
	void                            accumulateEnergy(energies& globalEnergies, const dofset::evaluation_time& when);
    void                            updateCurrentState(const dofset::evaluation_time& when);
    void                            gatherCurrentSolution(std::vector<double>& sol) const;
    
	
	// functions to modify the solution data
	void                            advanceInTime(const integrator& i, const double dt);
	void                            computeExplicitAcceleration();
    void                            computeExplicitEulerianAcceleration();
	void                            computeExplicitVelocity();
	void                            computeNodalDamping();
	void                            imposeInitialConditions();
	void                            imposeInitialRates();
	void                            incrementSolution(const integrator& i);
	void                            localizeSolutionIncrement(const longvector& delta,  partitionDofs pt=PARTITION_ALL);
	bool                            recompute(bool& ghostnodesappared);
	void                            resetNodalMasses();
	void                            resetNodalForces();
	void                            resetSolutionIncrements();
	void                            rewindSolution();
	void                            setNodalMasses(const longvector& v);
	void                            updateElementLabels(const int label);
	void                            updateNodeLabels(const int label);

    
    
private:
	
	static const int                                    ndm=3;
    size_t                                              nDofs;          // total # dofs*/
    size_t                                              dofsPerNode;    // maximum
    std::string                                         modelfilename;
    DOFnumberer*                                        theNumberer;
	
	double                                              mass;
    bool                                                _isInitialized;
    bool                                                _isParsed;
    bool                                                _mustPrint;
    
	int                                                 smallest_elmt;
	int                                                 smallest_node;
		
        
    
    // these are global dictionaries. The model does not own the
    // data but enables global retrieval
    std::map<std::string, elset*>                       elsets;
    std::map<std::string, nodeset*>                     nodesets;
    std::map<std::string, faceset*>                     facesets;
    
    std::list<initcondition*>                           initConditions;
    std::list<initrate*>                                initRates;
    std::list<loading*>                                 theLoads;
    std::list<constraint*>                              spconstraints;
    
    
#ifdef WITHMPI    
    std::map< int, std::vector<int> >                   tmpBoundaryNodes;
    std::list< std::pair<node*, std::vector<int> > >    boundaryConnectivity;
#endif
	
	friend class assembler;
#ifdef WITHMPI
	friend class MPIExplicitAssembler;
#endif
	friend class FEanalysis;
	friend class stepsolver;
	friend class contactpair;
	friend class centraldiff;
	friend class asynchronousAnalysis;
	friend class xpartitioner;
	
	// these functions are only for internal use
	void markConstrainedDofs();
	bool openModelFile(const commandLine &cl, char *message, std::ifstream &modelfile);
	void parse();

	
	// functions to read data from model file
	void fillWithFileCommands( const std::string& modelfilename);
	void scanBoundaryNodes(    const commandLine &cl, std::ifstream &modelfile);
	void scanConstraints(      const commandLine &cl, std::ifstream &modelfile, model &m);
	void scanDirectors(        const commandLine &cl);
	void scanElements(         const commandLine &cl, std::ifstream &modelfile);
	void scanInitialConditions(const commandLine &cl, std::ifstream &modelfile);
	void scanNodes(            const commandLine &cl, std::ifstream &modelfile);	
};	



inline double                           model :: getMass() const          {return mass;}
inline size_t                           model :: getNDofs()	const         {return nDofs;}
inline DOFnumberer&                     model :: theDOFNumberer()         {return *theNumberer;}
inline void                             model :: setNDofs(const int n)	  {nDofs = n;}
inline void                             model :: setStaggeredNumberer()	  {theNumberer = new splitDOFnumberer();}




double    computeAccelerationNorm(const model& theModel);
double    computeVelocityNorm(const model& theModel);







#endif
