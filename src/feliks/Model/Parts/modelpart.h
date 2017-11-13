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
 *  modelpart.h
 *  feliks
 *
 *  Created by Ignacio Romero on 4/11/08.
 *  Copyright 2008 Universidad Politécnica de Madrid. All rights reserved.
 *
 *
 * the modelpart is a class created to hold bodies and contactpairs, sharing the functions that are
 * common. A modelpart might not have any nodes, as in the case of a contactpair
 */


#pragma once
#ifndef _feliks_modelpart_h
#define _feliks_modelpart_h


#include <stdexcept>
#include <vector>
#include <list>
#include <map>

#include "Analysis/Dofs/dofset.h"
#include "Math/tensor.h"
#include "Model/Sets/nodeset.h"

class assembler;
class evalspot;
class elset;
class element;
class energies;
class integrator;
class longvector;
class model;
class node;
class resultdata;


namespace feliks
{
    namespace topology
    {
        class cell;
    }
    
    namespace geometry
    {
        class BRep;
    }
        
    namespace meshing
    {
        class solidmesh;
    }
}



class modelpart
{
    
public:
    modelpart(const std::string& name, model& m);
	virtual                     ~modelpart();
    
    // definition
    virtual void                initialize()=0;
    void                        markConstrainedElements();
    void                        rotateAndShift();
    void                        setName(const std::string& newname);
    
    
	// information
    virtual size_t              dimension() const = 0;
    std::string                 getName() const;
    int                         getNVisibleElements() const;
    bool                        hasNode(const int label) const;
    virtual double              maxEigenvalueEstimate() const = 0;
	size_t                      maxDofsPerNode() const;
    
	
	// make active or inactive
	void                        activate();
	void                        deactivate();   
	bool                        isActive() const;
	
    void                        scheduleForDeletion();
    bool                        isScheduledForDeletion() const;
	
    void                        setInitialized();
    bool                        isInitialized() const;
    
    
	// retrieve data
	elset&						getCompleteElset();
    element&					getElement(const int label) const;
    virtual double              getMass() const=0;
    nodeset&                    getExternalNodes();
    const nodeset&              getExternalNodes() const;
    node&                       getNode(const int label) const;
	virtual	void                result(resultdata& r){}
    
    // get data
    virtual std::set<node*>                     getNodesFrom0Cells(const set<feliks::topology::cell*>& cells) const = 0;
    virtual feliks::meshing::solidmesh&           getTheSolidmesh() = 0;
    virtual const feliks::meshing::solidmesh&     getTheSolidmesh() const = 0;
    
	
    void                        add(const element& e);
    void                        add(const evalspot& ev);
    void                        add(const node& nd);
    
    bool                        remove(element& e);
    bool                        remove(node& n);
    
    
	// operations on the data
    double                      accelerationNormSquared() const;
    void                        changeSignOfNodalForces();
    virtual void                computeExplicitAcceleration()=0;
	virtual void                computeExplicitVelocity()=0;

    virtual void                computeMass()=0;
	void                        localizeSolutionIncrement(const longvector& delta, partitionDofs pt=PARTITION_ALL);
	void                        replaceNode(node& deleted, node& survivor);
	void                        sortElements();
    void                        sortNodes();
    
    
    
    // virtual operations
    virtual void                advanceInTime(const integrator& i, const double dt) = 0;
    virtual bool                check() const;
    virtual void                parse()=0;
    virtual void                parseInformationReset()=0;
    virtual void                print(std::ostream&of = std::cout);
    
    virtual void                incrementSolution(const integrator& i) = 0;
    virtual bool                integrateDampingMatrix(assembler& as) const;
    virtual bool                integrateDEnergy() const;
    virtual bool                integrateEnergy(energies& energy, const dofset::evaluation_time& when) const;
    virtual bool                integrateMassMatrix(assembler& as) const;
    virtual void                updateCurrentState(const dofset::evaluation_time when) = 0;
    
    
    std::vector<node*>          nodes;
    std::vector<element*>       elements;
    std::vector<evalspot*>      evalspots;
    
    std::map<const feliks::topology::cell*, const node*> cell2node;
    std::map<const feliks::topology::cell*, const element*> cell2element;
    
    
protected:
    
    void                        colorDisjointElements();
    
    
    std::string                 _name;
	bool                        active;
    bool                        _elementsSorted;
	bool                        _nodesSorted;
    nodeset                     _externalNodes;
    
	model&                      linkedModel;
	bool                        _scheduledForDeletion;
    blue::ivector                     center, orientation;
    
	
private:
    modelpart();
    
    bool                        _isInitialized;
	elset*                      _theElements;
    double                      _mass;
    
    
	friend	class	assembler;
	friend	class	radicalLumperAssembler;
	friend	class	tbbResidualAssembler;
	friend	class	MPIImplicitAssembler;
    friend  class   meshlessbody;
	
};


inline void         modelpart :: activate() {active = true;}
inline void         modelpart :: deactivate(){active = false;}
inline std::string  modelpart :: getName() const {return _name;}
inline void         modelpart :: scheduleForDeletion(){ _scheduledForDeletion = true;}
inline bool         modelpart :: isActive() const	{return active;}
inline bool         modelpart :: isInitialized() const {return _isInitialized;}
inline bool         modelpart :: isScheduledForDeletion() const { return _scheduledForDeletion;}
inline void         modelpart :: setInitialized(){_isInitialized = true;}
inline void         modelpart :: setName(const std::string& newname) {_name = newname;}



void                printNodes(const modelpart& b, std::ostream &of);




#endif





