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
//
//  poorbody.h
//  feliks
//
//  Created by Romero Ignacio on 4/2/12.
//  Copyright (c) 2012 Universidad Politécnica de Madrid. All rights reserved.
//
//  a poorbody is a modelpart that seems like a body, but lacks the data
//  structures to be such. This is the the type of parts that one obtains
//  when the input data is a list of nodes and elements, for example,
//  in the old fashion FE style. A poorbody does not have any topological
//  features, which can be a limitation for refinement or for accessing
//  some of its boundaries. To mitigate this limitation, each poorbody
//  can define nodesets, elsets, and facesets, but always in a more
//  costly manner

#pragma once
#ifndef feliks_poorbody_h
#define feliks_poorbody_h


#include "Model/Parts/modelpart.h"
#include "Model/Sets/nodeset.h"
#include "Model/Sets/elset.h"
#include "Model/Sets/faceset.h"
#include "Math/tensor.h"

#include <iostream>
#include <list>
#include <string>
#include <map>
#include <set>


class commandLine;
class integrator;
class node;
class nodeset;


namespace feliks
{
    namespace topology
    {
        class cellcomplex;
    }
}



class poorbody : public modelpart
{
    
public:
    poorbody(const std::string& name, model& m);
                                            
    poorbody(const std::string& name, const std::string& geo, model& m, const size_t dim);
	virtual                                 ~poorbody();
    
    // scan the definition of the body from a commandline of the type 
	static poorbody&                        scan(const commandLine &cl, model &m);
	void                                    scanElements(const commandLine &cl, std::ifstream &modelfile);
    void                                    scanNodes(const commandLine &cl, std::ifstream &modelfile);
    
    virtual void                            advanceInTime(const integrator& integ, const double dt);
    virtual bool                            check() const;
    virtual double                          getMass() const;
    virtual void                            incrementSolution(const integrator& i);
    virtual void                            initialize();
    virtual void                            parse();
    virtual void                            parseInformationReset();
    virtual double                          maxEigenvalueEstimate() const;
    virtual void                            print(std::ostream& os=std::cout);
    
	
    virtual size_t                          dimension() const;
    const faceset&                          getExternalFaces()	const;
	faceset&                                getFaceset(const std::string& fname);
    std::set<element*>&                     getElementsTouchingNode(node& nd);
	void                                    resetSolutionIncrements();
    void                                    resetNodalMasses();
	void                                    resetNodalForces();
	void                                    rewindSolution();
    void                                    updateCurrentState(const dofset::evaluation_time when);
    
    
    virtual void                            computeMass();
    virtual void                            computeExplicitAcceleration();
	virtual void                            computeExplicitVelocity();
    void                                    computeExplicitEulerianAcceleration();

    
    // get data
    virtual std::set<node*>                     getNodesFrom0Cells(const set<feliks::topology::cell*>& cells) const;
    virtual feliks::meshing::solidmesh&           getTheSolidmesh();
    virtual const feliks::meshing::solidmesh&     getTheSolidmesh() const;

    
    std::map<std::string,faceset*>          theFacesets;    
protected:
    

    std::list<nodeset*>                     theNodesets;
    std::list<elset*>                       theElsets;    
    nodeset                                 allNodes;
    nodeset                                 externalNodes;
    faceset                                 externalFaces;
    void                                    computeExternalNodesAndFaces();
    int                                     remeshafterthesesteps;
    
    
private:
    
    std::map< int, std::pair<blue::ivector, node*> >          tmpNodes;
    std::string                             _geometricDescription;
    bool                                    _isParsed;
    bool                                    _hasVariableMass;
    bool                                    _massComputed;
    bool                                    _externalNodesFound;
    size_t                                  _dimension;
    double                                  mass;
    size_t                                  coordinationNumber;
    std::map<node*,std::set<element*> >     elementsTouchingNode;
    
    
    friend class meshlessbody;
};


#endif
