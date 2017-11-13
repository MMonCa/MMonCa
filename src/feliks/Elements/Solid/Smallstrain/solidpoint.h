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
 * solidpoint.h
 *
 * solid point, small strain
 * iro, september 2011
 */

#pragma once
#ifndef _feliks_solidpoint_h
#define _feliks_solidpoint_h

#include "Elements/evalspot.h"
#include "Elements/element.h"
#include "Elements/Interpolation/interpolation.h"
#include "Math/matrix.h"
#include "Math/tensor.h"
#include "Math/vector.h"
#include "Model/Node/node.h"
#include "Main/feliksinterface.h"

#include <iostream>
#include <vector>

class assembler;
class commandLine;
class energies;
class quadpoint;
class resultdata;
class smallStrainMaterial;
class smallStrainMP;




class solidEvalPoint : public evalspot
{
 public:
                                solidEvalPoint(const double                     vol, 
                                               const smallStrainMaterial&       mat,
                                               const std::vector<shapefun>&     sh,
                                               const feliks::elementState&        es);
    
                                solidEvalPoint(const double                     vol, 
                                               const smallStrainMaterial&       mat,
                                               const std::vector<shapefun>&     sh,
                                               const std::vector<shapefun>&     shBar);
    
    virtual                     ~solidEvalPoint();
    
    
    // dof information
    void                        getDofNumbers(vector<int>& theDofs) const;
    bool                        hasConstrainedDofs();
    virtual double              mass() const;

    
    bool                        check() const;
    blue::ivector               coordinates( const dofset::evaluation_time when=dofset::tna) const;
    blue::ivector               getReferenceCoordinates() const;
    void                        computeStrain(const dofset::evaluation_time when, blue::istensor& eps) const;

    bool                        genericIntegral(assembler& theAssembler);
    virtual double              getWaveVelocity() ;
    double                      getVolume() const;  
    bool                        gradient(resultdata& gr) const;
    
    virtual bool                integrateDampingMatrix(assembler& as);
    virtual bool                integrateEnergy(energies& ener, const dofset::evaluation_time& when);
    virtual bool                integrateDEnergy();
    virtual bool                integrateNodalMasses();
    virtual bool                integrateMassMatrix(assembler& theAssembler);
    
    
    // virtual functions to update info from time step to time step
    virtual void                commitCurrentState();
    virtual void                resetCurrentState();
	virtual void                updateCurrentState(const dofset::evaluation_time when);
    
    
    
    smallStrainMP               *theMaterialPoint;
    
private:
    
    bool                        mixedFormulation;
    double                      volume;
    std::vector<shapefun>       theShpBar;
    std::vector<Unode*>         theNodes;
    const feliks::elementState    *elState;
    
                                solidEvalPoint();
                                solidEvalPoint(const solidEvalPoint&);
    solidEvalPoint&             operator=(const solidEvalPoint&);
};


#endif


