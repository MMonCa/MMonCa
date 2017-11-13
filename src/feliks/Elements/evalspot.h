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
//  evalspot.h
//  feliks++
//
//  Created by Romero Ignacio on 12/8/11.
//  Copyright (c) 2011 Universidad Politécnica de Madrid. All rights reserved.
//

#pragma once
#ifndef feliks_element__evalspot_h
#define feliks_element__evalspot_h


#include "Elements/Interpolation/interpolation.h"

#include <vector>
#include "Math/feliksmath.h"

class assembler;




class evalspot
{
    
public:
    virtual                     ~evalspot(){}
    
    virtual bool                check() const=0;
    virtual blue::ivector             coordinates( const dofset::evaluation_time when=dofset::tna) const=0;
    virtual blue::ivector             getReferenceCoordinates() const=0;
    virtual void                getDofNumbers(std::vector<int>& theDofs) const=0;
    virtual bool                gradient(resultdata& gr) const =0;
    virtual bool                genericIntegral(assembler& as)=0;
    
    void                        TurnOnContact();
    void                        TurnOffContact();
    bool                        HasContact(){return contactflag;}
    virtual double              mass() const = 0;
    
    // virtual functions to update info from time step to time step
    virtual void                commitCurrentState()=0;
    virtual void                resetCurrentState()=0;
	virtual void                updateCurrentState(const dofset::evaluation_time when)=0;
    
    
    virtual bool                integrateDampingMatrix(assembler& as)=0;
    virtual bool                integrateEnergy(energies&, const dofset::evaluation_time& when)=0;
    virtual bool                integrateDEnergy()=0;
    virtual bool                integrateMassMatrix(assembler& as)=0;
    virtual bool                integrateNodalMasses()=0;
	void                        integrateResidualTangent(assembler& as);
    
    virtual void                readapt_nodes();
    virtual void                readapt_PositionNodes();

    
    std::vector<shapefun>&      getShapefunctions(){return theShp;}
    std::vector<shapefun>&      getPositionShapefunctions(){return thePositionShp;}
    const double&               getInitialSpacing(){return initialspacing;}
    blue::ivector&                    getLastConvSolution(){return lastconvsolution;}
    const double&               getSpacing(){return spacing;}
    blue::itensor&                    getDefGradientTn(){return Fn;}
    virtual const blue::ivector       getTauVelocity(const dofset::evaluation_time when);
    virtual void                updatestabilizedmass(){};

    virtual double              getWaveVelocity()  =0;
    virtual bool                hasConstrainedDofs()=0;
    void                        resetSpacing(){spacing = initialspacing;}
    void                        updateSpacing(const double& h){spacing = h;}
    void                        setLastConvSolution(blue::ivector& a){lastconvsolution = a;}

    
protected:
    evalspot();
    evalspot(const std::vector<shapefun>&   sh, const double inispacing);
    
    
private:
    std::vector<shapefun>       theShp;
    std::vector<shapefun>       thePositionShp;
    double                      initialspacing;
    double                      spacing;
    blue::itensor                     Fn;
    blue::ivector                     lastconvsolution;
    bool                        contactflag;
    
    friend class                fsolidEvalPoint;
    friend class                solidEvalPoint;
    friend class                fluidEvalPoint;
    friend class                navierstokesEvalPoint;
};








class shellSection;
class beamSection;


class evalFiber : public evalspot
{
    
private:
    shellSection     *theShellSection;
};



class evalsection : public evalspot
{
    
private:
    beamSection     *theBeamSection;
};


#endif
