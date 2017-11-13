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
 *  smallstrain.h
 *  feliks
 *
 *  Created by Ignacio Romero on 1/1/09.
 *  Copyright 2009. All rights reserved.
 *
 */


#pragma once
#ifndef _infinitesimal_h
#define _infinitesimal_h


#include "Materials/material.h"
#include "Math/tensor.h"
#include <iostream>


class commandLine;
class smallStrainMP;


class smallStrainMaterial : public material
{
	public:
                            smallStrainMaterial();
                            smallStrainMaterial(const std::map<std::string, std::pair<double, double> > &);
                            smallStrainMaterial(const commandLine& cl);

    virtual bool            check() const = 0;
    virtual smallStrainMP*	createMaterialPoint(const std::vector<double>* comp=0) const = 0;
    virtual double          density() const = 0;
	virtual void			print(std::ostream &of=std::cout) const = 0;
    virtual void            setRandom() = 0;
    
    std::map<std::string, std::pair<double, double> > eigenstrains_;
};




class smallStrainMP : public materialPoint
{
	
public:
                            smallStrainMP(const smallStrainMaterial &m);
    virtual                 ~smallStrainMP(){}

    
    virtual void			contractWithTangent(const blue::ivector &v1, const blue::ivector &v2, blue::itensor &T) const=0;
	virtual void			contractWithDeviatoricTangent(const blue::ivector &v1, const blue::ivector &v2, blue::itensor &T) const=0;
    virtual void            materialTangentTimesSymmetricTensor(const blue::istensor& M, blue::istensor& CM) const = 0;
    virtual double          plasticSlip() const = 0;
    virtual double          volumetricStiffness() const = 0;

    
    // energy
    virtual double          energyDissipationPotential() const = 0;
    virtual double          storedEnergy() const = 0;
    
    
    // stress
	virtual void            stress(blue::istensor& sigma) const =0;
    
    
    // bookkeeping
    virtual void            updateCurrentState(const double theTime, const blue::istensor& strain, const double temp)=0;
    virtual void            commitCurrentState()=0;
    virtual void            resetCurrentState()=0;
    
    
    // tests
    virtual void            setRandom()=0;
    bool                    testImplementation() const;
    
    // parent information
    double                          density() const { return theSmallStrainMaterial.density();}
    const smallStrainMaterial&      parentMaterial() const { return theSmallStrainMaterial;}

    
    
    
private:
    
	const smallStrainMaterial&	theSmallStrainMaterial;
};




#endif

