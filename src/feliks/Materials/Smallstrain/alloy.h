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
//  alloy.h
//  libfeliks
//
//  Created by Ignacio Romero on 5/7/13.
//  Copyright (c) 2013 Ignacio Romero. All rights reserved.
//

#ifndef __libfeliks__alloy__
#define __libfeliks__alloy__


#include "smallstrain.h"
#include "elastic.h"
#include "Math/tensor.h"

class commandLine;
class alloyMP;



class alloyMaterial : public smallStrainMaterial
{
public:
    alloyMaterial(const std::string& name,
                  const std::vector<elasticMaterial*>& mats);
	
	virtual bool            check() const;
    virtual smallStrainMP*  createMaterialPoint(const std::vector<double>* comp=0) const;
    virtual double          density() const;
	virtual void            print(std::ostream &of=std::cout) const;
    virtual void            setRandom();
    virtual bool            test();
    elasticMaterial&        operator[](const int i);
    const elasticMaterial&	operator[](const int i) const;
    

    
    
private:
    std::vector<elasticMaterial*>   theMaterials;
    
	friend class            alloyMP;
};





class alloyMP : public smallStrainMP
{
public:
    alloyMP(const alloyMaterial &m, const std::vector<double>* comp=0);
    virtual                 ~alloyMP();
    
    
    // three dimensional response
    virtual double          storedEnergy() const;
    virtual double          energyDissipationPotential() const;
    virtual void            stress(blue::istensor& sigma) const;
    virtual double          plasticSlip() const;
    
    
	virtual void            contractWithTangent(const blue::ivector &v1, const blue::ivector &v2, blue::itensor &T) const;
	virtual void            contractWithDeviatoricTangent(const blue::ivector &v1, const blue::ivector &v2, blue::itensor &T) const;
    virtual void            materialTangentTimesSymmetricTensor(const blue::istensor& M, blue::istensor& CM) const;
	virtual double          volumetricStiffness() const;
    
    
    // bookkeeping
    virtual void            updateCurrentState(const double theTime, const blue::istensor& strain, const double t);
    virtual void            commitCurrentState();
    virtual void            resetCurrentState();
    
    virtual void            setRandom(){}
    
    
private:
	const alloyMaterial&    theAlloyMaterial;
    double                  tn, tc;
    blue::istensor          En, Ec;
    double                  tempn, tempc;
    const std::vector<double>&    composition;
};




#endif /* defined(__libfeliks__alloy__) */
