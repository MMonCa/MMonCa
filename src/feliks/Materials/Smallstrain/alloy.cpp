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
//  alloy.cpp
//  libfeliks
//
//  Created by Ignacio Romero on 5/7/13.
//  Copyright (c) 2013 Ignacio Romero. All rights reserved.
//

#include "alloy.h"
#include <cmath>
#include <vector>

#include "Io/usercommand.h"
#include "Io/logger.h"
#include "Io/message.h"
#include "Math/tensor.h"


using namespace blue;


alloyMaterial :: alloyMaterial(const std::string& name,
              const std::vector<elasticMaterial*>& mats)
:
    theMaterials(mats)
{
    eigenstrains_ = mats[0]->eigenstrains_;
}




bool alloyMaterial :: check() const
{
    bool ret = true;
    std::vector<elasticMaterial*>::const_iterator iter = theMaterials.begin();
    while (iter != theMaterials.end() )
    {
        ret = ret && (*iter)->check();
        ++iter;
    }
    
    return ret;
}




/* an object of the type "alloy" creates a material point of type alloy. The
 material point holds information that is not part of the material itself but
 that is particular of the specific (physical) point.
 */
smallStrainMP* alloyMaterial :: createMaterialPoint(const std::vector<double>* comp) const
{
	smallStrainMP* mp = new alloyMP(*this, comp);
	return mp;
}




double alloyMaterial :: density() const
{
    double rho = 0.0;
    return rho;
}






/* this function is always called once the material is defined, so apart from
 printing its information, we take the opportunity to clean up some of its
 data, in particular, setting all the possible constants
 */
void alloyMaterial :: print(std::ostream &of) const
{
    of  << "\n   Alloy material ";

    std::vector<elasticMaterial*>::const_iterator iter = theMaterials.begin();
    while (iter != theMaterials.end() )
    {
        (*iter)->print(of);
        ++iter;
    }
}




void alloyMaterial :: setRandom()
{
    std::vector<elasticMaterial*>::const_iterator iter = theMaterials.begin();
    while (iter != theMaterials.end() )
    {
        (*iter)->setRandom();
        ++iter;
    }
}




bool alloyMaterial :: test()
{
    return true;
}



elasticMaterial&	alloyMaterial::operator[](const int i)
{
    return *theMaterials[i];
}




const elasticMaterial&	alloyMaterial::operator[](const int i) const
{
    return *theMaterials[i];
}

alloyMP :: alloyMP(const alloyMaterial &m, const std::vector<double>* comp) :
smallStrainMP(m),
theAlloyMaterial(m),
tn(0.0), tc(0.0),
composition(*comp)
{
    En.setZero();
    Ec.setZero();

    for (size_t a = 0; a < comp->size(); a++)
    {
    	double kk;
    	kk = composition[a];
    	double kk2 = kk;
    	kk2++;
    }
}




alloyMP :: ~alloyMP()
{
    
}




void alloyMP :: commitCurrentState()
{
    tn = tc;
    En = Ec;
}




/* Given the fourth order tensor of alloyities C, and two vectors v, w
 * compute the second order tensor T with components
 *   T_ij = C_ipjq v_p w_q
 *
 *  which, for linear isotropic materials is just
 *  T = lambda gradV otimes gradW + mu gradW otimes gradV + mu (gradV * gradW) One
 *
 *  Note that the result is not symmetric.
 */

void alloyMP :: contractWithTangent(const ivector &v1, const ivector &v2, itensor &T) const
{
    itensor Tsum;
    Tsum.setZero();
    
    for (size_t a = 0;  a < composition.size(); a++)
    {
        itensor A, B;
	
        T = itensor::identity();
        T *= v1.dot(v2) * theAlloyMaterial[a].mu;
	
        A = itensor::dyadic(v1, v2);
        B = A.transpose();
	
        A *= theAlloyMaterial[a].lambda;
        B *= theAlloyMaterial[a].mu;
	
        T += A;
        T += B;
     
        Tsum += composition[a]*T;
    }
    
    T = Tsum;
}




void alloyMP :: contractWithDeviatoricTangent(const ivector &v1, const ivector &v2, itensor &T) const
{
    itensor Tsum;
    Tsum.setZero();
    
    for (size_t a = 0;  a < composition.size(); a++)
    {
        const double&  mu = theAlloyMaterial[a].mu;
    
        T = mu * v1.dot(v2) * itensor::identity()
        +   mu *              itensor::dyadic(v2, v1)
        -   2.0/3.0*mu *      itensor::dyadic(v1, v2);
        
        Tsum += composition[a]*T;
    }
    T = Tsum;
}



double alloyMP :: energyDissipationPotential() const
{
    return storedEnergy();
}




void alloyMP :: materialTangentTimesSymmetricTensor(const istensor& M, istensor& CM) const
{
    istensor CMsum;
    CMsum.setZero();
    
    for (size_t a = 0;  a < composition.size(); a++)
    {
        CM = theAlloyMaterial[a].lambda * M.trace() * istensor::identity()
        +  2.0* theAlloyMaterial[a].mu * M;
        
        CMsum += CM*composition[a];
    }
    
    CM = CMsum;
}




double alloyMP :: plasticSlip() const
{
    return 0.0;
}



void alloyMP :: resetCurrentState()
{
    tc = tn;
    Ec = En;
}



double alloyMP :: storedEnergy() const
{
    double e = 0.0, u=0.0;
    
    for (size_t a = 0;  a < composition.size(); a++)
    {
        const double  tr = Ec.trace();
        e = theAlloyMaterial[a].mu*Ec.squaredNorm() + 0.5*theAlloyMaterial[a].lambda*tr*tr;
        u += e*composition[a];
    }
    
    return u;
}




void alloyMP :: stress(istensor& sigma) const
{
	sigma.setZero();
    
    for (size_t a = 0;  a < composition.size(); a++)
    {
    	const elasticMaterial& mm = theAlloyMaterial[a];

    	istensor s;
    	istensor s1, s2, s3;

    	s1 = (2.0*mm.mu)*Ec;
    	s2 = (mm.lambda*Ec.trace())*istensor::identity();
    	s3 = (3.0*mm.bulk*mm.thermal*(tempc-mm.temp0)) * istensor::identity();

        s = s1 + s2 + s3;
        
        sigma += composition[a]*s;
    }    
}




void alloyMP :: updateCurrentState(const double theTime, const istensor& strain, const double t)
{
    tc    = theTime;
    Ec    = strain;
    tempc = t;
}




double alloyMP :: volumetricStiffness() const
{
    double k=0.0;
    for (size_t a = 0;  a < composition.size(); a++)
    {
        k += composition[a] * theAlloyMaterial[a].bulk;
    }
    return k;
}

