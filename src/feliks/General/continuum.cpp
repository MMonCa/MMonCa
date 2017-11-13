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
//  continuum.cpp
//  feliks
//
//  Created by Romero Ignacio on 4/27/12.
//  Copyright (c) 2012 Universidad Politécnica de Madrid. All rights reserved.
//

#include <iostream>

#include "continuum.h"
#include "Math/tensor.h"

using namespace blue;


namespace feliks
{
    namespace cm
    {
        
        
        const istensor GreenLagrangeStrain(const itensor& F)
        {
            istensor C = rightCauchyGreenDeformation(F);
            
            return 0.5*( C - istensor::identity() );
        }
        
        
        
        const istensor rightCauchyGreenDeformation(const itensor& F)
        {
            return  istensor::tensorTransposedTimesTensor(F);
        }
        
        
        
        const istensor infinitesimalStrain(const itensor& gradu)
        {
            return istensor::symmetricPartOf(gradu);
        }
        
        
        
        const istensor leftCauchyGreenDeformation(const itensor& F)
        {
            istensor B = istensor::tensorTimesTensorTransposed(F);
            return B;
        }
        
        
        
        const istensor logarithmicStrain(const istensor& C)
        {
            ivector lambda2;
            ivector N[3];            
            C.spectralDecomposition(N, lambda2);
            
            istensor E;
            E.setZero();
            for (int a=0; a<3; a++)
                E.addScaledVdyadicV(0.5*log(lambda2(a)), N[a]);
            
            return E;
        }
        
        
        
        const istensor pushforwardStress(const itensor&F, const istensor& S)
        {
            return istensor::FSFt(F, S);
        }
        
        
        
        const istensor pullbackStrain(const itensor&F, const istensor& c)
        {
            return istensor::FtCF(F, c);
        }
        
        
    }
};

