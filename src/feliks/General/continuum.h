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
//  continuum.h
//  feliks
//
//  Created by Romero Ignacio on 4/27/12.
//  Copyright (c) 2012 Universidad Politécnica de Madrid. All rights reserved.
//
//
//  tensor/vector operatios frequently arising in continuum mechanics
//
//

#pragma once
#ifndef feliks_continuum_h
#define feliks_continuum_h

#include "Math/tensor.h"


namespace feliks
{
    namespace cm
    {
        // infinitesimal kinematics
        const blue::istensor infinitesimalStrain(const blue::itensor& gradu);
        
        
        // finite strain kinematics
        const blue::istensor GreenLagrangeStrain(const blue::itensor& F);
        const blue::istensor rightCauchyGreenDeformation(const blue::itensor& F);
        const blue::istensor leftCauchyGreenDeformation(const blue::itensor& F);
        const blue::istensor logarithmicStrain(const blue::istensor& C);
        
        const blue::istensor pushforwardStress(const blue::itensor&F, const blue::istensor& S);
        const blue::istensor pullbackStrain(const blue::itensor&F, const blue::istensor& c);        
    }
}




#endif
