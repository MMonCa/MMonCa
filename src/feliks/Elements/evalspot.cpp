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
//  evalspot.cpp
//  feliks++
//
//  Created by Romero Ignacio on 12/19/11.
//  Copyright (c) 2011 Universidad Politécnica de Madrid. All rights reserved.
//

#include <iostream>
#include <stdexcept>
#include "evalspot.h"
#include "Math/feliksmath.h"
#include "Analysis/assembler.h"
#include "Io/io.h"

using namespace blue;


evalspot :: evalspot() {}

evalspot :: evalspot( const std::vector<shapefun>&  sh, const double initspacing)
:
theShp(sh),
thePositionShp(sh),
initialspacing(initspacing),
spacing(initspacing),
contactflag(false)
{
    lastconvsolution[0]=0.5/initspacing;
    lastconvsolution[1]=0.5/initspacing;
    lastconvsolution[2]=0.5/initspacing;
    Fn=itensor::identity();
}




void evalspot :: integrateResidualTangent(assembler& theAssembler)
{	
	genericIntegral(theAssembler);
	
	extern int global_iteration;
	if (global_iteration == 0 && hasConstrainedDofs())
	{
		matrix dummytangent(0,0);
		BCassembler bca(theAssembler);
		genericIntegral(bca);		
	}
}

 const ivector   evalspot ::    getTauVelocity(const dofset::evaluation_time when)
{
    throw runtime_error("This evalspot type has no getTauVelocity function programmed.");
}



void  evalspot :: readapt_nodes()
{
    throw runtime_error("This evalspot type has no readapt_nodes function programmed.");
}

void  evalspot :: readapt_PositionNodes()
{
    throw runtime_error("This evalspot type has no readapt_PositionNodes function programmed.");
}




void evalspot:: TurnOnContact()
{
    contactflag = true;
}

void evalspot:: TurnOffContact()
{
    contactflag = false;
}
