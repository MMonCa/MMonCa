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
 *  surfaceloading.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 9/9/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "surfaceloading.h"

#include "loading.h"
#include "Analysis/Propfactors/propfactor.h"
#include "Model/model.h"
#include "Io/usercommand.h"
#include "Io/logger.h"


surfaceLoading :: surfaceLoading(const commandLine &cl) :
	loading(cl),
	theSurfaceName(),
	theSurface(0)
{
	theType = LOAD_SURFACE;
	
	// defaults and scan
	for (int j=1; j< cl.size(); j++)
	{
		const usercommand &uc = cl[j];
		if      ( uc.keyword() == "surface" )	theSurfaceName = uc.option();

	}
	
}


void surfaceLoading :: accumulateEnergy(energies& ener, const double time)
{
    cout << "Function surfaceLoading::accumulateEnergy not implemented" << endl;
}


void surfaceLoading :: assembleLoading(assembler& asb, const double time, const double f)
{
    cout << "Function surfaceLoading::assembleLoading not implemented" << endl;
}



void surfaceLoading :: initialize(model& theModel)
{
	loading::initialize(theModel);
	theSurface = &(theModel.getFaceset(theSurfaceName));
}



void surfaceLoading::print(std::ostream& of) const
{
    loading::print(of);
    of << "\n Surface loading on surface " << theSurfaceName;
}


void surfaceLoading :: transferLoadsToNodes(const double time)
{
    cout << "Function surfaceLoading::transfer not implemented" << endl;    
}


