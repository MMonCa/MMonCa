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
 *  nodesetloading.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 9/9/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "nodesetloading.h"

#include "loading.h"
#include "Analysis/assembler.h"
#include "Analysis/Propfactors/propfactor.h"
#include "Model/model.h"
#include "Model/Sets/nodeset.h"
#include "Io/usercommand.h"
#include "Io/logger.h"

#include <boost/foreach.hpp>

using namespace blue;

nodesetLoading :: nodesetLoading(const commandLine &cl) :
	loading(cl),
	theNodesetName(),
	theNodeset(0)
{
	theType = LOAD_NODESET;
	
	// defaults and scan
	for (int j=1; j< cl.size(); j++)
	{
		const usercommand &uc = cl[j];
		if      ( uc.keyword() == "nodeset" )	theNodesetName = uc.option();

	}
}




void nodesetLoading :: accumulateEnergy(energies& ener, const double time)
{
	if (factor == 0) factor = &(factorcombo::getPropfactorCombinationFromLabel(getPropfactorLabel()));
    
    double Vext = 0.0;
    
    BOOST_FOREACH(node* nd, *theNodeset)
    {
        double scale = factor->eval(time, nd->coordinates(dofset::tn1));
        ivector force(0.0, 0.0, 0.0);

        if (uloaded[0]) force[0] = fU[0]*scale;
        if (uloaded[1]) force[1] = fU[1]*scale;
        if (uloaded[2]) force[2] = fU[2]*scale;
        
        Vext += force.dot(nd->displacement(dofset::tn1));
    }
    
    ener.potential -= Vext;
    ener.energy    -= Vext;    
}




void nodesetLoading :: assembleLoading(assembler& asb, const double time, const double f)
{
    if (time <= 0.0 || f <= 0.0 || theNodeset == 0) return;
    
    BOOST_FOREACH( node* nd, *theNodeset)
    {
        double scale = factor->eval( time, nd->coordinates() );

        if (uloaded[0])  asb.assemble(*nd, 0, fU[0]*scale);
        if (uloaded[1])  asb.assemble(*nd, 1, fU[1]*scale);
        if (uloaded[2])  asb.assemble(*nd, 2, fU[2]*scale);

        if (ploaded   )  asb.assemble(*nd, 3, fP   *scale);

        if (rloaded[0])  asb.assemble(*nd, 3, fR[0]*scale);
        if (rloaded[1])  asb.assemble(*nd, 4, fR[1]*scale);
        if (rloaded[2])  asb.assemble(*nd, 5, fR[2]*scale);

        if (dloaded[0])  asb.assemble(*nd, 3, fD[0]*scale);
        if (dloaded[1])  asb.assemble(*nd, 4, fD[1]*scale);
    }
}



void nodesetLoading :: initialize(model& theModel)
{
	loading::initialize(theModel);
	theNodeset = &(theModel.getNodeset(theNodesetName));
}



void nodesetLoading :: print(std::ostream& of) const
{
    of << "\n Nodeset loading on nodeset " << theNodesetName;
    loading::print(of);
}


