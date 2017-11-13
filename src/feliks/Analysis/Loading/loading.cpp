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
 *  loading.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 4/5/10.
 *  Copyright 2010 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "loading.h"
#include "Analysis/Loading/surfaceloading.h"
#include "Analysis/Propfactors/propfactor.h"
#include "Model/model.h"
#include "Io/usercommand.h"
#include "Io/logger.h"

loading :: loading() :
	factor(0),
	factorlabel(0)
{
	fU.setZero();
	fR.setZero();
	fP = 0.0;
	fD.setZero();
}


loading :: loading(const int flabel) :
	factor(0),
	factorlabel(flabel)
{
	fU.setZero();
	fR.setZero();
	fP = 0.0;
	fD.setZero();
}


loading :: loading(const commandLine& cl) :
	factor(0),
	factorlabel(0)
{
	fU.setZero();
	fR.setZero();
	fD.setZero();
	fP = 0.0;
	
	uloaded[0] = uloaded[1] = uloaded[2] = false;
	rloaded[0] = rloaded[1] = rloaded[2] = false;
	dloaded[0] = dloaded[1] = false;
	ploaded    = false;
	
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		
        if      (uc.keyword() == "scaling")  factorlabel = uc.value();

		else if	(uc.keyword() ==  "fx")		fU[0]   = uc.value(), uloaded[0] = true;
		else if	(uc.keyword() ==  "fy")		fU[1]   = uc.value(), uloaded[1] = true;
		else if	(uc.keyword() ==  "fz")		fU[2]   = uc.value(), uloaded[2] = true;
		
		else if	(uc.keyword() ==  "fp")		fP      = uc.value(), ploaded    = true;
		
		else if	(uc.keyword() ==  "frotx")	fR[0]   = uc.value(), rloaded[0] = true;
		else if	(uc.keyword() ==  "froty")	fR[1]   = uc.value(), rloaded[1] = true;
		else if	(uc.keyword() ==  "frotz")	fR[2]   = uc.value(), rloaded[2] = true;
		
		else if	(uc.keyword() ==  "fdirx")	fD[0]   = uc.value(), dloaded[0] = true;
		else if	(uc.keyword() ==  "fdiry")	fD[1]   = uc.value(), dloaded[1] = true;
	}
}



factorcombo&  loading ::  getScalingFactor()
{
	if (factor == 0) factor = &(factorcombo::getPropfactorCombinationFromLabel(factorlabel));
	return *factor;
}


void loading :: initialize(model& theModel)
{
	factor = &(factorcombo::getPropfactorCombinationFromLabel(getPropfactorLabel()));
}



void loading :: print(std::ostream& of) const
{
    if (uloaded[0]) of << "\n       fx           : " << fU[0];
    if (uloaded[1]) of << "\n       fy           : " << fU[1];
    if (uloaded[2]) of << "\n       fz           : " << fU[2];
    of <<                 "\n       Scaling      : " << factorlabel;
}



void loading :: scan(const commandLine &cl, model& theModel)
{
	// Check that there are enough usercommands in the list. At least there must be 2
    if (cl.size() < 2)
    {
		logger :: mainlog << "\n Warning: loading cannot be defined. Incomplete definition.";
		logger :: mainlog << "\n Input is :";
        cl.print();
    }
	
	// get the second usercommand which must be of the form 'type = ...'
    const usercommand &uc  = cl[1];
    if (uc.keyword() == "surface")	
	{
		loading *ld = new surfaceLoading(cl);
		theModel.add(*ld);
	}

	else  
		logger::warnings << "No loading of type " << uc.option() << " is programmed." << endl;
	
}


