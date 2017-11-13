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
 *  nodesetloading.h
 *  feliks
 *
 *  Created by Ignacio Romero on 9/9/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */


#ifndef _nodesetloading_h
#define _nodesetloading_h


#include "loading.h"
#include <string>

class commandLine;
class model;
class nodeset;


class nodesetLoading : public loading
{
	
public:
                    nodesetLoading(const commandLine& uc);
	virtual         ~nodesetLoading(){}
	
	virtual void	accumulateEnergy(energies& ener, const double time);
	virtual void	assembleLoading(assembler& asb, const double time, const double f);
	virtual void	transferLoadsToNodes(const double time){};
	
	virtual void	initialize(model& theModel);
	virtual void	print(std::ostream& of=std::cout) const;
	

private:
	std::string		theNodesetName;
	nodeset*        theNodeset;
};


#endif
