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
 *  festatic.h
 *  feliks
 *
 *  Created by Ignacio Romero on 20/8/05.
 *
 */


#ifndef _festatic_h
#define _festatic_h
	
#include "Analysis/analysis.h"
#include "Analysis/Integrators/integrator.h"
#include "Io/usercommand.h"
#include <string>
#include <iostream>

namespace feliks
{
    class mmoncaInterface;
};



class staticFEanalysis: public FEanalysis
{
public:
                staticFEanalysis();
                staticFEanalysis(const feliks::mmoncaInterface& mi);
	virtual     ~staticFEanalysis();
	
	bool        check();
	void        info(std::ostream &of=std::cout);


private:
	integrator *theIntegrator;
	double     finalTime;
	int        nsteps;			// steps in the incremental solution
	bool       specificSolve();
	
	friend class interactiveFEanalysis;
};

#endif
