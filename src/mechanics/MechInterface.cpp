/*
 * MechInterface.cpp
 *
 *  Created on: Aug 18, 2011
 *
 * Author: ignacio.martin@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain
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

#include "NoMechanics.h"
#include "Uniform.h"
#include "io/Diagnostic.h"
#include "io/FileParameters.h"
#include "kernel/Domain.h"

using namespace Mechanics;

MechInterface::MechInterface(Kernel::Domain *p) : _pDomain(p)
{
	_dim   = Domains::global()->getFileParameters()->getInt("Mechanics/General/dimension");
	_every = Domains::global()->getFileParameters()->getFloat("Mechanics/General/every");
}

MechInterface * MechInterface::get(Kernel::Domain *p)
{
	std::string model = Domains::global()->mechModel();
	if(model == "Uniform")
		return new Uniform(p);
	else if(model == "None")
		return new NoMechanics(p);
	else
		ERRORMSG("Mechanics interface " << model << " not accepted");
	return 0;
}
