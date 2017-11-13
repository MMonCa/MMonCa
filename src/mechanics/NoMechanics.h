/*
 * NoMechanics.h
 *
 *  Created on: Aug 21, 2012
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

#ifndef NOMECHANICS_H_
#define NOMECHANICS_H_

#include "MechInterface.h"

namespace Mechanics {

class NoMechanics : public MechInterface {
public:
	NoMechanics(Kernel::Domain *p) : MechInterface(p) {}
	virtual ~NoMechanics() {}

	virtual void import() const {}
};

} /* namespace Mechanics */
#endif /* NOMECHANICS_H_ */
