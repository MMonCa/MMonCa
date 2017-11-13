/*
 * ChargedStates.h
 *
 *  Created on: Aug 23, 2013
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

#ifndef CHARGEDSTATES_H_
#define CHARGEDSTATES_H_

#include "kernel/StateManager.h"

namespace Kernel { class MeshElement; }

namespace OKMC {

class ChargedStates: public Kernel::StateManager
{
public:
	ChargedStates(Kernel::Domain *);
	virtual ~ChargedStates();

	virtual unsigned changeState      (Kernel::SubDomain *, Kernel::P_TYPE, const Kernel::MeshElement *pME) const;
	virtual bool     checkStateBarrier(Kernel::SubDomain *, Kernel::P_TYPE, unsigned, const Kernel::MeshElement *, const Kernel::MeshElement *) const;
	virtual double   bindingShift     (Kernel::P_TYPE, Kernel::P_POS, unsigned, const Kernel::MeshElement *) const;
	virtual bool     canInteract      (Kernel::P_TYPE, unsigned, Kernel::P_TYPE, unsigned, const Kernel::MeshElement *) const;
	virtual unsigned interactionState (Kernel::P_TYPE, unsigned, Kernel::P_TYPE, unsigned, const Kernel::MeshElement *, Kernel::P_TYPE) const;
	virtual unsigned breakUpState     (Kernel::P_TYPE, Kernel::P_POS, unsigned, Kernel::M_TYPE) const;
};

} /* namespace OKMC */
#endif /* CHARGESTATES_H_ */
