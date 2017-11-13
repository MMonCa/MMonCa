/*
 * ParticleType.h
 *
 *  Created on: Feb 24, 2011
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

#ifndef PARTICLETYPE_H_
#define PARTICLETYPE_H_

#include <map>
#include <ostream>
#include "Material.h"

namespace Kernel
{
	const unsigned char MAX_IMPURITIES = 9;
	const unsigned char MAX_PARTICLES = MAX_IMPURITIES*5                                                                                                                                                       ;
	const unsigned MAX_STATES = 7;
	const unsigned MID_STATE = 3;
	const unsigned char UNDEFINED_TYPE = MAX_PARTICLES + 1;
	const unsigned UNDEFINED_STATE = MAX_STATES + 1;
	const unsigned char IV_I = 0, IV_V = 1, IV_NONE = 2;
	typedef unsigned char P_TYPE;
	const P_TYPE V_TYPE = 0;
	enum P_POS { POS_0=0, POS_1=1, POS_I, POS_V, NO_POS, UNDEFINED_POS };
	struct ID
	{
		std::map<P_TYPE, unsigned> _pt;
		std::map<P_POS, unsigned> _pos;
		M_TYPE _mt;
		bool operator<(const ID &) const;
		bool operator==(const ID &) const;
	};
	std::ostream & operator<<(std::ostream &, const ID &id);
}

#endif /* PARTICLETYPE_H_ */
