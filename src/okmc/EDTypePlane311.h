/*
 * EDTypePlane311.h
 *
 *  Created on: Aug 14, 2013
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

#ifndef EDTYPEPLANE311_H_
#define EDTYPEPLANE311_H_

#include "EDType.h"

namespace OKMC {

class EDTypePlane311 : public EDType
{
	public:
		EDTypePlane311();

		Kernel::Coordinates surface(Kernel::RNG &, const Kernel::Coordinates *,
				const Kernel::Coordinates &center, const Kernel::ID &, const Particle *) const;
		const Particle * emitFrom (Kernel::RNG &, const Kernel::Coordinates &center, const std::vector<Particle *> &) const;
		ASPECT_TYPE getAspect() const { return EDType::aspect_plane311; }
		void 	sort311Template();
		Kernel::Coordinates getSortedCoords2Grow(unsigned size) {return _sortedTemplate311[size-1];};
	private:
		std::vector <Kernel::Coordinates >  _sortedTemplate311;
		const unsigned MAX_TEMPL_SIZE;  //only used to generate templates....
};

}

#endif /* EDTYPEPLANE311_H_ */
