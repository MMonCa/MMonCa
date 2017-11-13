/*
 * LatticeFCC1.cpp
 *
 *  Created on: Aug 17, 2012
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

#include "LatticeFCC1.h"
#include "domains/Global.h"
#include "io/ParameterManager.h"

using Kernel::Coordinates;

namespace LKMC {

LatticeFCC1::LatticeFCC1(Kernel::Domain *p, const LatticeParam *lp, Kernel::M_TYPE mt)
: Lattice(p, lp, mt)
{
}

LatticeFCC1::~LatticeFCC1()
{
}

/* Ignores values for alloys */
void LatticeFCC1::fill(const Kernel::Coordinates &mc, const Kernel::Coordinates &Mc,
		std::vector<LatticeInformation> &pos, bool onlyFrame) const
{
	int m[3], M[3];
	getIndices(mc, Mc, m, M);
	double a=_latticeParameter[0];
	for(int i=m[0]; i<=M[0]; ++i)
		for(int j=m[1]; j<=M[1]; ++j)
			for(int k=m[2]; k<=M[2]; ++k)
			{
				tryInBox(pos, mc, Mc, 1, 0, toRotated(Coordinates(a*(i    ), a*(j    ), a*(k    ))));
				if(onlyFrame)
					continue;
				tryInBox(pos, mc, Mc, 1, 0, toRotated(Coordinates(a*(i+.50), a*(j+.50), a*(k    ))));
				tryInBox(pos, mc, Mc, 1, 0, toRotated(Coordinates(a*(i+.50), a*(j    ), a*(k+0.50))));
				tryInBox(pos, mc, Mc, 1, 0, toRotated(Coordinates(a*(i    ), a*(j+.50), a*(k+0.50))));
			}
}

} /* namespace LKMC */
