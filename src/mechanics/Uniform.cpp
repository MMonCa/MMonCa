/*
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

#include "Uniform.h"

#include "io/Diagnostic.h"
#include "io/FileParameters.h"
#include "io/ParameterManager.h"

#include "lkmc/Lattice.h"

#include "kernel/Domain.h"
#include "kernel/Mesh.h"

namespace Mechanics {

Uniform::Uniform(Kernel::Domain *p) : MechInterface(p)
{
	_C = boost::numeric::ublas::zero_matrix<double>(6,6);
	_S = boost::numeric::ublas::zero_matrix<double>(6,6);

	// WARNING: Silicon parameters
	double C11 = Domains::global()->getFileParameters()->getFloat("Silicon/Mechanics/C11"); // Pa
	double C12 = Domains::global()->getFileParameters()->getFloat("Silicon/Mechanics/C12"); // Pa
	double C44 = Domains::global()->getFileParameters()->getFloat("Silicon/Mechanics/C44"); // Pa

	double S11 = Domains::global()->getFileParameters()->getFloat("Silicon/Mechanics/S11");
	double S12 = Domains::global()->getFileParameters()->getFloat("Silicon/Mechanics/S12");
	double S44 = Domains::global()->getFileParameters()->getFloat("Silicon/Mechanics/S44");

	_C(0,0) = C11;
	_C(1,1) = C11;
	_C(2,2) = C11;
	_C(3,3) = C44;
	_C(4,4) = C44;
	_C(5,5) = C44;
	_C(1,2) = C12;
	_C(2,1) = C12;
	_C(1,3) = C12;
	_C(3,1) = C12;
	_C(2,3) = C12;
	_C(3,2) = C12;

	_S(0,0) = S11;
	_S(1,1) = S11;
	_S(2,2) = S11;
	_S(3,3) = S44;
	_S(4,4) = S44;
	_S(5,5) = S44;
	_S(0,1) = S12;
	_S(0,2) = S12;
	_S(1,0) = S12;
	_S(1,2) = S12;
	_S(2,0) = S12;
	_S(2,1) = S12;
}

void Uniform::import() const
{
	ublas::vector<double> stress = ublas::zero_vector<double>(6);
	//ublas::vector<double> strain = ublas::zero_vector<double>(6);

	stress(0) = Domains::global()->getFileParameters()->getFloat("Mechanics/Uniform/stress.xx");
	stress(1) = Domains::global()->getFileParameters()->getFloat("Mechanics/Uniform/stress.yy");
	stress(2) = Domains::global()->getFileParameters()->getFloat("Mechanics/Uniform/stress.zz");

	LOWMSG("Assigning to elements...");

	if(Domains::global()->getDomains() != 1)
			ERRORMSG("Sorry, Uniform interface not supported for " << Domains::global()->getDomains() << " domains");
	Kernel::Domain *pDomain = Domains::global()->getDomain(0);

	for(Kernel::Mesh::iterator mit = pDomain->_pMesh->begin();	mit != pDomain->_pMesh->end(); ++mit)
		if(Domains::global()->PM()->getMaterialName(mit->getMaterial()) != "Gas")
		{
			boost::numeric::ublas::matrix<double> m(3, 3);

			pDomain->_pLat[mit->getMaterial()]->getInvR(m);
			ublas::matrix<double> t(6, 6); // tensor conversion

			for (unsigned i = 0; i < 3; ++i) {
				t(i, 0) = m(i, 0) * m(i, 0);
				t(i, 1) = m(i, 1) * m(i, 1);
				t(i, 2) = m(i, 2) * m(i, 2);
				t(i, 3) = 2. * m(i, 1) * m(i, 2);
				t(i, 4) = 2. * m(i, 0) * m(i, 2);
				t(i, 5) = 2. * m(i, 0) * m(i, 1);
			}

			t(3, 0) = m(1, 0) * m(2, 0);
			t(3, 1) = m(1, 1) * m(2, 1);
			t(3, 2) = m(1, 2) * m(2, 2);
			t(3, 3) = m(1, 1) * m(2, 2) + m(1, 2) * m(2, 1);
			t(3, 4) = m(1, 0) * m(2, 2) + m(1, 2) * m(2, 0);
			t(3, 5) = m(1, 0) * m(2, 1) + m(1, 1) * m(2, 0);

			t(4, 0) = m(2, 0) * m(0, 0);
			t(4, 1) = m(2, 1) * m(0, 1);
			t(4, 2) = m(2, 2) * m(0, 2);
			t(4, 3) = m(0, 1) * m(2, 2) + m(0, 2) * m(2, 1);
			t(4, 4) = m(0, 0) * m(2, 2) + m(0, 2) * m(2, 0);
			t(4, 5) = m(0, 0) * m(2, 1) + m(0, 1) * m(2, 0);

			t(5, 0) = m(0, 0) * m(1, 0);
			t(5, 1) = m(0, 1) * m(1, 1);
			t(5, 2) = m(0, 2) * m(1, 2);
			t(5, 3) = m(0, 1) * m(1, 2) + m(0, 2) * m(1, 1);
			t(5, 4) = m(0, 0) * m(1, 2) + m(0, 2) * m(1, 0);
			t(5, 5) = m(0, 0) * m(1, 1) + m(0, 1) * m(1, 0);

			mit->stress() = stress;
			mit->strain() = ublas::prod(_S, mit->stress());
		}
	LOWMSG("Done");
}

}
