/*
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

#ifndef UNIFORM_H_
#define UNIFORM_H_

#include "MechInterface.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace Mechanics {

class Uniform : public MechInterface {
public:
	Uniform(Kernel::Domain *);
	virtual ~Uniform() {}

	virtual void import() const;

private:
	// boost::numeric::ublas::matrix<double> _C_cryst;
	boost::numeric::ublas::matrix<double> _C;
	boost::numeric::ublas::matrix<double> _S;
};

} /* namespace Mechanics */

#endif /* UNIFORM_H_ */
