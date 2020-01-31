/*
 * DefectIterator.h
 *
 *  Created on: Jun 20, 2012
 *      Author: ignacio.martin@imdea.org
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

#ifndef DEFECTDATAITERATOR_H_
#define DEFECTDATAITERATOR_H_

#include <memory>

namespace OKMC { class Defect; }
namespace Kernel { class Event;
					template <typename T_EVENT, typename T_RATE> class SelectorIterator; }

namespace Domains {

class Global;

class DefectIterator {
public:
	DefectIterator(const Kernel::SelectorIterator<Kernel::Event, double> &);
	DefectIterator(const DefectIterator &);
	DefectIterator & operator=(const DefectIterator &);
	~DefectIterator();

	DefectIterator & operator++();
	bool             operator!=(const DefectIterator &m) const;
	bool             operator==(const DefectIterator &m) const { return !operator!=(m); }

	const OKMC::Defect * operator*()  const;
	const OKMC::Defect * operator->() const;

private:
	std::shared_ptr<Kernel::SelectorIterator<Kernel::Event, double> > _pSi;

	unsigned _nDomain;
	unsigned _nSubDomain;
	unsigned _nLevel;

	friend class Global;
};

} /* namespace Domains */
#endif /* DEFECTDATAITERATOR_H_ */
