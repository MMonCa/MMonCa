/*
 * MeshElementIterator.h
 *
 *  Created on: Jun 20, 2012
 *      Author: ignacio.martin@imdea.org
 *
 *      Iterator to obtain mesh data information
 *      It is "pseudo", because it does not contain a pointer to a value
 *      but it is, somehow, the value itself.
 *
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

#ifndef MESHDATAITERATOR_H_
#define MESHDATAITERATOR_H_

namespace Kernel { class MeshElement; }

namespace Domains {

class Global;

class MeshElementIterator {
public:
	MeshElementIterator() : _pME(0)  {}
	MeshElementIterator(const MeshElementIterator &m);
	MeshElementIterator & operator=(const MeshElementIterator &m);

	MeshElementIterator & operator++();
	bool                  operator!=(const MeshElementIterator &m) const { return m._pME != _pME; }
	bool                  operator==(const MeshElementIterator &m) const { return !operator!=(m); }

	const Kernel::MeshElement * operator*()  const { return _pME; }
	const Kernel::MeshElement * operator->() const { return _pME; }
	Kernel::MeshElement       * modify() { return _pME; }

private:
	Kernel::MeshElement *_pME;
	unsigned _nElement;
	unsigned _nDomain;

	friend class Domains::Global;
};

} /* namespace Domains */

#endif /* MESHDATAITERATOR_H_ */
