/*
 * IDContainer.h
 *
 *  Created on: Oct 21, 2014
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
 *
 *
 *  This class is a container, where the "specifier"
 *  is Kernel::ID
 *
 *  It uses hashes internally.
 *
 */

#ifndef IDContainer_H_
#define IDContainer_H_

#include "ParticleType.h"
#include <vector>
#include <map>
#include <cassert>

namespace Kernel {

template <typename T_CONTAINER>
class IDContainer {
public:
	unsigned cluster2hash(const ID &) const;
	const ID & hash2cluster(unsigned) const;
	unsigned addHash(const ID &);
	unsigned invalidHash() const { return _inverseHash.size(); }
	std::vector<T_CONTAINER> _map; //container using the hashes

private:
	//   size              map      hash                  bigger impurity
	// example: I23C2P4, impurity is I, size is 23
	std::vector< std::map<ID,unsigned> > _hashTable[MAX_IMPURITIES];
	std::vector<ID> _inverseHash;
};

template <typename T_CONTAINER>
unsigned IDContainer<T_CONTAINER>::addHash(const ID &id) {
	std::map<P_TYPE, unsigned>::const_iterator it = id._pt.begin();
	unsigned size = it->second;
	P_TYPE pt = it->first;
	++it;

	for (; it != id._pt.end(); ++it)
		if (it->second > size) {
			size = it->second;
			pt = it->first;
		}
	//OK, I have size, pt and mt now.
	while (_hashTable[pt].size() <= size)
		_hashTable[pt].push_back(std::map<ID, unsigned>());
	std::map<ID, unsigned> &theList = _hashTable[pt][size];
	if (theList.find(id) == theList.end()) {
		theList[id] = _inverseHash.size();
		_map.push_back(T_CONTAINER());
		_inverseHash.push_back(id);
	}
	return theList[id];
}

template <typename T_CONTAINER>
unsigned IDContainer<T_CONTAINER>::cluster2hash(const ID &id) const {
	std::map<P_TYPE, unsigned>::const_iterator it = id._pt.begin();
	assert(it != id._pt.end());
	unsigned size = it->second;
	P_TYPE pt = it->first;
	++it;

	for (; it != id._pt.end(); ++it)
		if (it->second > size) {
			size = it->second;
			pt = it->first;
		}
	//OK, I have size, pt and mt now.
	if (_hashTable[pt].size() <= size)
		return _inverseHash.size(); //nohash

	const std::map<ID, unsigned> &theList = _hashTable[pt][size];
	std::map<ID, unsigned>::const_iterator it2 = theList.find(id);
	if (it2 != theList.end())
		return it2->second;
	return _inverseHash.size();
}

template <typename T_CONTAINER>
const ID & IDContainer<T_CONTAINER>::hash2cluster(unsigned hash) const {
	static ID theMap;
	if (hash < _inverseHash.size())
		return _inverseHash[hash];
	return theMap;
}

} /* namespace Kernel */

#endif /* IDHASH_H_ */
