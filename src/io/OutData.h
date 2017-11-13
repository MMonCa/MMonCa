/*
 * OutData.h
 *
 *  Created on: Aug 6, 2012
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

#ifndef OUTDATA_H_
#define OUTDATA_H_

#include "kernel/Coordinates.h"
#include <ostream>
#include <vector>

namespace IO {


template <typename T>
struct OutDataPoint
{

	OutDataPoint() {}
	OutDataPoint(unsigned dim, const Kernel::Coordinates &c, const T&data) :
		_dim(dim), _coords(c), _data(data) {}
	unsigned _dim;
	Kernel::Coordinates _coords;
	T _data;
};

template <typename T>
struct OutDataConc
{

	OutDataConc() {}
	OutDataConc(unsigned dim, const Kernel::Coordinates &m, const Kernel::Coordinates &M, const T&data) :
		_dim(dim), _cm(m), _cM(M), _data(data) {}
	unsigned _dim;
	Kernel::Coordinates _cm, _cM;
	T _data;
};

template <typename T>
inline std::ostream & operator << (std::ostream &is, const OutDataPoint<T> & odp)
{
	if(odp._dim > 0)
		is << odp._coords._x << " ";
	if(odp._dim > 1)
		is << odp._coords._y << " ";
	if(odp._dim > 2)
		is << odp._coords._z << " ";
	is << odp._data;
	return is;
}

template <typename T>
inline std::ostream & operator << (std::ostream &is, const OutDataConc<T> & odp)
{
	if(odp._dim > 0)
		is << (odp._cm._x + odp._cM._x)/2. << " ";
	if(odp._dim > 1)
		is << (odp._cm._y + odp._cM._y)/2. << " ";
	if(odp._dim > 2)
		is << (odp._cm._z + odp._cM._z)/2. << " ";
	is << odp._data;
	return is;
}

template <typename T>
class OutDataVectorP : public std::vector<OutDataPoint<T> >
{
public:
	void push(unsigned dim, const Kernel::Coordinates &c, const T & data)
	{ std::vector<OutDataPoint<T> >::push_back(OutDataPoint<T>(dim, c, data)); }
};

template <typename T>
class OutDataVectorC : public std::vector<OutDataConc<T> >
{
public:
	void push(unsigned dim, const Kernel::Coordinates &m, const Kernel::Coordinates &M, const T & data)
	{ std::vector<OutDataConc<T> >::push_back(OutDataConc<T>(dim, m, M, data)); }
};

template <typename T>
inline std::ostream & operator << (std::ostream &is, const OutDataVectorP<T> & odv)
{
	for(typename OutDataVectorP<T>::const_iterator it=odv.begin(); it!=odv.end(); ++it)
		is << *it << std::endl;
	return is;
}

template <typename T>
inline std::ostream & operator << (std::ostream &is, const OutDataVectorC<T> & odv)
{
	for(typename OutDataVectorC<T>::const_iterator it=odv.begin(); it!=odv.end(); ++it)
		is << *it << std::endl;
	return is;
}

}
#endif /* OUTDATA_H_ */
