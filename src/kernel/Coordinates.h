/* Author: ignacio.martin@imdea.org
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

#ifndef COORDINATES_H
#define COORDINATES_H
#include <cassert>
#include <istream>
#include <ostream>
#include <iomanip>
#include <cmath>

namespace Kernel
{
	template <typename T>
	struct CoordinatesT
	{
		using type = T;
		
		CoordinatesT() {}
		CoordinatesT(T x, T y, T z) : _x(x), _y(y), _z(z) {}
		T _x, _y, _z;
		T operator[](int i) const { switch(i)
		{ case 0: return _x; case 1: return _y; case 2: return _z; default: assert(false); } return _x; }
		T & operator[](int i) { switch(i)
		{ case 0: return _x; case 1: return _y; case 2: return _z; default: assert(false); } return _x; }
		CoordinatesT<T> & operator+=(T f) { _x+=f; _y+=f; _z+=f; return *this; }
		CoordinatesT<T> & operator-=(T f) { _x-=f; _y-=f; _z-=f; return *this; }
		bool isInto(const CoordinatesT<T> &m, const CoordinatesT<T> &M) const
		{ return ( _x >= m._x && _x < M._x && _y >= m._y && _y < M._y && _z >= m._z && _z < M._z); }
		bool isIntoEqual(const CoordinatesT<T> &m, const CoordinatesT<T> &M) const
		{ return ( _x >= m._x && _x <= M._x && _y >= m._y && _y <= M._y && _z >= m._z && _z <= M._z); }
		bool operator==(const CoordinatesT<T> &c) const { return c._x == _x && c._y == _y && c._z == _z; }
		bool operator!=(const CoordinatesT<T> &c) const { return !(c == *this); }

		bool operator<(const CoordinatesT<T> &c) const {
			if (_x < c._x) return true;
			if (_x > c._x) return false;
			if (_y < c._y) return true;
			if (_y > c._y) return false;
			if (_z < c._z) return true;
			if (_z > c._z) return false;
			return false;
		}


		void operator+=(const CoordinatesT &c) { _x += c._x; _y += c._y; _z += c._z; }
		void operator-=(const CoordinatesT &c) { _x -= c._x; _y -= c._y; _z -= c._z; }
		void operator*=(T f) { _x*=f; _y*=f; _z*=f; }
		void operator/=(T f) { _x/=f; _y/=f; _z/=f; }
		CoordinatesT operator+(const CoordinatesT &c) const { return CoordinatesT(_x+c._x, _y+c._y, _z+c._z); }
		CoordinatesT operator-(const CoordinatesT &c) const { return CoordinatesT(_x-c._x, _y-c._y, _z-c._z); }
		CoordinatesT operator*(T f) const { return CoordinatesT(_x*f, _y*f, _z*f); }
		CoordinatesT operator/(T f) const { return CoordinatesT(_x/f, _y/f, _z/f); }
		T            operator*(const CoordinatesT &c) const { return _x*c._x + _y*c._y + _z*c._z; }
		T abs() const { return std::sqrt(*this * *this); }
		CoordinatesT product(const CoordinatesT &c) const { return CoordinatesT(_y*c._z - _z*c._y, _z*c._x - _x*c._z, _x*c._y - _y*c._x); }

		static bool intersection(const CoordinatesT &m1, const CoordinatesT &M1, CoordinatesT &m2, CoordinatesT &M2);
	};
	template <typename T>
	inline std::ostream & operator << (std::ostream &os, const CoordinatesT<T> & c)
	{
		os << "(" <<
				std::setprecision(30) << c._x << ", " <<
				std::setprecision(30) << c._y << ", " <<
				std::setprecision(30) << c._z << ")";
		return os;
	}
	template <typename T>
	inline std::istream & operator >> (std::istream &is, CoordinatesT<T> & c)
	{
		char ch;
		is >> ch >> c._x >> ch >> c._y >> ch >> c._z >> ch;
		return is;
	}
	//sets the intersection in m2 and M2
	template <typename T>
	inline bool CoordinatesT<T>::intersection(const CoordinatesT &m1, const CoordinatesT &M1, CoordinatesT &m2, CoordinatesT &M2)
	{
		for(unsigned i=0; i<3; ++i)
		{
			m2[i] = std::max(m1[i], m2[i]);
			M2[i] = std::min(M1[i], M2[i]);
			if(M2[i] < m2[i])
				return false;
		}
		return true;
	}
	typedef CoordinatesT<float> Coordinates;

	template<typename T>
	struct VtkOrderLess final {
		bool operator()(CoordinatesT<T> const& aLhs, CoordinatesT<T> const& aRhs) {
			if (aLhs._z < aRhs._z) return true;
			if (aLhs._z > aRhs._z) return false;
			if (aLhs._y < aRhs._y) return true;
			if (aLhs._y > aRhs._y) return false;
			if (aLhs._x < aRhs._x) return true;
			if (aLhs._x > aRhs._x) return false;
			return false;
		}
	};
}
	 
#endif
