/* Finite Element Method Module
 *
 * Author: ignacio.romero@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain,
 *      and       Technical University of Madrid (UPM), Madrid, Spain
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
/* feliksutil.h
 *
 * iro, april 2009 
 *
 */

#ifndef _feliksutil_h
#define _feliksutil_h


#include <vector>
#include <utility>
#include <stdint.h>
#include <cassert>
#include <algorithm>


#ifdef WITHTBB
#include "tbb/parallel_for_each.h"
#endif


// memory utilities
inline bool isAligned(void* data, int alignment)
{
	// check that the alignment is a power of two
	assert((alignment & (alignment-1)) == 0); 
	return ((uintptr_t)data & (alignment-1)) == 0;
}

// use DATA_ALIGN to align data in memory for static variables
// the usage is: DATA_ALIGN(int buffer[16], 128)
//
#if defined(VISUAL_STUDIO)
#define DATA_ALIGN(declaration, alignment) __declspec(align(alignment)) declaration
#elseif defined(GCC)
#define DATA_ALIGN(declaration, alignment) declaration __attribute__ ((aligned (alignment)))
#else
#define DATA_ALIGN(declaration, alignment)
#endif



// cartesian product of two lists/std::vectors...
// usage: cartesianProduct< double, double >
//
// example
// std::vector<int> v1, v2;
// v1.push_back(11); v1.push_back(12); v2.push_back(21); v2.push_back(22);
// cartesianProduct<int,int> v( v1, v2);
//
// std::vector< pair<int,int> >::iterator iter = v.begin();
// while ( iter != v.end() )
// {
//		cout << "\n (" << iter->first << " , " << iter->second << ")" << flush;
//		++iter;
//	}
	

template<class T1, class T2> 
class cartesianProduct
	{
	public:
		cartesianProduct(std::vector<T1> cont1, std::vector<T2> cont2)
		{
			typename std::vector<T1>::iterator iter1 = cont1.begin();
			while ( iter1 != cont1.end() )
			{
				typename std::vector<T2>::iterator iter2 = cont2.begin();
				while ( iter2 != cont2.end() )
				{
					data.push_back( std::pair<T1,T2>( *iter1, *iter2 ) );
					iter2++;
				}
				iter1++;
			}
		}
		
		typename std::vector< std::pair<T1,T2> >::iterator begin(){return data.begin();}
		typename std::vector< std::pair<T1,T2> >::iterator end() {return data.end();}
		
		
	private:
		std::vector< std::pair<T1,T2> > data;
	};



template<class InputIterator, class Function>
void feliks_for_each(InputIterator first, InputIterator last, Function f)
{
#ifdef WITHTBB
    tbb::parallel_for_each(first, last, f);
#else
    std::for_each(first, last, f);
#endif
}





#endif
