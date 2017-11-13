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
/*
 *  cell.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 10/3/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "cell.h"
#include "chain.h"
#include "cellcomplex.h"
#include <map>
#include <set>
#include "boost/foreach.hpp"

using namespace feliks::topology;


#ifdef DEBUG
int cell::_last_label[4] = {0,0,0,0};
#endif


cell :: cell(const int n) :
	_dimension(n),
	_boundary(n-1),
	_coboundary(n+1)
{
#ifdef DEBUG
	_label = cell::_last_label[n]++;
#endif
}

cell :: ~cell()
{
}



// checks if cell *this belongs to the cells of dimension this->dimension
// contained in mother cell
bool cell :: belongsToClosedCell(const cell& mother) const
{
    int mdim = mother.dimension();
    int cdim = this->dimension();
    
    if ( mdim < cdim ) return false;
    
    std::map< int , std::set<const cell*> > m;
    
    std::set<const cell *> s; 
    s.insert( &mother );
    m[ mdim ] = s;
    
    for (int p= mdim; p>cdim; p--)
    {
        std::set<const cell *> sp; 
        BOOST_FOREACH(const cell* c, m[p])
        {
            chain::const_iterator iter = c->boundary().begin();
            while ( iter != c->boundary().end() )
            {
                sp.insert( (*iter).first );
                ++iter;
            }
            m[p-1] = sp;
        }
    }
    
    return m[cdim].find( this ) != m[cdim].end() ;
}


chain& cell :: boundary()
{
	return _boundary;
}


const chain& cell :: boundary() const
{
	return _boundary;
}


cochain& cell :: coboundary()
{
	return _coboundary;
}


const cochain& cell :: coboundary() const
{
	return _coboundary;	
}




void cell :: boundary2coboundary()
{
	std::map<cell *, int>::iterator iter = _boundary.begin();
	while ( iter != _boundary.end() )
	{
		(*iter).first->coboundary()[this] = (*iter).second;
		++iter;
	}	
}



void cell :: coboundary2boundary()
{	
	std::map<cell *, int>::iterator iter = _coboundary.begin();
	while ( iter != _coboundary.end() )
	{
		(*iter).first->boundary()[this] = (*iter).second;
		++iter;
	}
}


const int cell :: dimension() const
{
	return _dimension;
}


void cell :: print(std::ostream& of) const
{
	of << _dimension << "-cell ";

#ifdef DEBUG
    of << _label;
#endif
}


// touching cells of dimension d
std::set<cell* >  cell::touchingCells(const int d) const
{    
    // two sets of cells to extract coboundaries
    std::set<cell *> low;
    std::set<cell *> high;
    
        
    // fill high with cells in coboundary
    {
        const cochain& co = this->coboundary();
        std::map<cell *,int>::const_iterator iter = co.begin();
        while ( iter != co.end() )
        {
            high.insert( (*iter).first );
            ++iter;
        }
    }
    
    
    
    // sequence of coboundaries
    for (size_t a = this->dimension()+2; a<=d; a++)
    {
        low = high;
        high.clear();
        
        std::set<cell *>::iterator citer = low.begin();
        while (citer != low.end())
        {
            // put the coboundary of citer in high
            {
                const cochain& co = (*citer)->coboundary();
                std::map<cell *,int>::const_iterator iter = co.begin();
                while ( iter != co.end() )
                {
                    high.insert( (*iter).first );
                    ++iter;
                }
            }
            
            ++citer;
        }
    }
        
            
    return high;
}




