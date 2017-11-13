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
 *  cellcomplex.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 10/3/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "cellcomplex.h"
#include "cell.h"
#include "chain.h"

#include <algorithm>
#include <vector>
#include <set>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "boost/foreach.hpp"

using namespace std;
using namespace feliks::topology;


static int ccrandom(int lower, int upper);


cellcomplex::cellcomplex(const int dim) :
dimension(dim)
{
	resize(dim+1);
    cellTable.resize(dim+1);
}



cellcomplex::~cellcomplex()
{
	std::vector< set < cell *> >::iterator iter = this->begin();
	while ( iter != this->end() )
	{
		std::set< cell *>::iterator siter = iter->begin();
        while ( siter != iter->end() )
		{
			delete *siter;
			++siter;
		}
		(*iter).clear();
		++iter;
	}
    
    for (int i=0; i< cellTable.size(); i++) cellTable[i].clear();
    cellTable.clear();
}



size_t cellcomplex::getDimension() const
{
	return dimension;
}


const vector< vector< cell*> >&
cellcomplex::getTheCellTable() const
{
    return cellTable;
}


// get the 0-cells that are contained in the k-cell c or in any
// of the (k-1)-cells, (k-2)-cells ... contained in c.
//
std::set< feliks::topology::cell* >
cellcomplex:: get0Cells(const feliks::topology::cell& c) const
{
    std::map< int, std::set<feliks::topology::cell*> > v;
    
    // fill p-level with the cell itself
    std::set<feliks::topology::cell*> s;
    s.insert(const_cast<feliks::topology::cell*>(&c));
    v[c.dimension()] = s;
    
    for (int p = c.dimension(); p>0; p--)
    {
        std::set<feliks::topology::cell*> s;
        BOOST_FOREACH( feliks::topology::cell* c, v[p])
        {
            feliks::topology::chain::const_iterator iter = c->boundary().begin();
            while ( iter != c->boundary().end() )
            {
                s.insert( (*iter).first );
                ++iter;
            }
        }
        v[p-1] = s;
    }
    
    return v[0];
}




// get the d-cells that are contained in the k-cell c or in any
// of the (k-1)-cells, (k-2)-cells ... contained in c.
//
std::set< feliks::topology::cell* >
cellcomplex:: getSubCells(const feliks::topology::cell& c, const int d) const
{
    std::map< int, std::set<feliks::topology::cell*> > v;
    std::set<feliks::topology::cell*> s;
    
    // quick return for trivial case
    if ( d > c.dimension() ) return s;
    
    // fill p-level with the cell itself
    s.insert(const_cast<feliks::topology::cell*>(&c));
    v[c.dimension()] = s;
    
    for (int p = c.dimension(); p>d; p--)
    {
        std::set<feliks::topology::cell*> s;
        BOOST_FOREACH( feliks::topology::cell* c, v[p])
        {
            feliks::topology::chain::const_iterator iter = c->boundary().begin();
            while ( iter != c->boundary().end() )
            {
                s.insert( (*iter).first );
                ++iter;
            }
        }
        v[p-1] = s;
    }
    
    return v[d];
}




void cellcomplex::print(std::ostream& of) const
{
    for (int i=0; i<= dimension; i++)
    {
        of << "\n ";
        std::set< cell *>::iterator siter = (*this)[i].begin();
        while (siter != (*this)[i].end() )
        {
            (*siter)->print(of);
            of << ", ";
            ++siter;
        }
    }
    of << std::flush;
}



void cellcomplex::printChains(std::ostream& of) const
{
    of << "\n";
    for (int i=0; i<= dimension; i++)
    {
        of << "\n " << i << "Cells:" ;
        std::set< cell *>::iterator siter = (*this)[i].begin();
        while (siter != (*this)[i].end() )
        {
            (*siter)->print(of);
            of << "\n    Boundary : ";
            (*siter)->boundary().print(of);
            of << "\n    coboundary : ";
            (*siter)->coboundary().print(of);
            of << "\n";
            ++siter;
        }
    }
    
}


chain cellcomplex::randomChain(const int d) const
{
	assert( d <= dimension);
	chain ch;
	int   nd_cells = (*this)[d].size();
	int   nch_cells = ccrandom(0, nd_cells );
	
	for (int i = 0; i< nch_cells; i++)
	{
		int pick = ccrandom(0, nd_cells);
		std::set<cell *>::const_iterator iter = (*this)[d].begin();
		advance(iter, pick);
		ch[ *iter ] = ccrandom(0,10);
	}
	return ch;
}


cochain cellcomplex::randomCochain(const int d) const
{
	assert( d <= dimension);
	cochain ch;
	int   nd_cells = (*this)[d].size();
	int   nch_cells = ccrandom(0, nd_cells );
	
	for (int i = 0; i< nch_cells; i++)
	{
		int pick = ccrandom(0, nd_cells);
		std::set<cell *>::const_iterator iter = (*this)[d].begin();
		advance(iter, pick);
		ch[ *iter ] = ccrandom(0,10);
	}
	return ch;
}




bool feliks::topology::runCellComplexTests(const cellcomplex& cc)
{
    bool ret = true;
    
	// test 1: boun boun = 0
	if (cc.getDimension() > 1)
	{
		vector< set<cell *> >::const_iterator iter = cc.begin() + 2;
		while ( iter != cc.end() )
		{
			set< cell *>::iterator citer = (*iter).begin();
			
			while ( citer != (*iter).end() )
			{
				chain ch = (*citer)->boundary();
				chain Dch = ch.boundary();
				
				if ( !Dch.isZero() )
                {
					cout << "\nTest failed: Cellcomplex 1";
                    ret = false;
                }
				++citer;
			}
			++iter;
		}
	}
	
	// test 2: coboun o coboun = 0
	if (cc.getDimension() > 1)
	{
		vector< set<cell *> >::const_iterator iter = cc.begin();
		while ( iter != cc.end()-2 )
		{
			set< cell *>::iterator citer = (*iter).begin();
			
			// test1: boundary o boundary = 0
			while ( citer != (*iter).end() )
			{
				chain ch = (*citer)->boundary();
				chain Dch = ch.boundary();
				
				if ( !Dch.isZero() )
                {
					cout << "\nTest failed: Cellcomplex 2";
                    ret = false;
                }
				++citer;
			}
			++iter;
		}
	}
	
	
	// test 3 < coch, boun ch> = < coboun coch, ch >
	{
		for (int i=0; i<cc.getDimension()-1; i++)
		{
			chain   ch = cc.randomChain(i+1);
			cochain co = cc.randomCochain(i);
			
			int t1 = co( ch.boundary() ) - ch( co.coboundary() );
			
			if ( t1 != 0)
			{
				cout << "\nTest failed: Cellcomplex 3";
                ret = false;
			}
		}
	}
	
	// test 4, boundary is homomorphism
	{
		for (int i=1; i<cc.getDimension(); i++)
		{
			chain ch1 = cc.randomChain(i);
			chain ch2 = cc.randomChain(i);
			
			
			chain test = (ch1+ch2).boundary() - ch1.boundary() - ch2.boundary();
			if ( !test.isZero())
            {
				cout << "\nTest failed: Cellcomplex 4";
                ret = false;
            }
		}
	}
	
	
	// test 5, coboundary is homomorphism
	{
		for (int i=0; i<cc.getDimension()-1; i++)
		{
			cochain ch1 = cc.randomCochain(i);
			cochain ch2 = cc.randomCochain(i);
			
			
			cochain test = (ch1+ch2).coboundary() - ch1.coboundary() - ch2.coboundary();
			if ( !test.isZero())
            {
				cout << "\nTest failed: Cellcomplex 5";
                ret = false;
            }
		}
	}
    
    
    // test 6, determine correctly the belonging
    {
        cell* cM  = *(cc[cc.getDimension()].begin());
        cell* cM1 = cM->boundary().begin()->first;
        
        if (  !(cM1->belongsToClosedCell(*cM)) )
        {
            cout << "\nTest failed: Cellcomplex 6";
            ret = false;
        }
        
        
        // test that cM2, the daughter of cM1, belongs to closedCell(cM)
        cell* cM2;
        if (cM1->dimension() > 0 && !(cM1->boundary().empty()) )
        {
            cM2 = cM1->boundary().begin()->first;
            if (  !(cM2->belongsToClosedCell(*cM)) )
            {
                cout << "\nTest failed: Cellcomplex 6, second level";
                ret = false;
            }
            

            // test that cM3, the daughter of cM2, belongs to closedCell(cM)
            cell* cM3;
            if (cM2->dimension() > 0 && !(cM2->boundary().empty()))
            {
                cM3 = cM2->boundary().begin()->first;
                if (  !(cM3->belongsToClosedCell(*cM)) )
                {
                    cout << "\nTest failed: Cellcomplex 6, third level";
                    ret = false;
                }
                
            }
        }
    }
    
	return ret;
}


static int ccrandom(int lower, int upper)
{
	int r;
	
	if (upper - lower == 0)
		r = lower;
	else
		r = lower +  rand() % (upper-lower);
	
	return r;
}


