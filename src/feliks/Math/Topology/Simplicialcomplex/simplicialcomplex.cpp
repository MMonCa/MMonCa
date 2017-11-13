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
 *  simplicialcomplex.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 10/10/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "simplicialcomplex.h"
#include "../topology.h"
#include <cassert>
#include <utility>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include "boost/foreach.hpp"

using namespace feliks::topology;
using namespace std;



connectivityTable::connectivityTable() 
{
}



connectivityTable::connectivityTable(const int size)
{
    (*this).resize( size );
}



void connectivityTable::print(std::ostream& of) const
{
    int tet = 1;
    BOOST_FOREACH( vector<int> v, *this)
    {
        of << "\n (" << std::setw(3) << tet++ << ") --> [";
        BOOST_FOREACH(int i, v)  of << " " << setw(4) << i;
        of << " ]";
    }
    of << std::flush;
}



void abstractSimplex::removeComponent(const int i)
{
    assert(i< (*this).size() );
    abstractSimplex::iterator iter = (*this).begin();
    advance(iter,i);
    (*this).erase(iter);
}


void abstractSimplex::print(std::ostream& of) const
{
    abstractSimplex::iterator iter = (*this).begin();
    of << "[";
    while ( iter != (*this).end() )
    {
        (*iter)->print(of);
        of << ", ";
        ++iter;
    }
    of << "]";
}


simplicialcomplex :: simplicialcomplex(const int dim) :
cellcomplex(dim),
isCTbuilt(false)
{
    
}



const connectivityTable& simplicialcomplex :: getConnectivity() const
{
    if (!isCTbuilt) const_cast<simplicialcomplex&>(*this).buildConnectivity();
    return connectivity;
}


simplicialcomplex::~simplicialcomplex()
{
    
}




void simplicialcomplex :: buildConnectivity()
{
    std::map< feliks::topology::cell* , std::set< feliks::topology::cell*> > bigset;
    int d = this->dimension;
    
    // form set with cells of dimension d
    BOOST_FOREACH( cell* c, (*this)[d])
    {
        bigset[c].insert(c);
    }
    
    for (int p=d; p>0; --p)
    {
        // form set with cells of dimension p-1
        BOOST_FOREACH( cell* c, (*this)[d])
        {
            std::set< cell*>& smallset = bigset[c];
            BOOST_FOREACH( feliks::topology::cell* k, smallset )
            {
                const chain& ch = k->boundary();
                std::map<cell*, int>::const_iterator iter = ch.begin();
                while ( iter != ch.end() )
                {
                    smallset.insert( (*iter).first );
                    ++iter;
                }
                smallset.erase( k );
            }
        }
    }
    
    // the elements of bigset contain the sets of 0-cells
    // now they need to be numbered uniquely
    int maxlabel = -1;
    std::map< cell*, int> labels;
    {
        std::map< cell*, std::set<cell*> >::iterator iter = bigset.begin();
        while ( iter != bigset.end() )
        {
            BOOST_FOREACH(cell* c, (*iter).second)
            {
                if (labels.find(c) == labels.end()) labels[c] = ++maxlabel;
            }
            ++iter;
        }        
    }
    
    // write connectivity as array of integers
    connectivity.resize( bigset.size() );
    {
        std::map< cell*, std::set<cell*> >::iterator iter = bigset.begin();
        int cellnumber = 0;
        while ( iter != bigset.end() )
        {
            BOOST_FOREACH(cell* c, (*iter).second)
            {
                connectivity[cellnumber].push_back( labels[c] );
            }
            ++iter;
            ++cellnumber;
        }        
    }
}




void simplicialcomplex :: populate(const connectivityTable& ct)
{
    connectivity = ct;
    isCTbuilt    = true;
    
    if (ct.empty()) return;
    
    // verify that the vertex list starts from 1 and has no holes
    std::vector<int> full;
    int maxv, minv;
    full.reserve( ct.size() * ct[0].size() );
    {
        connectivityTable::const_iterator citer = ct.begin();
        while ( citer != ct.end() )
        {
            full.insert( full.end() , (*citer).begin() , (*citer).end() );
            ++citer;
        }
        std::sort( full.begin() , full.end() );
        full.erase( std::unique( full.begin() , full.end() ), full.end() );
        
        std::vector<int>::iterator mi = std::max_element( full.begin() , full.end() );
        minv = full[0];
        maxv = *mi;
        
        if ( minv != 1 )
            throw std::runtime_error("bad connectivity table");
        else if ( maxv != full.size() )
            throw std::runtime_error("bad connectivity table");
    }
    
    
    // create unique vertex list (this sets an ordering in stl)
    std::map<int, cell*> label2cell;
    for (int i=1; i<=maxv; i++)
        label2cell[i] = new cell(0);
    
    
    // store vertices by label in cellTable[0]
    cellTable.resize( this->getDimension() + 1 );
    cellTable[0].resize( label2cell.size() );
    for (int i=0; i< label2cell.size(); i++)
        cellTable[0][i] = label2cell[i+1];
    
    // create abstract connectivity table as set< set < cell * > >
    set< abstractSimplex > abstract_connectivity;
    {
        connectivityTable::const_iterator citer = ct.begin();
        while ( citer != ct.end() )
        {
            abstractSimplex tet;
            vector<int>::const_iterator viter = (*citer).begin();
            while (viter != (*citer).end() )
            {
                tet.insert( label2cell[ *viter ] );
                ++viter;
            }
            
            // in the abstract_connectivity, the tets are ordered according
            // to internal rules, not to insertion sequence
            abstract_connectivity.insert( tet );
            
            ++citer;   
        }
    }
    
    
    // create abstract simplicial complex
    vector< set< abstractSimplex > > S;
    S.resize(dimension+1);
    {
        S[dimension] = abstract_connectivity;
        for (int p=dimension; p>0; p--)
        {
            set< abstractSimplex >::iterator piter = S[p].begin();
            while (piter != S[p].end() )
            {
                for (int i=0; i<=p; i++)
                {
                    abstractSimplex  simplex = *piter;
                    simplex.removeComponent(i);
                    S[p-1].insert(simplex);
                }
                ++piter;
            }
        }
    }
    
    // print for debugging
    if (false)
    {
        cout << "\n Abstract simplicial complex\n";
        for (int p=0; p<=dimension; p++)
        {
            int k=1;
            set< abstractSimplex >::iterator piter = S[p].begin();
            while (piter != S[p].end() )
            {
                cout << "(" << k++ << "): ";
                (*piter).print(); cout << "\n";
                ++piter;
            }
        }
    }
    
    // allocate maps from S to K
    vector< map< abstractSimplex , cell *> > g;
    g.resize(dimension+1);
    
    
    // initialize g map for 0-cells
    {
        set< abstractSimplex >::iterator piter = S[0].begin();
        while (piter != S[0].end() )
        {
            abstractSimplex vertex = *piter;
            cell*           c      = *(vertex.begin());
            (*this)[0].insert(c);
            g[0][vertex] = c;
            ++piter;
        }
    }
    
    
    // define simplicial complex
    for (int p=1; p<=dimension; p++)
    {
        set<abstractSimplex>::iterator piter = S[p].begin();
        while (piter != S[p].end() )
        {
            cell* c = new cell(p);
            (*this)[p].insert(c);
            g[p][*piter] = c;
            cellTable[p].push_back(c);
            
            ++piter;
        }
    }
    
    
    // print for debugging
    if (false)
    {
        cout << "\n\nSimplicial complex";
        (*this).print();
    }
    
    
    // print for debugging
    if (false)
    {
        cout << "\n\nAbstract mapping g\n";
        for (int p=0; p<=dimension; p++)
        {
            set<abstractSimplex>::iterator piter = S[p].begin();
            while (piter != S[p].end() )
            {
                (*piter).print();
                cout << " --> ";
                g[p][*piter]->print(); cout << "\n";
                
                ++piter;
            }
        }
    }
    
    
    // define boundary and coboundary operator
    for (int p=1; p<=dimension; p++)
    {
        set<abstractSimplex >::iterator piter = S[p].begin();
        while (piter != S[p].end() )
        {
            cell* sigma = g[p][*piter];
            
            int sign = -1;
            for (int i=0; i<=p; i++)
            {
                sign *= -1;
                
                abstractSimplex simplex = *piter;
                simplex.removeComponent(i);
                cell* f = g[p-1][simplex];
                
                sigma->boundary()[f] = f->coboundary()[sigma] = sign;
            }
            ++piter;
        }
    }
}




void simplicialcomplex :: populate(const std::set< abstractSimplex >& abstract_connectivity)
{
    isCTbuilt    = false;
    if (abstract_connectivity.empty()) return;
    
    // allocate cellTable
    cellTable.resize( this->getDimension() + 1);
    
    
    // create abstract simplicial complex
    vector< set< abstractSimplex > > S;
    S.resize(dimension+1);
    {
        S[dimension] = abstract_connectivity;
        for (int p=dimension; p>0; p--)
        {
            set< abstractSimplex >::iterator piter = S[p].begin();
            while (piter != S[p].end() )
            {
                for (int i=0; i<=p; i++)
                {
                    abstractSimplex  simplex = *piter;
                    simplex.removeComponent(i);
                    S[p-1].insert(simplex);
                }
                ++piter;
            }
        }
    }
    
    // print for debugging
    if (false)
    {
        cout << "\n Abstract simplicial complex\n";
        for (int p=0; p<=dimension; p++)
        {
            set< abstractSimplex >::iterator piter = S[p].begin();
            while (piter != S[p].end() )
            {
                (*piter).print(); cout << "\n";
                ++piter;
            }
        }
    }
    
    // allocate maps from S to K
    vector< map< abstractSimplex , cell *> > g;
    g.resize(dimension+1);
    
    
    // initialize g map for 0-cells
    {
        set< abstractSimplex >::iterator piter = S[0].begin();
        while (piter != S[0].end() )
        {
            abstractSimplex vertex = *piter;
            cell*           c      = *(vertex.begin());
            (*this)[0].insert(c);
            g[0][vertex] = c;
            cellTable[0].push_back(c);
            ++piter;
        }
    }
    
    
    // define simplicial complex
    for (int p=1; p<=dimension; p++)
    {
        set<abstractSimplex>::iterator piter = S[p].begin();
        while (piter != S[p].end() )
        {
            cell* c = new cell(p);
            (*this)[p].insert(c);
            g[p][*piter] = c;
            cellTable[p].push_back(c);

            ++piter;
        }
    }
    
    
    // print for debugging
    if (false)
    {
        cout << "\n\nSimplicial complex";
        (*this).print();
    }
    
    
    // print for debugging
    if (false)
    {
        cout << "\n\nAbstract mapping g\n";
        for (int p=0; p<=dimension; p++)
        {
            set<abstractSimplex>::iterator piter = S[p].begin();
            while (piter != S[p].end() )
            {
                (*piter).print();
                cout << " --> ";
                g[p][*piter]->print(); cout << "\n";
                
                ++piter;
            }
        }
    }
    
    
    // define boundary and coboundary operator
    for (int p=1; p<=dimension; p++)
    {
        set<abstractSimplex >::iterator piter = S[p].begin();
        while (piter != S[p].end() )
        {
            cell* sigma = g[p][*piter];
            
            int sign = -1;
            for (int i=0; i<=p; i++)
            {
                sign *= -1;
                
                abstractSimplex simplex = *piter;
                simplex.removeComponent(i);
                cell* f = g[p-1][simplex];
                
                sigma->boundary()[f] = f->coboundary()[sigma] = sign;
            }
            ++piter;
        }
    }
    
    
    // form simple connectivity table
    {
        std::map< cell*, int> c2i;
        int label = 1;
        std::set<cell*>::iterator  it = (*this)[0].begin();
        while ( it != (*this)[0].end() )
        {
            c2i[ *it ] = label++;
            ++it;
        }
        
        std::set< abstractSimplex >::iterator eiter = abstract_connectivity.begin();
        while ( eiter != abstract_connectivity.end() )
        {
            std::vector<int> el;
            abstractSimplex::iterator viter = eiter->begin();
            while ( viter != eiter->end() )
            {
                el.push_back( c2i[*viter]);
                ++viter;
            }
            connectivity.push_back(el);
            
            ++eiter;
        }
        
    isCTbuilt = true;
    }
}





bool feliks::topology::runSimplicialComplexTests()
{
    bool ret = true;
    
    // test 1: a triangle
    {
        connectivityTable tri(1);
        
        tri[0].push_back(1);
        tri[0].push_back(2);
        tri[0].push_back(3);
        
        simplicialcomplex sc(2);
        sc.populate(tri);
        if (!runCellComplexTests(sc))
        {
            cout << "\nTest failed: simplicialcomplex 1";
            ret = false;
        }
    }
    
    
    // test 2 a tetrahedron
    {
        connectivityTable tet;
        tet.resize(1);
        
        tet[0].push_back(1);
        tet[0].push_back(2);
        tet[0].push_back(3);
        tet[0].push_back(4);
        
        simplicialcomplex sc(3);
        sc.populate(tet);
        if (!runCellComplexTests(sc))
        {
            cout << "\nTest failed: simplicialcomplex 2";
            ret = false;
        }
    }
    
    
    // a cube
    /*
     *          
     *   5----------8
     *   |\         |\
     *   | \        | \
     *   |  \       |  \
     *   |   6----------7
     *   |   |      |   | 
     *   1---|------\4  | 
     *    \  |       \  |
     *     \ |        \ |
     *      \|         \|
     *       2----------3
     *
     * tets:
     
     4 5 2 7     // center one
     6 2 5 7     // 6 opposite to 4
     2 3 4 7     // 3 opposite to 5
     4 7 8 5     // 8 opposite to 2
     5 4 2 1     // 1 opposite to 7
     
     *
     */
    {
        connectivityTable meshedCube;
        meshedCube.resize(5);
        
        meshedCube[0].push_back(4);
        meshedCube[0].push_back(5);
        meshedCube[0].push_back(2);
        meshedCube[0].push_back(7);
        
        meshedCube[1].push_back(6);
        meshedCube[1].push_back(2);
        meshedCube[1].push_back(5);
        meshedCube[1].push_back(7);
        
        meshedCube[2].push_back(2);
        meshedCube[2].push_back(3);
        meshedCube[2].push_back(4);
        meshedCube[2].push_back(7);
        
        meshedCube[3].push_back(4);
        meshedCube[3].push_back(7);
        meshedCube[3].push_back(8);
        meshedCube[3].push_back(5);
        
        meshedCube[4].push_back(5);
        meshedCube[4].push_back(4);
        meshedCube[4].push_back(2);
        meshedCube[4].push_back(1);
        
        simplicialcomplex sc(3);
        sc.populate(meshedCube);
        if (!runCellComplexTests(sc))
        {
            cout << "\nTest failed: simplicialcomplex 3";
            ret = false;
        }
    }
    
    if (true)
    {
        feliks::topology::cell* v1 = new cell(0);
        feliks::topology::cell* v2 = new cell(0);
        feliks::topology::cell* v3 = new cell(0);
        feliks::topology::cell* v4 = new cell(0);
        feliks::topology::cell* v5 = new cell(0);
        feliks::topology::cell* v6 = new cell(0);
        feliks::topology::cell* v7 = new cell(0);
        feliks::topology::cell* v8 = new cell(0);
        
        std::set< abstractSimplex > ct;
        abstractSimplex t0;
        t0.insert(v4);
        t0.insert(v5);
        t0.insert(v2);
        t0.insert(v7);
        ct.insert(t0);
        
        abstractSimplex t1;
        t1.insert(v6);
        t1.insert(v2);
        t1.insert(v5);
        t1.insert(v7);
        ct.insert(t1);
        
        abstractSimplex t2;
        t2.insert(v2);
        t2.insert(v3);
        t2.insert(v4);
        t2.insert(v7);
        ct.insert(t2);
        
        abstractSimplex t3;
        t3.insert(v4);
        t3.insert(v7);
        t3.insert(v8);
        t3.insert(v5);
        ct.insert(t3);
        
        abstractSimplex t4;
        t4.insert(v5);
        t4.insert(v4);
        t4.insert(v2);
        t4.insert(v1);
        ct.insert(t4);
        
        
        simplicialcomplex sc(3);
        sc.populate(ct);
        if (!runCellComplexTests(sc))
        {
            cout << "\nTest failed: simplicialcomplex 4";
            ret = false;
        }        
    }
    
    return ret;
}



