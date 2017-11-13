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
 *  simplicialcomplex.h
 *  feliks
 *
 *  Created by Ignacio Romero on 10/10/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#pragma once
#ifndef _simplicialcomplex_h
#define _simplicialcomplex_h

#include "Math/Topology/Cellcomplex/cellcomplex.h"
#include <vector>
#include <set>
#include <iostream>


namespace feliks
{
	namespace topology
	{
        class cell;
        
		class connectivityTable : public std::vector< std::vector< int > >
		{
            public:
                                            connectivityTable();
                                            connectivityTable(const int size);
            
            void                            print(std::ostream& of=std::cout) const;			
		};
        
		
		
        class abstractSimplex : public std::set<cell *>
        {
            public:
                void                        removeComponent(const int i);
                void                        print(std::ostream& of=std::cout) const;
        };
		
        
        
		class simplicialcomplex : public cellcomplex
		{
            public:
                                            simplicialcomplex(const int dimen);
                virtual                     ~simplicialcomplex();
            
                void                        populate(const connectivityTable&);
                void                        populate(const std::set< abstractSimplex >&);
                const connectivityTable&    getConnectivity() const;
            
            private:
                void                        buildConnectivity();
                connectivityTable           connectivity;
                bool                        isCTbuilt;
		};
        
        
        bool runSimplicialComplexTests();
	}
}


#endif

