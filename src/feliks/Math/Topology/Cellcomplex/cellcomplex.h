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
 *  cellcomplex.h
 *  feliks
 *
 *  Created by Ignacio Romero on 10/3/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#pragma once
#ifndef _cellcomplex_h
#define _cellcomplex_h

#include <vector>
#include <set>
#include <map>
#include <iostream>

#include "cell.h"



namespace feliks
{
	namespace topology
	{		
		class cellcomplex : public std::vector< std::set<cell *> >
		{
		public:
			virtual                             ~cellcomplex();

            size_t                              getDimension() const;
            const std::vector< std::vector<cell*> >&  getTheCellTable() const;
            std::set< cell* >                   get0Cells(const cell& ) const;
            std::set< cell* >                   getSubCells(const cell&c, const int dim) const;
            void                                print(std::ostream& of=std::cout) const;
            void                                printChains(std::ostream& of=std::cout) const;
			chain                               randomChain(const int d) const;
			cochain                             randomCochain(const int d) const;  

						
		protected:
			cellcomplex(const int dimen);
            std::vector< std::vector< cell*> >  cellTable;
			size_t                              dimension;
			
		private:
			cellcomplex();
		};
		
				
		bool runCellComplexTests(const cellcomplex&);
	}	
}
#endif

