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
 *  cell.h
 *  feliks
 *
 *  Created by Ignacio Romero on 10/3/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */


#pragma once
#ifndef _topology_cell_h

#include "chain.h"
#include <iostream>
#include <set>

namespace feliks
{	
	namespace topology
	{		
		class cell
		{
			
		public:
                                cell(const int dim);
                                ~cell();
			
            bool                belongsToClosedCell(const cell&) const;
			const chain&        boundary() const;
			chain&              boundary();
			void                boundary2coboundary();
			cochain&            coboundary();
			const cochain&      coboundary() const;
			void                coboundary2boundary();
			const int           dimension() const;
			void                print(std::ostream& of=std::cout) const;
            std::set<cell* >    touchingCells(const int d) const;

            
            
		private:
                                cell();
                                cell(const cell&);
            
			int                 _dimension;
			chain               _boundary;
			cochain             _coboundary;
            
#ifdef DEBUG
			int                 _label;
			static int          _last_label[4];
#endif
		};
	}
}

#endif
