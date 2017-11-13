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
 *  chain.h
 *  feliks
 *
 *  Created by Ignacio Romero on 10/9/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#ifndef _chain_h
#define _chain_h

#include <vector>
#include <set>
#include <map>
#include <iostream>


namespace feliks
{
	namespace topology
	{		
		
		class cell;
		class cochain;
		
		// virtual class for chain and cochain
		class baseChain : public std::map<cell *,int>
		{
			
		public:
                                baseChain(const baseChain&);
                                baseChain& operator=(const baseChain&);
			virtual             ~baseChain();
			
			bool				isZero() const;
			virtual void		print(std::ostream& of=std::cout) const;
			
                
		protected:
                                baseChain(const int d);
                                baseChain();
			int					dimension;
		};
		
		
        
		
		class chain : public baseChain
		{
		public:
                            chain();
                            chain(const int d);
			chain			boundary();
			int				operator()(const cochain& co);
			void			print(std::ostream& of=std::cout) const;
			
			chain&			operator+=(const chain &v);
			chain&			operator-=(const chain &v);
			chain&			operator*=(const int a);
			chain			operator*(const int a);
			chain			operator+(const chain &right);
			chain			operator-(const chain &right);
		};
		
		
        
        
		class cochain : public baseChain
		{
		public:
                            cochain();
                            cochain(const int d);
			cochain			coboundary();
			int				operator()(const chain& co);
			cochain			operator+(const cochain &right);
			cochain&		operator+=(const cochain &v);
			cochain			operator-(const cochain &right);
			cochain&		operator-=(const cochain &v);
			cochain			operator*(const int a);
			cochain&		operator*=(const int a);
			

			void			print(std::ostream& of=std::cout) const;
		};
	}
}

#endif

