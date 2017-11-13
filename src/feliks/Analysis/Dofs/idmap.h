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
 *  idmap.h
 *  feliks
 *
 *  Created by Ignacio Romero on 22/6/09.
 *
 *  
 */


#ifndef _idmap_h
#define _idmap_h


#include <iostream>
#include <iomanip>

// labels for nodal degrees of freedom
#define NOTDEFINED     -2
#define CONSTRAINED    -1

class IDmap
	{
	
	public:
		IDmap();
		IDmap(const int ndof);
		
		int			getNDOF() const						{return _ndof;}
		int			getDof(const int l) const			{return _map[l];}
		
		bool		constrain(const int l);
		bool		isConstrained(const int l) const	{return _map[l] == CONSTRAINED;}
		bool		isUndefined() const					{return isUndefined(0) && isUndefined(1) && isUndefined(2);}
		bool		isUndefined(const int l) const		{return _map[l] == NOTDEFINED;}
		void		reset();
		void		setSize(const int ndofs)			{_ndof = ndofs;}
		
		int&		operator()(const int l)				{return _map[l];}
		const int&  operator()(const int l)	const		{return _map[l];}

		friend std::ostream&  operator<<(std::ostream &os, const IDmap &mm);
				
		
	private:
		int _ndof;
		int _map[3];
		
		friend class DOFnumberer;
	};


#endif
