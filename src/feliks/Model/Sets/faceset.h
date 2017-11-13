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
 *  faceset.h
 *  feliks
 *
 *  Created by Ignacio Romero on 31/3/10.
 *  Copyright 2010 Universidad Politécnica de Madrid. All rights reserved.
 *
 */


#ifndef _faceset_h
#define _faceset_h

#include <string>
#include <vector>
#include <iostream>

#include "Model/elmtface.h"

class faceset
{

public:
	faceset();
	faceset(const std::string& n);
	~faceset();
	
	void		add(const elmtface &f);
	void		add(const faceset &fs);
	void		eliminateDuplicatedFaces();
	bool		empty() const;
	const std::string&	getName() const;
	nodeset&	getNodes();
	void		info(std::ostream& of=std::cout) const;
	void		print(std::ostream& of=std::cout) const;
	void		sort();
	void		setName(const std::string& n);
	size_t		size() const;
	const elmtface&	operator[](const int k);

	
	// iterators to enable stl-like operations with facesets
	std::vector<elmtface>::iterator begin();
	std::vector<elmtface>::iterator end();
	
	std::vector<elmtface>::const_iterator begin() const;
	std::vector<elmtface>::const_iterator end() const;
	
	
private:	
	std::string						name;
	bool							sorted;
	std::vector<elmtface>			faces;
	nodeset							nodesInFaceset;
};



inline std::vector<elmtface>::iterator			faceset :: begin()			{return faces.begin();}
inline std::vector<elmtface>::iterator			faceset :: end()			{return faces.end();}
inline std::vector<elmtface>::const_iterator		faceset :: begin() const	{return faces.begin();}
inline std::vector<elmtface>::const_iterator		faceset :: end() const		{return faces.end();}
inline bool									faceset :: empty() const	{return faces.empty();}
inline size_t								faceset :: size() const		{return faces.size();}
inline const elmtface&						faceset :: operator[](const int k){return faces[k];}
inline const std::string&					faceset	:: getName() const	{return name;}

#endif

