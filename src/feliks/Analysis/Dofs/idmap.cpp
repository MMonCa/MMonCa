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
 *  IDmap.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 22/6/09.
 *  Copyright 2009 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include "idmap.h"

#include <iostream>
#include <iomanip>


IDmap :: IDmap()
{
	_map[0] = _map[1] = _map[2] = NOTDEFINED; 
}

IDmap :: IDmap(const int ndof)
	: _ndof(ndof)
{ 
	_map[0] = _map[1] = _map[2] = NOTDEFINED; 
}


bool IDmap :: constrain(const int l)
{
	bool ret(false);
	
	if ( !isConstrained(l) )
	{
		_map[l] = CONSTRAINED;
		ret = true;
	}
	return ret;
}



// reset all dof numbers, leaving those that are known to be
// constrained 
void IDmap :: reset()
{
	if (_map[0] > CONSTRAINED ) _map[0] = NOTDEFINED;
	if (_map[1] > CONSTRAINED ) _map[1] = NOTDEFINED;
	if (_map[2] > CONSTRAINED ) _map[2] = NOTDEFINED;
}


std::ostream&  operator<<(std::ostream &os, const IDmap &mm)
{
	for (int k=0; k<mm.getNDOF(); k++)
	{
		
		if      (mm.isUndefined(k)   )	os << "    UNDEFINED ";
		else if (mm.isConstrained(k) )	os << "  CONSTRAINED ";
		else							os << "   " << std::setw(9) << mm(k);
	}
	return os;
}
