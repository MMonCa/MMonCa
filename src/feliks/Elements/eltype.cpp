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
/* eltype.cpp
 *
 *
 *  FELIKS
 *
 *  Created by Ignacio Romero on Nov 2001.
 *  Copyright (c) 2003 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include "Elements/eltype.h"
#include "Elements/Interpolation/shapefun.h"
#include "Io/io.h"

#include "Math/Topology/topology.h"
#include "Math/tensor.h"

#include "Model/Node/node.h"

#include "Elements/Solid/Smallstrain/solid.h"
#include <string>
#include <sstream>
#include <ostream>
#include <iomanip>
#include <stdexcept>

using namespace std;
using namespace blue;


// initialization of the global list of element types
map <int, eltype*> eltype :: allEltypes;
int                eltype :: largestEltypeLabel=-1;




eltype :: eltype() :
	_name("unknown eltype"), 
    formulation("unknown formulation"),
    explanation("no explanation") ,
	_materialLabel(-1), 
	_geometry(NOPLOT),
	geoNames(0),
	geoDimensions(0)
{
}




eltype :: eltype(commandLine &cl) :
	_name("unknown eltype"), 
	formulation("unknown formulation"),
	explanation("no explanation") ,
	_materialLabel(-1), 
	_geometry(NOPLOT),
	geoNames(0),
	geoDimensions(0)
{
	// fill up material data if there is any
	for (int j=1; j< cl.size(); j++)
	{
		// first retrive the property from the input data
		const usercommand &uc = cl[j];
		
		// check if it is the material assignment
		if ( uc.keyword() == "material") _materialLabel = (int) uc.value();
	}
}




eltype :: ~eltype()
{
	geoDimensions.clear();
}




void eltype :: addToList()
{
	// insert it in the map
	pair< map<int,eltype*>::iterator, bool> insertCheck = allEltypes.insert( make_pair(_label, this) );
	
	// check that the label has never been used before
	if ( !insertCheck.second )
		logger :: mainlog << "\nElement type " << _label << " has alredy been defined. Use a different label.\n";
	
	// this is the definition when everything goes ok	
	// update counters
	largestEltypeLabel = (largestEltypeLabel > _label) ? largestEltypeLabel : _label;
}


evalspot* eltype :: createEvalspot(const double volume, const std::vector<shapefun>& theShp, const double initialspacing) const
{
    throw std::runtime_error("in createEvalspot");
    return 0;
}



node* eltype ::	createNode(const int label, const ivector& coor) const 
{
	node* ret(0);
	
	if			( nodetype() == _Pnode)     ret = new Pnode (label, coor);
	else if		( nodetype() == _Unode)     ret = new Unode (label, coor);
	else if		( nodetype() == _UPnode)	ret = new UPnode(label, coor);
	else if		( nodetype() == _UDnode)	ret = new UDnode(label, coor);
	
	return ret;
}


node* eltype ::	createNode(const int label, const ivector& coor, const feliks::topology::cell *the0Cell) const 
{
	node* ret(0);
	
	if			( nodetype() == _Pnode)     ret = new Pnode (label, coor, the0Cell);
	else if		( nodetype() == _Unode)     ret = new Unode (label, coor, the0Cell);
	else if		( nodetype() == _UPnode)	ret = new UPnode(label, coor, the0Cell);
	else if		( nodetype() == _UDnode)	ret = new UDnode(label, coor, the0Cell);
	
	return ret;
}




int	eltype :: getDofsPerNode() const
{	
	int ret = 0;
	
	if			( nodetype() == _Pnode)     ret = 1;
	else if		( nodetype() == _Unode)     ret = 3;
	else if		( nodetype() == _UPnode)	ret = 4;
	else if		( nodetype() == _URnode)	ret = 6;
	else if		( nodetype() == _UDnode)	ret	= 5;
	else if		( nodetype() == _EMPTYnode)	ret	= 0;
	
	return ret;
}



size_t eltype :: getNEltypes()
{
	return allEltypes.size();
}




void eltype :: freeEltypes()
{
	map<int, eltype*>::iterator iter = allEltypes.begin();
	
	// go one by one deallocating the eltypes
	while (iter != allEltypes.end())
	{
		delete iter->second;
		++iter;
	}
	
	// then erase the map itself
	allEltypes.erase( allEltypes.begin(), allEltypes.end());
}




eltype&  eltype ::  getEltypeFromLabel(int label)
{
	map<int, eltype*>::iterator iter = allEltypes.find(label);
	
	if (iter == allEltypes.end()) 
    {
        logger::mainlog << "\nERROR : Element type number " << label << " not defined\n" << endl;
        throw std::runtime_error("Error in getEltypeFromLabel");
    }
	return *(iter->second);
}




void eltype :: listEltypes(ostream &of)
{
	stringstream title;
	title << "E l e m e n t   t y p e s  (";
	title << allEltypes.size() << ")";
	string str(title.str());

	of  << "\n\n";
	printCentered(of, str);
	
	
	if (allEltypes.size() == 0)
	{
		of << "\n No element types defined in this model" << endl;
		return;
	}
	
	map<int, eltype*>::iterator iter = allEltypes.begin();
	while (iter != allEltypes.end() )
    {
		iter->second->print(of);
		++iter;
	}
	of << "\n";	
}




void eltype :: print(ostream &of)
{
	of	<< "\n\nElement type " << setw(3) << _label << " : "    << _name
    << "\n     "  << formulation
    << "\n     "  << explanation;
	if (_materialLabel > 0)
		of << "\n     Linked to material with label " << _materialLabel;
	else
		of << "\n     Element type is not linked to any material (warning).";
	of.precision(6);
	for (int k=0; k< geoDimensions.size(); k++)
		of  << "\n\t"   << setw(15) << geoNames[k] << " = " 
        << setw(10) << geoDimensions[k] << flush;
}




/* this fills up one entry of the elementtypelabel (the one defined by the input file)
 * table with an elementtype.
 * The commandlist must be of the form "eltype = 1, type = beam, thickness = 0.4, rotation = yes... "
 */
void eltype :: scan(const commandLine &ccl)
{	
	
}




