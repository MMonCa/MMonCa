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
/* material.cpp
 *
 * iro
 * june 2000
 *
 */


#include "material.h"
#include <string>
#include <map>
#include <cstdlib>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <stdexcept>

#include "Main/feliks.h"
#include "Io/logger.h"
#include "Io/message.h"
#include "Io/usercommand.h"

#include "Materials/Smallstrain/smallstrainlib.h"


std::map<int,material*> material :: allMaterials;
int                     material :: largestMaterialLabel = 0;


 

material :: material() :
	_name("unknown name"),
    _label(0)
{
}



material :: material(const commandLine &cl) :
	_name("unknown name"),
    _label(0)
{
	for (size_t k=0; k < cl.size(); k++)
	{
		const usercommand &uc = cl[k];
			
		if      (uc.keyword() == "name")         _name    = uc.option();
        else if (uc.keyword() == "label")        _label   = uc.integerValue();
	}
}




material :: ~material()
{}




void material :: freeMaterials()
{
    std::map<int, material*>::iterator iter = allMaterials.begin();
	
	while( iter != allMaterials.end() )
	{
		delete iter->second;
		++iter;
	}
}




/* retrieves a material structure given its number in the static list of materials */
material& material :: getMaterial(int label)
{
	std::map<int,material*>::iterator iter = allMaterials.find(label);
	
	if (iter == allMaterials.end())
		ErrorMessage("Material not defined");
    return *(iter->second);
}



material& material :: getMaterial(const std::string& name)
{
	material *s(0);
	std::map<int,material*>::iterator iter = allMaterials.begin();
	
	while( iter != allMaterials.end() )
	{
		if ( iter->second->name() == name) 
		{
			s = iter->second;
			break;
		}
		++iter;
	}
	
	if (s == 0) throw std::runtime_error("Material not found in global list.");
	return *s;
	
}





size_t material :: getNMaterials()
{
	return allMaterials.size();
}





/* Purpose:   prints in screen all the materials in the material_list
*	if initialized is false, it means there are no materials in the analysis
*/
void material :: listMaterials(std::ostream &of)
{
	of  << "\n\n\n";
	std::stringstream title;
	title << "M a t e r i a l s  (" << allMaterials.size() << ")";
	std::string str(title.str());
	printCentered(of, str);
	
    std::map<int, material*>::iterator iter = allMaterials.begin();
	while( iter != allMaterials.end() )
	{
		of << "\n\nMaterial #" << iter->first << ": " << iter->second->name();
		(iter->second)->print(of);
		++iter;
	}
}

	


materialPoint :: materialPoint()
{}



materialPoint :: materialPoint(material &m)
{
}

