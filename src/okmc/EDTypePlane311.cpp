/*
 * EDTypePlane311.cpp
 *
 *  Created on: Aug 14, 2013
 *
 * Author: ignacio.martin@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain
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

#include "EDTypePlane311.h"
#include "Particle.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"

using Kernel::Coordinates;
using std::map;

namespace OKMC {

EDTypePlane311::EDTypePlane311() : 	MAX_TEMPL_SIZE(2000)
{
	sort311Template();
}

Coordinates EDTypePlane311::surface(Kernel::RNG &, const Kernel::Coordinates c[2],
		const Kernel::Coordinates &center, const Kernel::ID &id, const Particle *pPart) const
{
	map<Kernel::P_POS, unsigned>::const_iterator it=id._pos.find(Kernel::POS_I);
	unsigned size = (it == id._pos.end()? 0 : it->second);
	if(size == 0)
		size = 1;  //correctly initialize.
	if(pPart == 0 || Domains::global()->PM()->isIorV(id._mt, pPart->getPType()) == Kernel::IV_I)
		return c[0]*_sortedTemplate311[size-1]._x+c[1]*_sortedTemplate311[size-1]._y;
	else
		return pPart->getCoordinates() - center;
}

const Particle * EDTypePlane311::emitFrom (Kernel::RNG &rng, const Kernel::Coordinates &center,
				const std::vector<Particle *> &parts) const
{
	return parts[rng.rand()*parts.size()];
}

void EDTypePlane311::sort311Template()
{
	std::vector<Kernel::Coordinates> genericTemplate;
	Kernel::Coordinates c;
	for (int i = -100; i < 100; i++)
		for (int j = -100; j < 100; j++)
		{
			c._x = float(i);
			c._y = float(j);
			c._z = 0;
			genericTemplate.push_back(c);
		}

	_sortedTemplate311.clear();
	std::vector <Kernel::Coordinates> template311;
	template311.clear();
	Kernel::Coordinates prop;
	//resizing Generic Template
	for (unsigned i = 0; i < genericTemplate.size(); i++)
	{
		prop._x = 0.543*sqrt(22.)/4.*genericTemplate[i]._x;
		prop._y = 0.543/sqrt(2.)*genericTemplate[i]._y;
		prop._z = genericTemplate[i]._z;
		template311.push_back(prop);
	}

	if(template311.size() != genericTemplate.size())
		ERRORMSG("Error while rearranging regular template for 311's");

	float maxX=0., maxY=0.,minX=0.,minY=0.;
	for(unsigned i = 0; i <template311.size(); i++)
	{
		if(maxX <template311[i]._x)
			maxX = template311[i]._x;
		if(maxY <template311[i]._y)
			maxY = template311[i]._y;
		if(minX >template311[i]._x)
			minX = template311[i]._x;
		if(minY >template311[i]._y)
			minY = template311[i]._y;

	}
	//Computing geometric center of new template
	float centerX = (maxX-fabs(minX))/2., centerY=(maxY-fabs(minY))/2.;
	float minDCX = 500000000., minDCY = 500000000.; //min distance to center.
	unsigned idx = 5635465;  //absurd number. If it does not change, it means that center is not found
	for(unsigned i = 0; i <template311.size(); i++)
		if(minDCX >= fabs(template311[i]._x-centerX) && minDCY >= fabs(template311[i]._y-centerY) )
		{
			minDCX = fabs(template311[i]._x-centerX);
			minDCY = fabs(template311[i]._y-centerY);
			idx = i;
		}
	if(idx==5635465)
		ERRORMSG("Failed while computing center of generic template in sortTemplate311");
	_sortedTemplate311.push_back(template311[idx]);  //first atom in the 311 is the center of the template and it will construct around it

	float yDist = 0.543/sqrt(2.), xDist = 0.543*sqrt(22.)/4.;  //distance between atoms in cols and rows respectively (x & y)
	float newXDist = xDist, newYDist=yDist;		//"surface" to fill.

	//Builds the template in order depending on "surfaces". Pushes back a new atom if the distance in x and y is proper
	//and it doesn't already exist, maintaining (aprox) the ratio between length and width:
	//W=sqrt(0.5*L) as in "Solid State Electronics 52 (2008) 1430-1436. Martin-Bragado et al.
	while(_sortedTemplate311.size() < template311.size())
	{
		for(unsigned i = 0; i <template311.size(); i++)
		{
			bool bPush = true; //controls if the atom exists in sortedTemplate so it's not pushed into the vector

			if( (fabs(_sortedTemplate311[0]._x-template311[i]._x) <= newXDist + 0.01 ) &&
					(fabs(_sortedTemplate311[0]._x-template311[i]._x) >= newXDist - 0.01) &&
					(fabs(_sortedTemplate311[0]._y-template311[i]._y) < newYDist + 0.01) )
			{
				for(unsigned j = 0; j < _sortedTemplate311.size(); j++)
					if(_sortedTemplate311[j]._x == template311[i]._x  &&_sortedTemplate311[j]._y == template311[i]._y  )
					{
						bPush = false;
						break;
					}
				if(bPush)
					_sortedTemplate311.push_back(template311[i]);
			}
			else if (fabs(_sortedTemplate311[0]._y-template311[i]._y) <= newYDist + 0.01 )
			{
				if (fabs(_sortedTemplate311[0]._y-template311[i]._y) >= newYDist - 0.01)
					if(fabs(_sortedTemplate311[0]._x-template311[i]._x) < newXDist + 0.01)
					{
						for(unsigned j = 0; j < _sortedTemplate311.size(); j++)
							if(_sortedTemplate311[j]._x == template311[i]._x  &&_sortedTemplate311[j]._y == template311[i]._y  )
							{
								bPush = false;
								break;
							}
						if(bPush)
							_sortedTemplate311.push_back(template311[i]);
					}
			}
		}

		//Computes the maximum and minimum of the current template to check the next surface to fill
		float minXinTemplate=5000., maxXInTemplate = 0., minYinTemplate=5000., maxYInTemplate = 0.;
		for(unsigned i = 0; i <_sortedTemplate311.size(); i++)
		{
			if(minXinTemplate > _sortedTemplate311[i]._x)
				minXinTemplate =_sortedTemplate311[i]._x;
			if(maxXInTemplate < _sortedTemplate311[i]._x)
				maxXInTemplate =_sortedTemplate311[i]._x;
			if(minYinTemplate > _sortedTemplate311[i]._y)
				minYinTemplate =_sortedTemplate311[i]._y;
			if(maxYInTemplate < _sortedTemplate311[i]._y)
				maxYInTemplate =_sortedTemplate311[i]._y;
		}

		if( (maxXInTemplate - minXinTemplate) < sqrt(0.5*(maxYInTemplate - minYinTemplate)))
		{
			newXDist += xDist;
		}
		else if((maxYInTemplate - minYinTemplate) < 0.5*pow(float(_sortedTemplate311.size()), 2./3.))
		{
			newYDist += yDist;
		}
		else
		{
			ERRORMSG("Failed to grow 311 template according to ratio between length and width");
		}

		if(_sortedTemplate311.size() > MAX_TEMPL_SIZE)
			break;
	}

}

}
