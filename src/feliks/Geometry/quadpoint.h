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
 * quadrature points function
 *
 * i. romero, feb 2002
 *
 */



#pragma once
#ifndef _quadpoint_h
#define _quadpoint_h
	

#include "Elements/eltype.h"


class quadpoint
{

public:

	double      coor[3];
	double      weight;

	double      getWeight() const;
    void        print() const;
};


inline double quadpoint :: getWeight() const { return weight;}





void    fullQuadraturePoints(       const size_t nnode, const geometryT geo, size_t *ngp, quadpoint *q);
size_t  numberOfSPRQuadraturePoints(const size_t nnode, const geometryT geo);
size_t  numberOfQuadraturePoints(   const size_t nnode, const geometryT geo);
void    quadraturePoints(const size_t nnode, const geometryT geo, const size_t ngp, quadpoint *q);
void    SPRQuadraturePoints(        const size_t nnode, const geometryT geo, const bool corner, size_t *ngp, quadpoint *q);




#endif



