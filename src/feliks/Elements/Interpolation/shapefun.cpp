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
//
//  shapefun.cpp
//  feliks++
//
//  Created by Romero Ignacio on 12/25/11.
//  Copyright (c) 2011 Universidad Politécnica de Madrid. All rights reserved.
//

#include "shapefun.h"
#include <iostream>
#include "Model/Node/node.h"
#include "Math/Topology/topology.h"


shapefun::shapefun() :
    value(0.0),
    parentNode(0),
    secondDflag(false)
{
    dxyz.setZero();
    dd.setZero();
}


shapefun::shapefun(node* const p) :
    value(0.0),
    parentNode(p),
    secondDflag(false)
{
    dxyz.setZero();
    dd.setZero();
}


shapefun::shapefun(const shapefun& N) :
    value(N.value),
    parentNode(N.parentNode),
    dxyz( N.dxyz ),
    secondDflag( N.secondDflag)
{
    if (secondDflag) dd = N.dd;
}



shapefun&  shapefun::operator=(const shapefun& N)
{
    value       = N.value;
    parentNode  = N.parentNode;
    dxyz        = N.dxyz;
    secondDflag = N.secondDflag;
    if (secondDflag) dd = N.dd;
    return *this;
}


shapefun::~shapefun()
{
    
}


void shapefun ::  print(std::ostream& of) const
{
    of << "\n\n Value at Gauss point  = " << std::scientific << std::right << value;
    of << "\n Derivative w/r to x   = " << dxyz[0];
    of << "\n Derivative w/r to y   = " << dxyz[1];
    of << "\n Derivative w/r to z   = " << dxyz[2];    
}
