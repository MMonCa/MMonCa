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
//  shapefun.h
//  feliks++
//
//  Created by Romero Ignacio on 12/25/11.
//  Copyright (c) 2011 Universidad Politécnica de Madrid. All rights reserved.
//

#pragma once
#ifndef feliks___shapefun_h
#define feliks___shapefun_h

#include "Math/tensor.h"
#include <iostream>

class node;

class shapefun
{
    
public:
                            shapefun();
                            shapefun(node* const p);
                            shapefun(const shapefun&);
    virtual                 ~shapefun();
    shapefun&               operator=(const shapefun& s);

    virtual void            print(std::ostream& of=std::cout) const;

    
    node                    *parentNode;
    double                  value;			// value of the shape function at the ev. point
    blue::ivector                 dxyz;			// derivatives w/r to global coordinates
    blue::itensor                 dd;
    bool                    secondDflag;
};



#endif
