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
//  elmttasks.cpp
//  feliks++
//
//  Created by Romero Ignacio on 12/16/11.
//  Copyright (c) 2011 Universidad Politécnica de Madrid. All rights reserved.
//

#include <iostream>
#include "elmttasks.h"


elmtTasks :: elmtTasks() 
{
	reset();
}



void elmtTasks :: reset()
{
    residualtangent	 = false;
    tangent          = false;
	static_tangent   = false;
	dynamic_tangent  = false;
    gradient         = false;
    static_residual  = false;
    dynamic_residual = false;
    mass_matrix      = false;
	damping_matrix	 = false;
    error_estimate   = false;
    energy           = false;
	eigmax           = false;
	lump_mass        = false;
	L2_projection    = false;
	L2_projection_matrix = false;
}

