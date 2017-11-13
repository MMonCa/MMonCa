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
//  elmttask.h
//  feliks++
//
//  Created by Romero Ignacio on 12/16/11.
//  Copyright (c) 2011 Universidad Politécnica de Madrid. All rights reserved.
//

#pragma once
#ifndef feliks___elmttask_h
#define feliks___elmttask_h

// this is type that contains the information about which task(s)
// must be performed when the element functions are called 

class elmtTasks
{
	
public:    
    bool residualtangent;		// Compute both residual and tangent
    bool tangent;               // Compute element tangent
    bool static_tangent;        // Compute the static part of the tangent
    bool dynamic_tangent;       // Compute the dynamic part of the tangent
    bool gradient;              // Compute element gradients at integration point
    bool static_residual;		// Compute the static part of the element residual
    bool dynamic_residual;		// Compute the dynamic part of the element residual
    bool mass_matrix;			// Compute consistent mass
    bool damping_matrix;		// Compute consistent mass
    bool error_estimate;		// Error estimation
    bool energy;         		// Energy, momenta, etc
	bool eigmax;				// maximum eigenvalue for explicit transients
	bool lump_mass;				// lump nodal masses
	bool L2_projection;			// RHS for L2 projection
	bool L2_projection_matrix;
    
	elmtTasks();
	void reset();
};


#endif
