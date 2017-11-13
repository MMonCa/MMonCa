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
/* initcondition.h
 *
 * ignacio romero
 *
 * the data is identical to that of pointLoads, so we reuse all the functions. 
 * Other options would be
 *  i) to rewrite everything identical to pointload.h and pointload.c
 * ii) to eliminate the constraint data type and just consider pointloads over
 * constrained dofs differently.
 *
 * I have chosen to have a different datatype for clarity but identical implementation.
 * Data types and functions are then just alias to the ones in pointLoad.
 */


#ifndef _initcondition_h
#define _initcondition_h
	
#include "Analysis/Loading/pointload.h"
#include <iostream>

class commandLine;
class model;
class node;

class initcondition : public pointLoad {

public:
	
	initcondition();
	initcondition(node& aNode, int dofn, double F, int flabel) : pointLoad(aNode, dofn, F, flabel) {}
	
	static void scan(commandLine &cl, std::ifstream &meshfile, model &m);

};	

#endif
