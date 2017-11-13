/*
 * Element.h
 *
 *  Created on: Apr 17, 2013
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

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <string>

struct Element
{
	Element(double m, const std::string &s, unsigned e) : _mass(m), _name(s), _element(e) {}
	Element() : _mass(0), _element(0) {}
	double _mass;
	std::string _name;
	unsigned _element;
	Kernel::M_TYPE _mt;
};


#endif /* ELEMENT_H_ */
