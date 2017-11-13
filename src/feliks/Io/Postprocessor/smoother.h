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
 *  smoother.h
 *  feliks
 *
 *  Created by Ignacio Romero on 11/11/08.
 *  Copyright 2008.
 *
 */

#ifndef _smoother_h
#define _smoother_h

#include <iostream>
#include <string>

class model;
class node;
class resultdata;
class body;
class poorbody;



class smoother
{
    
public:
    virtual                 ~smoother(){};
	virtual std::string		getName() const=0;
	virtual bool			smoothToNode(resultdata& r, node& nd, body& bd)=0;
	virtual bool			smoothToNode(resultdata& r, node& nd, poorbody& bd)=0;
};



// Superconvergent Patch Recovery 
class SPRsmoother : public smoother
{
    
public:
    void	printNodalGradients(const model& m, const std::string& gradname, std::ostream& of=std::cout);
    bool	smoothToNode(resultdata& r, node& nd, body& md);
    bool	smoothToNode(resultdata& r, node& nd, poorbody& md);
    virtual std::string	getName() const {return "SPR smoother";}
};


// simple lumping scheme 
class lumpsmoother : public smoother
{	
public:
    bool	smoothToNode(resultdata& r, node& nd, body& md);
    bool	smoothToNode(resultdata& r, node& nd, poorbody& md);
    virtual std::string	getName() const {return "Lump smoother";}
};


#endif

