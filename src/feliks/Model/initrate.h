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
/* initrate.h
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


#ifndef _initrate_h
#define _initrate_h
	
#include "Model/Sets/nodeset.h"


#include <iostream>
#include <memory>
#include <string>
#include <list>



class node;
class commandLine;
class model;




class initrate
{

public:
    virtual                     ~initrate(){}
    
    virtual void                impose() const = 0;
    virtual void                initialize(model& m) = 0;
    virtual void                print(std::ostream& of=std::cout) const = 0;
};	




class nodesetInitrate: public initrate
{
	
public:
                                nodesetInitrate(const commandLine &cl);
    virtual                     ~nodesetInitrate(){}
    
    
    virtual void                impose() const;
    virtual void                initialize(model& m);
    virtual void                print(std::ostream& of=std::cout) const;
    
	
protected:
    nodeset*                    theNodeset;
    std::string                 theNodesetName;
    double                      u[3], p, r[3], d[2];
};




class manifoldInitrate : public nodesetInitrate
{
    
public:
                                manifoldInitrate(const commandLine& cl);
    virtual                     ~manifoldInitrate();
    
    virtual void                initialize(model& m);
    virtual void                print(std::ostream& of=std::cout) const;
    
    
private:
    
    std::string                 bodyname;
    int                         dimension, mnfldLabel;
    feliks::topology::cell*       theBRepConstrainedCell;
    nodeset                     privateNodeset;
};




#endif
