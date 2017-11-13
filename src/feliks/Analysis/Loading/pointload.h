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
/* pointload.h
 *
 * single point external load
 *
 * ignacio romero
 *
 * converted to c++ in jan 2006
 */

	

#ifndef _pointload_h
#define _pointload_h

#include <iostream>
#include "Analysis/Loading/loading.h"


class energies;
class node;
class commandLine;
class model;

class pointLoad : public loading{
	
protected:
    int     dof;				// local dof of the node. Starts from 0;
    double  value;				// the reference value of the load.
    int     nodelabel;			// label of the node where the load is applied
    node    *nd;
		
public:
	enum    loadtype{ displacement, pressure, rotation, director};

	
	pointLoad(node& aNode, int dofn, double F, int flabel);
	virtual ~pointLoad(){}
	
	static void scan(commandLine &cl, std::ifstream &meshfile, model &m);

	void	accumulateEnergy(energies& ener, const double time);
	void	assembleLoading(assembler& asb, const double time, const double f);
	void	initialize(model& theModel);
	
	void    replaceNode(node& deleted, node& survivor);
	void    transferLoadsToNodes(const double time);
	
	void    print(std::ostream &of=std::cout) const;
	friend  std::ostream&  operator<<(std::ostream &os, const pointLoad &load);

	int     getDof() const       {return dof;}
	node   &getNode() const      {return *nd;}
	int     getNodeLabel() const {return nodelabel;}
	double  getValue()	const	{return value;}	
	
private:
	pointLoad();
};





#endif
