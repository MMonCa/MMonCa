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
/* initrate.cpp
 *
 * initial rates.
 *
 * numbering of local degrees of freedom starts from 0.
 *
 * ignacio romero
 * march 2007
 */

#include "initrate.h"
#include "Io/io.h"
#include "Analysis/Propfactors/propfactor.h"

#include "Model/model.h"
#include "Model/Node/node.h"
#include "Model/Sets/nodeset.h"
#include "Model/Parts/body.h"

#include "Geometry/Solidmodel/BRep.h"
#include "Math/Topology/topology.h"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>


#include <boost/foreach.hpp>

#define DEBUG_SCANNER 0

using namespace std;





manifoldInitrate :: manifoldInitrate(const commandLine& cl)
:
nodesetInitrate(cl),
theBRepConstrainedCell(0)
{
    for (int k=1; k< cl.size(); k++)
    {
        const usercommand  &uc = cl[k];
        
        if      (uc.keyword() == "body")         bodyname   = uc.option();
        else if (uc.keyword() == "dimension")    dimension  = uc.integerValue();
        else if (uc.keyword() == "label")        mnfldLabel = uc.integerValue();
    }
}




manifoldInitrate :: ~manifoldInitrate()
{
    
}




void manifoldInitrate :: initialize(model& m)
{
}




void manifoldInitrate :: print(std::ostream& of) const
{
    of << "\n\n Initrate on a body boundary";
    of << "\n Name of the body   : " << bodyname;
    of << "\n Boundary dimension : " << dimension;
    of << "\n Boundary label     : " << mnfldLabel;
    
    if (u[0] != 0.0) of << "\n       ux imposed   : " << u[0];
	if (u[1] != 0.0) of << "\n       uy imposed   : " << u[1];
	if (u[2] != 0.0) of << "\n       uz imposed   : " << u[2];
	if (p    != 0.0) of << "\n       p  imposed   : " << p;
	if (r[0] != 0.0) of << "\n     rotx imposed   : " << r[0];
	if (r[1] != 0.0) of << "\n     roty imposed   : " << r[1];
	if (r[2] != 0.0) of << "\n     rotz imposed   : " << r[2];
	if (d[0] != 0.0) of << "\n     dirx imposed   : " << d[0];
	if (d[1] != 0.0) of << "\n     diry imposed   : " << d[1];
}






nodesetInitrate :: nodesetInitrate(const commandLine &cl) :
theNodeset(0),
theNodesetName()
{
    u[0] = u[1] = u[2] = 0.0;
    p    = 0.0;
    r[0] = r[1] = r[2] = 0.0;
    d[0] = d[1] = 0.0;

    
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		
		if      (uc.keyword() == "nodeset")	theNodesetName = uc.option();
        else if	(uc.keyword() ==  "ux")		u[0]   = uc.value();
		else if	(uc.keyword() ==  "uy")		u[1]   = uc.value();
		else if	(uc.keyword() ==  "uz")		u[2]   = uc.value();
		
		else if	(uc.keyword() ==  "p")		p      = uc.value();
		
		else if	(uc.keyword() ==  "rotx")	r[0]   = uc.value();
		else if	(uc.keyword() ==  "roty")	r[1]   = uc.value();
		else if	(uc.keyword() ==  "rotz")	r[2]   = uc.value();
		
		else if	(uc.keyword() ==  "dirx")	d[0]   = uc.value();
		else if	(uc.keyword() ==  "diry")	d[1]   = uc.value();
	}
}




void nodesetInitrate :: initialize(model& m)
{
	theNodeset = &(m.getNodeset(theNodesetName));
    m.add(*this);
}





void nodesetInitrate :: impose() const
{
	for (int i=0; i<3; i++)
    {
		if (u[i] != 0.0)
        {
			BOOST_FOREACH(node* n, *theNodeset)
                n->setInitialRate(i, u[i] );
        }
    }
    
    if (p != 0.0)
    {
        BOOST_FOREACH(node* n, *theNodeset)
        {
            n->constrainDof(3);
        }
    }
    
	for (int i=0; i<3; i++)
	{
        if (r[i] != 0.0)
        {
			BOOST_FOREACH(node* n, *theNodeset)
                n->constrainDof(3+i);
        }
    }
    
	for (int i=0; i<2; i++)
    {
		if (d[i] != 0.0)
        {
			BOOST_FOREACH(node* n, *theNodeset)
                n->constrainDof(3+i);
        }
    }
}



void nodesetInitrate :: print(std::ostream& of) const
{
	of << "\n\n Initial rate on nodeset " << theNodesetName;
    if (u[0] != 0.0) of << "\n       ux imposed   : " << u[0];
	if (u[1] != 0.0) of << "\n       uy imposed   : " << u[1];
	if (u[2] != 0.0) of << "\n       uz imposed   : " << u[2];
	if (p    != 0.0) of << "\n       p  imposed   : " << p;
	if (r[0] != 0.0) of << "\n     rotx imposed   : " << r[0];
	if (r[1] != 0.0) of << "\n     roty imposed   : " << r[1];
	if (r[2] != 0.0) of << "\n     rotz imposed   : " << r[2];
	if (d[0] != 0.0) of << "\n     dirx imposed   : " << d[0];
	if (d[1] != 0.0) of << "\n     diry imposed   : " << d[1];
}





