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
 *  initcondition.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 20/4/07.
 *
 */

#include "initcondition.h"
#include "Io/logger.h"
#include "Model/Node/node.h"
#include "Io/message.h"
#include "Analysis/Propfactors/propfactor.h"
#include "Model/model.h"
#include "Model/Sets/nodeset.h"
#include "Io/message.h"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <boost/foreach.hpp>

#define DEBUG_SCANNER 0

using namespace std;


#define DEBUG_SCANNER 0

/* the format for the boundary conditions is
initial_conditions (or ic)
...
...
...
[blank line]

or
initial_conditions, nodeset = [label]
dof1 value factor ...
[blank line]
*/
void initcondition :: scan(commandLine &cl, ifstream &meshfile, model &m)
{
	string oneLine;
	// input of standard bc, just a list of bc
    if (cl.size() < 2)
    {
        // scan boundary conditions definitions  until blank line or EOF
        while ( getline(meshfile, oneLine) )
        {
			// if the line has only spaces, it's the end of the block
			if ( oneLine.find_first_not_of(" ") == string::npos) break;
			
            // set the default values, so it does not take the last one if missing
			int nlabel		 = 0;
            int dofnumber    = -1;
            double bcvalue   = 0.0;
            int factorlabel  = 0;
			
			stringstream strm(oneLine);
			strm >> nlabel >> dofnumber >> bcvalue >> factorlabel;
            initcondition *ic = new initcondition(m.getNode(nlabel), dofnumber, bcvalue, factorlabel);
			if (DEBUG_SCANNER == 1) cout << "adding ic on node: " << nlabel << ", dof " << dofnumber << endl;
            if (ic != NULL) m.add(*ic);
        }
    }
	
	
	// second way of defining boundary conditions, via nodeset
    else
    {
        usercommand uc = cl[1];
        if ( uc.keyword() == "nodeset")
        {
            // get the nodeset on which the bc are applied
            nodeset &nds = m.getNodeset( uc.option() );
			
			// now read the bc and add them to bc list
			while ( getline(meshfile, oneLine) )
			{
				if ( oneLine.find_first_not_of(" ") == string::npos) break;
				
                // set the default values, so it does not take the last one if missing
                int    dofnumber    = -1;
                double bcvalue      = 0.0;
                int    factorlabel  = 0;
				
				stringstream strm(oneLine);
				strm >>  dofnumber >> bcvalue >> factorlabel;
                BOOST_FOREACH(node* nd, nds)
				{
					initcondition *ic = new initcondition( *nd, dofnumber, bcvalue, factorlabel);
					if ( Debug_Level() > 0 ) 
                        logger::mainlog << "\n Adding initial condition for node " << nd->getLabel();
					if (ic != NULL) m.add(*ic);
				}
			}
		}
	}
}
