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
/* pointload.cpp
 *
 * single point load.
 *
 * numbering of local degrees of freedom starts from 0.
 *
 * ignacio romero
 * june 2000, translated to C++ in jan 2006
 */


#include "Analysis/Loading/pointload.h"
#include "Analysis/Loading/loading.h"
#include "Analysis/assembler.h"
#include "Analysis/Propfactors/propfactor.h"
#include "Analysis/energies.h"
#include "Io/logger.h"
#include "Model/Node/node.h"
#include "Model/Sets/nodeset.h"
#include "Model/model.h"
#include "Math/vector.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <set>

#define DEBUG_SCANNER 0

using namespace std;



// note that when defining loads, the numbering of the dofs starts from 1, like the user tells
// we only initialize the labels. The pointers to the node and propfactor must be filled later
pointLoad :: pointLoad(node& aNode, int dofn, double f, int flabel) :
	loading(flabel),
	nd(&aNode)
{
	theType = LOAD_POINT;
     if (dofn < 1)
     {
		 nodelabel   = -1;
		 dof         = -1;
		 logger::warnings << "\n\tWarning: local degree of freedom should be larger or equal to 1";
		 logger::warnings << "\n\t         Boundary condition / load ignored.";
     }

     else
     {
         nodelabel  = nd->getLabel();
         dof        = dofn-1;
         value      = f;
     }
}




void pointLoad :: accumulateEnergy(energies& ener, const double time)
{
    cout << "\n accumulate Energy in pointLoad should not be called";
	if (factor == 0) factor = &(factorcombo::getPropfactorCombinationFromLabel(getPropfactorLabel()));

	double f = factor->eval(time, nd->coordinates(dofset::tn1) ) * value;
    if (dof <3)
    {
        ener.potential -= nd->displacement(dofset::tn1)[dof] * f;
        ener.energy    -= nd->displacement(dofset::tn1)[dof] * f;
    }
}



void pointLoad :: assembleLoading(assembler& asb, const double time, const double f)
{
	if (time <= 0.0 || f <= 0.0) return;
	asb.assemble(*nd, dof, f*value*factor->eval( time, nd->coordinates() ));
}



void pointLoad :: initialize(model& theModel)
{
	factor = &(factorcombo::getPropfactorCombinationFromLabel(getPropfactorLabel()));
	nd->setHasLoad();	
}




/* prints the information of a point load.
* The local degree of freedom is incremented by one, since it is more natural for
* the user
*/
ostream&  operator<<(ostream& os, const pointLoad &load)
{
	os  << endl << setw(6) << load.nodelabel 
		<< setw(3) << load.dof+1 
		<< load.value << load.getPropfactorLabel() << flush;
	return os;
}




/* prints the information of a point load.
 * The local degree of freedom is incremented by one, since it is more natural for
 * the user
 */
void pointLoad :: print(ostream &of) const
{
	of  << endl 
		<< setw(4) << nodelabel << "    "
		<< setw(2) << dof+1     << "   "
		<< setw(9) << value     << "   "
		<< setw(3) << getPropfactorLabel() << flush;
}



void pointLoad :: replaceNode(node& deleted, node& survivor)
{
	if (nd == &deleted) nd = &survivor;
}



void pointLoad :: scan(commandLine &cl, ifstream &meshfile, model &m)
{
	string oneLine;
	
	
	// input of standard pointload, just a list of bc
    if (cl.size() < 2)
    {
        // scan pointLoad definitions  until blank line or EOF
        while ( commandLine::scanCommandLine(meshfile, oneLine) && oneLine.size() > 0)
        {
            // set the default values, so it does not take the last one if missing
			int    nlabel		= 0;
            int    dofnumber    = -1;
            double bcvalue      = 0.0;
            int    factorlabel  = 0;
			
			stringstream strm(oneLine);
			strm >> nlabel >> dofnumber >> bcvalue >> factorlabel;
            pointLoad *pl = new pointLoad(m.getNode(nlabel), dofnumber, bcvalue, factorlabel);
			if (DEBUG_SCANNER == 1) 
                cout << oneLine << endl 
                    << "adding pointload node: " << nlabel << ", dof " << dofnumber 
                    << ", val " << bcvalue << ", sc " << factorlabel << endl;
            if (pl != NULL) m.add(*pl);
        }
    }
	
	
	// second way of defining loads, via nodeset
    else
    {
        usercommand uc = cl[1];
        if ( uc.keyword() == "nodeset")
        {
            // get the nodeset on which the point are applied
			nodeset *nds;
			nds = &(m.getNodeset(uc.option()));
			
			// now read the bc and add them to bc list
			while ( commandLine::scanCommandLine(meshfile, oneLine)  && oneLine.size() > 0)
			{
                // set the default values, so it does not take the last one if missing
                int    dofnumber    = -1;
                double bcvalue      = 0.0;
                int    factorlabel  = 0;
				
				stringstream strm(oneLine);
				strm >> dofnumber >> bcvalue >> factorlabel;
                std::set<node*>::iterator iter = nds->begin();
                while ( iter != nds->end() )
                {
					pointLoad *pl = new pointLoad( **iter, dofnumber, bcvalue, factorlabel);
					if (pl != NULL) m.add(*pl);
                    ++iter;
				}
			}
		}
	}
}





void   pointLoad :: transferLoadsToNodes(const double time)
{
	// pointloads are not initialized at the beginning and the pointer to the factor
	// is left to null. The first time we use them, we initilize the pointer, for faster
	// access
	if (factor == 0) cout << "Error in transferLoadsToNodes pointloads" << endl;
	nd->assembleForceDof(dof, factor->eval(time, nd->coordinates() ) * value);
}






