/*
 * LatticeBCC1Param.cpp
 *
 *  Created on: Aug 17, 2012
 *      Author: ignacio.martin@imdea.org
 */

#include "LatticeBCC1Param.h"
#include "io/ArrheniusAlloys.h"
#include "io/ParameterManager.h"
#include "io/FileParameters.h"
#include <string>
#include <vector>
#include <map>

using std::vector;
using std::map;
using std::string;

namespace LKMC {

LatticeBCC1Param::LatticeBCC1Param(const IO::ParameterManager *pPM, const IO::FileParameters *pFP, Kernel::M_TYPE mt)
: LatticeParam(pPM, pFP, mt)
{
	_type = BCC;

	string base             = pPM->getMaterialName(mt) + "/Epitaxy/";
	//Epitaxy
	pFP->loadProcedure(Domains::global()->getTcl(), base+"activation", 0);
	map<string, IO::Arrhenius> actMap = pFP->getArrheniusProc(Domains::global()->getTcl(), base+"activation");
	//process information
	for(unsigned i=0; i<=FIRST; ++i)
		for(unsigned j=0; j<=SECOND; ++j)
			_activation[i][j] = IO::Arrhenius(0,5);
	for(map<string, IO::Arrhenius>::iterator it = actMap.begin(); it!=actMap.end(); ++it)
	{
		string index = it->first;
		vector<string> tokens;
		pPM->getTokens(index, ',', tokens);
		if(tokens.size() != 2)
			ERRORMSG(base + "activation" << ": Needs two tokens in " << index);
		std::stringstream ss;
		ss << tokens[0] << " " << tokens[1];
		unsigned nei1, nei2;
		ss >> nei1 >> nei2;
		//if(ss.failbit)
			//ERRORMSG(base + "activation" << ": Syntax error in " << index);
		if(nei1 > FIRST || nei2 > SECOND)
			ERRORMSG(base + "activation" << ": Incorrect number of neighbors in " << index);
		_activation[nei1][nei2] = it->second;
	}
}

LatticeBCC1Param::~LatticeBCC1Param() {

}

} /* namespace Domains */
