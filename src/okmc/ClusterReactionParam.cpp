/*
 * ClusterReactionParam.cpp
 *
 *  Created on: Feb 25, 2014
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


#include "ClusterReactionParam.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include "io/FileParameters.h"
#include "ClusterParam.h"

#include <vector>

using std::vector;
using std::string;
using Kernel::M_TYPE;

namespace OKMC {

int ClusterInter::operator()(unsigned size0, unsigned size1) const
{
	if(size0 == size1 && _defectType[0] != -1)
		return _defectType[0];
	if(fabs(float(size0) - float(size1)) / std::max(size0,size1) < _approx && _defectType[1] != -1)
		return _defectType[1];
	if(size0 < size1)
		return _defectType[2];
	if(size0 > size1)
		return _defectType[3];
	return -1;
}

ClusterReactionParam::ClusterReactionParam(const ClusterParam * pClPar, const IO::FileParameters * pPar)
{
	IO::ParameterManager *pPM = Domains::global()->PM();
	for(M_TYPE mt=0; mt < pPM->getNMaterials(); ++mt)
	{
		string base_param = Domains::global()->PM()->getMaterialName(mt) + '/' + "Models/interaction.result";
		vector<string> reactions = pPar->getStrings(base_param, 5);
		unsigned noDef = pClPar->defectSize(mt);
		for(vector<string>::iterator it=reactions.begin(); it!=reactions.end(); ++it)
		{
			vector<string> tokens;
			IO::ParameterManager::getTokens(*it, ' ', tokens);
			if(tokens.size() != 5)
				ERRORMSG(base_param << ". Wrong syntax in " << *it);
			unsigned defect[3];
			for(unsigned i=0; i<3; ++i)
			{
				defect[i] = pClPar->getDefectNumber(mt, tokens[2*i]);
				if(defect[i] == noDef)
					ERRORMSG(base_param << ". Defect not recognized: " << tokens[2*i]);
			}

			while(_interactions[mt].size() <= defect[0])
				_interactions[mt].push_back(vector<ClusterInter>());
			while(_interactions[mt][defect[0]].size() <= defect[1])
				_interactions[mt][defect[0]].push_back(ClusterInter());
			while(_interactions[mt].size() <= defect[1])
				_interactions[mt].push_back(vector<ClusterInter>());
			while(_interactions[mt][defect[1]].size() <= defect[0])
				_interactions[mt][defect[1]].push_back(ClusterInter());

			// ~= could be ~= or ~=,0.05, etc...
			vector<string> tokens_operator;  //for the two incoming species
			IO::ParameterManager::getTokens(tokens[1], ',', tokens_operator);
			if(tokens_operator[0] == "+")
				for(unsigned i=0; i<4; ++i)
				{
					_interactions[mt][defect[0]][defect[1]]._defectType[i] = defect[2];
					if(defect[1] != defect[0])
					{
						if(_interactions[mt][defect[1]][defect[0]]._defectType[i] != -1)
							warning(base_param, mt, defect[0], defect[1], " + ", " + ");
						_interactions[mt][defect[1]][defect[0]]._defectType[i] = defect[2];
					}
				}
			else if(tokens_operator[0] == "==")
			{
				_interactions[mt][defect[0]][defect[1]]._defectType[0] = defect[2];
				if(defect[1] != defect[0])
				{
					if(_interactions[mt][defect[1]][defect[0]]._defectType[0] != -1)
						warning(base_param, mt, defect[0], defect[1], " == ", " == ");
					_interactions[mt][defect[1]][defect[0]]._defectType[0] = defect[2];
				}
			}
			else if(tokens_operator[0] == "~=")
			{
				if(tokens_operator.size() != 2)
					ERRORMSG(base_param << ". Wrong syntax in " << *it << ". It should be like \"~=,.05\" with .05 being the approximate proportion");
				_interactions[mt][defect[0]][defect[1]]._defectType[1] = defect[2];
				_interactions[mt][defect[0]][defect[1]]._approx = std::atof(tokens_operator[1].c_str());
				if(defect[1] != defect[0])
				{
					if(_interactions[mt][defect[1]][defect[0]]._defectType[1] != -1)
						warning(base_param, mt, defect[0], defect[1], " ~= ", " ~= ");
					_interactions[mt][defect[1]][defect[0]]._defectType[1] = defect[2];
					_interactions[mt][defect[1]][defect[0]]._approx = std::atof(tokens_operator[1].c_str());
				}

			}
			else if(tokens_operator[0] == "<")
			{
				_interactions[mt][defect[0]][defect[1]]._defectType[2] = defect[2];
				if(defect[1] != defect[0])
				{
					if(_interactions[mt][defect[1]][defect[0]]._defectType[3] != -1)
						warning(base_param, mt, defect[0], defect[1], " < ", " > ");
					_interactions[mt][defect[1]][defect[0]]._defectType[3] = defect[2];
				}
			}
			else if(tokens_operator[0] == ">")
			{
				_interactions[mt][defect[0]][defect[1]]._defectType[3] = defect[2];
				if(defect[1] != defect[0])
				{
					if(_interactions[mt][defect[1]][defect[0]]._defectType[2] != -1)
						warning(base_param, mt, defect[0], defect[1], " > ", " < ");
					_interactions[mt][defect[1]][defect[0]]._defectType[2] = defect[2];
				}
			}
			else
				ERRORMSG(base_param << ". Operator not recognized: " << tokens[1]);
		}
		//check the reactions
		for(unsigned def0 = 0; def0 < pClPar->defectSize(mt); ++def0)
			for(unsigned def1=0; def1 <= def0; ++def1)
			{
				if(pClPar->getParams(mt, def0)->_interactMC.find(def1) != pClPar->getParams(mt, def0)->_interactMC.end()) //if allowed
				{
					if(_interactions[mt].size() <= def0 || _interactions[mt][def0].size() <= def1)
						ERRORMSG(base_param << ". Interactions not defined for " << pPM->getEDName(mt, def0) << "+" << pPM->getEDName(mt, def1));
					if(_interactions[mt][def0][def1](50, 50) == -1)
						ERRORMSG(base_param << ". Interaction == or ~= not defined for " << pPM->getEDName(mt, def0) << "+" << pPM->getEDName(mt, def1));
					if(_interactions[mt][def0][def1](50, 1) == -1)
						ERRORMSG(base_param << ". Interaction > not defined for " << pPM->getEDName(mt, def0) << "+" << pPM->getEDName(mt, def1));
					if(_interactions[mt][def0][def1](1, 50) == -1)
						ERRORMSG(base_param << ". Interaction == or ~= not defined for " << pPM->getEDName(mt, def0) << "+" << pPM->getEDName(mt, def1));
					if(_interactions[mt][def0][def1](1,50) != _interactions[mt][def1][def0](50, 1))
						ERRORMSG(base_param << ". Interaction " <<
								pPM->getEDName(mt, def0) << " < " << pPM->getEDName(mt, def1) << " does not correspond to " <<
								pPM->getEDName(mt, def1) << " > " << pPM->getEDName(mt, def0));
					if(_interactions[mt][def0][def1](50,50) != _interactions[mt][def1][def0](50, 50))
						ERRORMSG(base_param << ". Interaction " <<
								pPM->getEDName(mt, def0) << " == " << pPM->getEDName(mt, def1) << " does not correspond to " <<
								pPM->getEDName(mt, def1) << " == " << pPM->getEDName(mt, def0));
					if(_interactions[mt][def0][def1](10000,9999) != _interactions[mt][def1][def0](9999, 10000))
						ERRORMSG(base_param << ". Interaction " <<
								pPM->getEDName(mt, def0) << " ~= " << pPM->getEDName(mt, def1) << " does not correspond to " <<
								pPM->getEDName(mt, def1) << " ~= " << pPM->getEDName(mt, def0));
				}
			}
	}
}

void ClusterReactionParam::warning(string param, M_TYPE mt, unsigned def0, unsigned def1, string op1, string op2)
{
	IO::ParameterManager *pPM = Domains::global()->PM();
	WARNINGMSG(param << ". Both " <<
			pPM->getEDName(mt, def0) << op1 << pPM->getEDName(mt, def1) <<
			" and " <<
			pPM->getEDName(mt, def1) << op2 << pPM->getEDName(mt, def0) <<
			" are implemented.");
}

} /* namespace OKMC */
