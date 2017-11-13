/*
 * Param.cpp
 *
 *  Created on: Feb 15, 2011
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

#include "ParamCmd.h"
#include "kernel/Domain.h"
#include "io/FileParameters.h"
#include "domains/Global.h"
#include "io/ParameterManager.h"

using std::string;

namespace IO {

ParamCmd::ParamCmd(Tcl_Interp *tcl, int argc, const char *argv[]) : Command(tcl, argc, argv)
{
}

int ParamCmd::operator ()()
{
	if(specified("set"))
	{
		string type = getString("type");
		string key  = getString("key");
		string value= getString("value");
		bool   bNew = specified("new");
		string index= specified("index")? getString("index") : "";

		if(!bNew && !Domains::global()->getFileParameters()->specified(key, type))
			ERRORMSG("Parameter " << type << " "  << key << " does not exist in database");
		Domains::global()->getFileParameters()->setCommand(type, key, value, index);
		return TCL_OK;
	}
	if(specified("add"))
	{
		string type = getString("type");
		string key  = getString("key");
		string value= getString("value");
		string index= specified("index")? getString("index") : "";
		if(!Domains::global()->getFileParameters()->specified(key, type))
			ERRORMSG("Parameter " << type << " "  << key << " does not exist in database");
		Domains::global()->getFileParameters()->addCommand(type, key, value, index);
		return TCL_OK;
	}
	if(specified("get"))
	{
		string type   =                       getString("type");
		string key    =                       getString("key");
		string index  =  (specified("index")? getString("index") : "");

		if(!Domains::global()->getFileParameters()->specified(key, type))
			ERRORMSG("Parameter " << type << " "  << key << " does not exist in database");

		Tcl_AppendResult(_pTcl, Domains::global()->getFileParameters()->getCommand(type, key, index).c_str(), 0);
		return TCL_OK;
	}
	if(specified("unset"))
	{
		string type   =                       getString("type");
		string key    =                       getString("key");
		string index  = (specified("index")?  getString("index") : "");

		if(!Domains::global()->getFileParameters()->specified(key, type))
			ERRORMSG("Parameter " << type << " " << key << " does not exist in database");

		Domains::global()->getFileParameters()->unsetCommand(type, key, index);
		return TCL_OK;
	}
	if(specified("get.reaction"))
	{
		string material = getString("material");
		string index = getString("index");
		std::vector<string> tokens;  //for the two incoming species
		IO::ParameterManager::getTokens(index, '+', tokens);
		if(tokens.size() != 2)
			ERRORMSG("get.reactions: Two reactants needed for " << index);
		Kernel::M_TYPE mt = Domains::global()->PM()->getMaterialNumber(material);
		if(mt == Kernel::MAX_MATERIALS)
			ERRORMSG("get.reactions: Incorrect material " << material);

		Tcl_AppendResult(_pTcl, Domains::global()->isReaction(mt, tokens[0], tokens[1])? "true" : "false", 0);
		return TCL_OK;
	}
	ERRORMSG("param needs either 'set', 'add', 'get.reaction' or 'get' options");
	return TCL_ERROR;
}

}
