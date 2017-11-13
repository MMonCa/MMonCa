/*
 * ClusterCmd.cpp
 *
 *  Created on: May 22, 2012
 *      Author: ignacio.martin@imdea.org
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

#include "ClusterCmd.h"
#include <string>
#include <vector>
#include <tcl.h>
#include "ParameterManager.h"

using std::string;

namespace IO {

int ClusterCmd::operator()(Tcl_Interp *p, int argc, const char *argv[])
{
	if(argc < 2)
	{
		Tcl_AppendResult(p, "cluster: Need an option", 0);
		return TCL_ERROR;
	}
	if(string(argv[1]) == "index")
	{
		if(argc != 4)
		{
			Tcl_AppendResult(p, "Not enough parameters. Use:  cluster index Species Cluster", 0);
			return TCL_ERROR;
		}
		string species(argv[2]);
		string cluster(argv[3]);

		string::size_type i=cluster.find(species);
		if(i == string::npos)
		{
			Tcl_AppendResult(p, "0", 0);
			return TCL_OK;
		}
		string index;
		for(string::iterator it=cluster.begin()+i+species.size(); it != cluster.end(); ++it)
			if(isdigit(*it))
				index += *it;
			else
				break;
		if(index.size())
			Tcl_AppendResult(p, index.c_str(), 0);
		else
			Tcl_AppendResult(p, "1", 0);
		return TCL_OK;
	}
	if(string(argv[1]) == "emit")
	{
		if(argc != 3)
		{
			Tcl_AppendResult(p, "Not enough parameters. Use:  cluster emit string", 0);
			return TCL_ERROR;
		}
		std::vector<string> tokens;
		ParameterManager::getTokens(argv[2], ',', tokens);
		Tcl_AppendResult(p, tokens[1].c_str(), 0);
		return TCL_OK;
	}
	return TCL_ERROR;
}

} /* namespace IO */
