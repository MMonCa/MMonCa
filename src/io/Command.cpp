/*
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

#include "Command.h"
#include "Diagnostic.h"
#include "AnnealCmd.h"
#include "CascadeCmd.h"
#include "InitCmd.h"
#include "InsertCmd.h"
#include "LowMsgCmd.h"
#include "SaveCmd.h"
#include "TestCmd.h"
#include "ExtractCmd.h"
#include "ParamCmd.h"
#include "ProfileCmd.h"
#include "ReportCmd.h"
#include "RestartCmd.h"
#include "ClusterCmd.h"
#include <cstring>
#include <cstdlib>

using std::vector;
using std::map;
using std::string;
using std::pair;
namespace IO
{

Command::Command(Tcl_Interp *tcl, int argc, const char *argv[], bool bPrint) : _pTcl(tcl)
{
	std::stringstream ss;
	std::string parameter;
	for(int i=1; i<argc; ++i)
		parameter = parameter + argv[i] + ' ';
	
	//------ line for command
	for(unsigned i=0; i< 79 - std::strlen(argv[0]) - 7; ++ i)
		ss << "-";
	ss << " " << argv[0] << " -----\n";
	ss << argv[0];
	unsigned count=0;
	bool quotes = false;
	string arg;
	for(string::iterator it=parameter.begin(); it!=parameter.end(); ++it)
	{
		switch(*it)
		{
		case '{':
			if(count != 0)
				arg += '{';
			count++;
			break;
		case '}':
			count--;
			if(count != 0)
				arg += '}';
			break;
		case '"':
			if(count == 0)
				quotes = !quotes;
			else
				arg += '"';
			break;
		case ' ':
			if(arg.size() && count == 0 && quotes == false)
			{
				vector<string> tokens;
				separateTokens(arg.c_str(), '=', tokens);
				if(tokens[0].size() == 0)
					ERRORMSG("Parameter " << arg << " is not well defined");
				if(tokens.size() == 1)
				{
					insert(tokens[0], "", "");
					_used[tokens[0]] = false;
					ss << " " << tokens[0];
				}
				else if(tokens.size() == 2)
				{
					insert(tokens[0], "", tokens[1]);
					_used[tokens[0]] = false;
					if(tokens[1].size() < 100)
						ss << " " << tokens[0] << "='" << tokens[1] << "'";
					else
						ss << " " << tokens[0] << "='" << tokens[1][0] << tokens[1][1] << tokens[1][2] << " (...)'";
				}
				else
					ERRORMSG("Parameter '" << arg << "' contains more than 2 tokens");
				arg.clear();
				break;
			}
		default:
			arg += *it;
			break;
		}
	}

	if(specified("msg.low"))
		Domains::global()->setVerbosity(0);
	else if(specified("msg.medium"))
		Domains::global()->setVerbosity(1);
	else if(specified("msg.high"))
		Domains::global()->setVerbosity(2);

	ss << std::endl;
	for(int i=0; i< 79; ++ i)
		ss << "-";
	ss << std::endl;
	if(bPrint && !specified("no.print"))
		LOWMSG(ss.str());
}

Command::~Command()
{
	for(map<string, bool>::iterator it=_used.begin(); it!= _used.end(); ++it)
		if(it->second == false)
			WARNINGMSG("Parameter " << it->first << " not used");
}

int Command::Anneal(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	AnnealCmd anneal(interp, argc, argv);
	return anneal();
}

int Command::Cluster(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	ClusterCmd cluster;
	return cluster(interp, argc, argv);
}

int Command::Save(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	SaveCmd save(interp, argc, argv);
	return save();
}

int Command::Cascade(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	CascadeCmd cascade(interp, argc, argv);
	return cascade();
}

int Command::Init(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	InitCmd init(interp, argc, argv);
	return init();
}

int Command::Insert(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	InsertCmd insert(interp, argc, argv);
	return insert();
}

int Command::Exit(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	Domains::global()->deleteGlobal();
	int exit_code;
	if(argc == 1)
		exit_code = 0;
	else
		exit_code = atoi(argv[1]);
	exit(exit_code);
	return TCL_OK;
}

int Command::Extract(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	ExtractCmd extract(interp, argc, argv);
	return extract();
}

int Command::License(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	LOWMSG("");
	LOWMSG(" OKMC: Copyright (C) 2014 IMDEA Materials Institute, Getafe, Madrid, Spain");
	LOWMSG("       Main author: I. Martin-Bragado. ignacio.martin@imdea.org");
	LOWMSG(" LKMC: Copyright (C) 2014 IMDEA Materials Institute, Getafe, Madrid, Spain");
	LOWMSG("       Main author: I. Martin-Bragado. ignacio.martin@imdea.org");
	LOWMSG(" FEM:  Copyright (C) 2014 IMDEA Materials Institute, Getafe, Madrid, Spain and");
	LOWMSG("       Copyright (C) 2014 Technical University of Madrid (UPM), Madrid, Spain");
	LOWMSG("       Main author: I. Romero. ignacio.romero@imdea.org");
	LOWMSG("       FEM code based on FELIKS");
	LOWMSG(" BCA:  Copyright (C) 2000-2014 University of Valladolid (UVA), Valladolid, Spain and");
	LOWMSG("       Copyright (C)      2014 IMDEA Materials Institute, Getafe, Madrid, Spain.");
	LOWMSG("       Main author: J. Hernandez-Mangas. jesus.hernandez.mangas@tel.uva.es");
	LOWMSG("       BCA code based on Ion Implant Simulator (IIS, IIS-UVA)");
	LOWMSG("");
	LOWMSG("Licensed under the Apache License, Version 2.0 (the \"License\");");
	LOWMSG("you may not use this file except in compliance with the License.");
	LOWMSG("You may obtain a copy of the License at");
	LOWMSG("");
	LOWMSG("http://www.apache.org/licenses/LICENSE-2.0");
	LOWMSG("");
	LOWMSG("Unless required by applicable law or agreed to in writing, software");
	LOWMSG("distributed under the License is distributed on an \"AS IS\" BASIS,");
	LOWMSG("WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.");
	LOWMSG("See the License for the specific language governing permissions and");
	LOWMSG("limitations under the License.");

	return TCL_OK;
}

int Command::LowMsg(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	LowMsgCmd lmsg(interp, argc, argv);
	return lmsg();
}

int Command::Param(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	ParamCmd setParam(interp, argc, argv);
	return setParam();
}

int Command::Profile(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	ProfileCmd profile(interp, argc, argv);
	return profile();
}

int Command::Report(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	ReportCmd report(interp, argc, argv);
	return report();
}

int Command::Restart(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	RestartCmd restart(interp, argc, argv);
	return restart();
}

int Command::Test(ClientData, Tcl_Interp *interp, int argc, const char *argv[])
{
	TestCmd test(interp, argc, argv);
	return test();
}

//only takes into account the first occurrence of the token.
// vec.size() <= 2 always.
void Command::separateTokens(const char *str, char delimiter, vector<string> &vec)
{
	vec.clear();
	if(*str)
		vec.push_back("");
	char p=*str;
	while(p)
	{
		if(p == delimiter && vec.size() < 2)
			vec.push_back("");
		else
			vec.back() += p;
		p=*++str;
	}
}

}
