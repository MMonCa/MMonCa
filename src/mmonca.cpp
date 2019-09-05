/*
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

#include <iostream>
#include <cstdlib>
#include <fenv.h>
#include <tcl.h>
#include <sys/utsname.h>
#include "io/Diagnostic.h"
#include "io/Command.h"
#include "version.h"

using IO::Command;

int Tcl_AppInit(Tcl_Interp *interp)
{
	Domains::Global::setTcl(interp);
	if(Tcl_Init(interp) == TCL_ERROR)
		return TCL_ERROR;
		
	//Register commands
	Tcl_CreateCommand(interp, "anneal",     Command::Anneal,     (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "exit",       Command::Exit,       (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "extract",    Command::Extract,    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "cascade",    Command::Cascade,    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "cluster",    Command::Cluster,    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "init",       Command::Init,       (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "insert",     Command::Insert,     (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "license",    Command::License,    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "lowmsg",     Command::LowMsg,     (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "param",      Command::Param,      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "profile",    Command::Profile,    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "report",     Command::Report,     (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "restart",    Command::Restart,    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "save",       Command::Save,       (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateCommand(interp, "test",       Command::Test,       (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
	//start ...
	return TCL_OK;
}

int main(int argc, char *argv[])
{
	if(argc >= 2)
	{
		Domains::Global::setInputFileName(std::string(argv[1]));
		Domains::Global::setLogFileName(std::string(argv[1])+".log");
	}
	else
	{
		Domains::Global::setInputFileName("MMonCa");
		Domains::Global::setLogFileName("MMonCa.log");
	}

	utsname buf;
    uname(&buf);


	LOWMSG("#   #        #   #  ###  #   #      ###  ###    ");
	LOWMSG("## ##        ## ## #   # ##  #     #    #   #      Modular");
	LOWMSG("# # #        # # # #   # # # #     #    #####        MC  ");
	LOWMSG("#   #        #   # #   # #  ##     #    #   #      Simulator");
	LOWMSG("#   # odular #   #  ###  #   # te   ### #   # rlo ");
	LOWMSG("Version: " << VERSION);
	LOWMSG("Compiled on " << __DATE__ << " " __TIME__ << " for " << buf.machine << "-" << buf.sysname);
	LOWMSG("for " << buf.version);
	LOWMSG("        Contact: imartin2@ucam.edu");
	LOWMSG("        https://github.com/imartinbragado/MMonCa");
	LOWMSG("OKMC: (C) 2011-2014 IMDEA Materials Institute.");
	LOWMSG("OKMC: (C) 2015-2018 Ignacio Martin, Ignacio Dopico");
	LOWMSG("LKMC: (C) 2011-2014 IMDEA Materials Institute.");
	LOWMSG("LKMC: (C) 2015-2018 Ignacio Martin, Ignacio Dopico");
	LOWMSG("BCA:  (C) 2014 Universidad de Valladolid.");
#if MMONCA_FELIKS
	LOWMSG("FEM:  (C) 2014 IMDEA Materials Institute and");
	LOWMSG("      (C) 2014 Technical University of Madrid (UPM)");
	LOWMSG("      Module based on FELIKS.");
#endif
	LOWMSG(" For licensing details, write \"license\"");

	if(argc > 2 )
		ERRORMSG("usage: " << argv[0] << " [script]");
#ifdef __linux__
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

	Tcl_FindExecutable(argv[0]);
	Tcl_Main(argc, argv, Tcl_AppInit);

	return 0;
}
