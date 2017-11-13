/*
 * ReportCmd.cpp
 *
 *  Created on: Jun 14, 2012
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

#include "ReportCmd.h"
#include "domains/MCClient.h"

namespace IO {

ReportCmd::ReportCmd(Tcl_Interp *interp, int argc, const char *argv[]) : Command(interp, argc, argv)
{
}

int ReportCmd::operator ()()
{
	if(specified("insertions") || specified("all"))
		Domains::global()->client()->endInsert();
	if(specified("events") || specified("all"))
		Domains::global()->eventReport();
	if(specified("reactions") || specified("all"))
		Domains::global()->reactionReport();
	if(specified("reactions.interface"))
		Domains::global()->reactionInterfaceReport();
	if(specified("defects") || specified("all"))
		Domains::global()->defectReport();
	if(specified("mesh") || specified("all"))
		Domains::global()->meshReport();
	if(specified("domains"))
		Domains::global()->domainReport();

	return TCL_OK;
}

} /* namespace IO */
