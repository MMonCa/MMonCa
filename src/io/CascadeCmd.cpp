/*
 * CascadeCmd.cpp
 *
 *  Created on: May 23, 2011
 *      Author: ignacio.martin@imdea.org
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

#include "CascadeCmd.h"
#include "cascade/FromFile.h"
#include "io/ParameterManager.h"
#include "kernel/Domain.h"
#include "kernel/EventLog.h"
#include "kernel/ReactionLog.h"

#include "lkmc/Lattice.h"
#include "io/SaveCmd.h"
#include "kernel/RateManager.h"

using std::string;

namespace IO {


CascadeCmd::CascadeCmd(Tcl_Interp *p, int argc, const char *argv[]) : Command(p, argc, argv)
{

}

int CascadeCmd::operator()()
{
	const float fluence   =  getFloat("fluence");
	const float flux      =  specified("flux") ? getFloat("flux") : 0; //0 means "instantaneous"
	const bool  periodic  =  specified("periodic");
	const float temp      =  specified("temp") ? getFloat("temp") : -1; //-1 means previous
	const bool  voluminic =  specified("voluminic");
	const bool  bReact    = !specified("do.not.react");
	const bool  bCorrectX =  specified("correct.for.surface");

	if(specified("file"))
	{
		string file   = getString("file");
		string format = specified("format")?  getString("format")  : "A:B:C:D";
		string defect = specified("defects")? getString("defects") : "";
		bool bDisplace = !specified("do.not.shift");

		 std::vector<std::string> txts;
		 IO::ParameterManager::getTokens(format, ':', txts);
		 if(txts.size() != 4)
			 ERRORMSG("implant: format should be a string with 4 fields separated by ':'");
		 if(txts[0].size() != 1 || !std::isupper(txts[0][0]))
			 ERRORMSG("First format column must be a single uppercase character...");
		 Cascade::FromFile fromFile;
		 fromFile(file, txts, fluence, flux, defect, bDisplace, periodic, bReact, temp, voluminic, bCorrectX);
	}
	return TCL_OK;
}

}
