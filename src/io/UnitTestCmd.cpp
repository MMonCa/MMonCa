/*
 * UnitTestCmd.cpp
 *
 *  Created on: Nov 2021
 *
 * Author: ignacio.martin
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

#include "UnitTestCmd.h"
#include "electrostatics/FermiDirac.h"
namespace IO {

UnitTestCmd::UnitTestCmd(Tcl_Interp *interp, int argc, const char *argv[]) : Command(interp, argc, argv)
{

}

int UnitTestCmd::operator ()()
{
	if(!test_FermiDirac())
		return TCL_ERROR;

	return TCL_OK;
}

bool UnitTestCmd::test_FermiDirac()
{
	for (double x=0; x>-1e4; x-=100)
		LOWMSG("mhalf " << x << "->" << Electrostatics::FermiDirac::mhalf(x));
	for (double x=0; x<1e4; x+=100)
		LOWMSG("mhalf " << x << "->" << Electrostatics::FermiDirac::mhalf(x));
	for (double x=0; x>-1e4; x-=100)
		LOWMSG("phalf " << x << "->" << Electrostatics::FermiDirac::phalf(x));
	for (double x=0; x<1e4; x+=100)
		LOWMSG("phalf " << x << "->" << Electrostatics::FermiDirac::phalf(x));
	
	return Electrostatics::FermiDirac::phalf(-800) == 0 && 
		   Electrostatics::FermiDirac::mhalf(-800) == 0;
}

}
