/*
 * TestCmd.cpp
 *
 *  Created on: Jul 6, 2011
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

#include "TestCmd.h"
#include <sstream>

namespace IO {

TestCmd::TestCmd(Tcl_Interp *interp, int argc, const char *argv[]) : Command(interp, argc, argv)
{

}

int TestCmd::operator ()()
{
	const std::string tag = (specified("tag")? getString("tag") : "");
	float total=0;
	bool bTrue = !specified("not");
	unsigned nFields=0;
	if(specified("array"))
	{
		std::stringstream ss;
		ss << getString("array");
		float value  = getFloat("value");
		float relerr = getFloat("error");
		float init   = getFloat("init");
		float end    = getFloat("end");

		float x, v, max_err=0;
		ss >> x;
		while(!ss.eof())
		{
			float oldx = x;
			ss >> v >> x;
			if(oldx < init || oldx >= end)
				continue;
			total += v;
			nFields++;
			float rel_err = (value != 0? fabs((v-value)/value) : 1.);
			if(value == 0 && v == 0)
				rel_err = 0;
			LOWMSG(tag << " " << oldx << " " << v << " -> (" << value << " " << rel_err*100 << "% )");
			if(rel_err > max_err)
				max_err = rel_err;
		}
		WARNINGMSG(tag << ": Requested error is " << (bTrue? " < ": " > ") <<
				relerr*100 << "% " << " maximum error is " << max_err*100 << "%");
		std::stringstream out;
		if( (max_err > relerr) == bTrue)
		{
			LOWMSG("Test FAILED... aborting");
			out << tag << " Test failed: " << max_err << (bTrue? " > ": " < ") << relerr;
			Tcl_AppendResult(_pTcl, out.str().c_str(), 0);
			return TCL_ERROR;
		}
		out << total/nFields;
		Tcl_AppendResult(_pTcl, out.str().c_str(), 0);
		LOWMSG("Test PASSED... continuing");
		return TCL_OK;
	}
	if(specified("array.2D"))
	{
		std::stringstream ss;
		ss << getString("array.2D");
		float value  = getFloat("value");
		float relerr = getFloat("error");
		Kernel::Coordinates init   = getCoordinates("init");
		Kernel::Coordinates end    = getCoordinates("end");

		float x, y, v, max_err=0;
		ss >> x >> y;
		while(!ss.eof())
		{
			float oldx = x;
			float oldy = y;
			ss >> v >> x >> y;
			if(oldx < init._x || oldx >= end._x || oldy < init._y || oldy >= end._y)
				continue;
			total += v;
			nFields++;
			float rel_err = (value != 0? fabs((v-value)/value) : 1.);
			if(value == 0 && v == 0)
				rel_err = 0;
			LOWMSG(tag << " " << oldx << " " << oldy << " " << v << " -> (" << value << " " << rel_err*100 << "% )");
			if(rel_err > max_err)
				max_err = rel_err;
		}
		WARNINGMSG(tag << ": Requested error is " << (bTrue? " < ": " > ") <<
				relerr*100 << "% " << " maximum error is " << max_err*100 << "%");
		std::stringstream out;
		if( (max_err > relerr) == bTrue)
		{
			LOWMSG("Test FAILED... aborting");
			out << tag << " Test failed: " << max_err << (bTrue? " > ": " < ") << relerr;
			Tcl_AppendResult(_pTcl, out.str().c_str(), 0);
			return TCL_ERROR;
		}
		out << total/nFields;
		Tcl_AppendResult(_pTcl, out.str().c_str(), 0);
		LOWMSG("Test PASSED... continuing");
		return TCL_OK;
	}
	if(specified("arrays.2D"))
	{
		std::stringstream ss0, ss1;
		std::stringstream out;

		ss0 << getString("one");
		ss1 << getString("two");
		float relerr = getFloat("error");

		float x[2], y[2], v[2], max_err=0;
		ss0 >> x[0] >> y[0];
		ss1 >> x[1] >> y[1];
		while(!ss0.eof())
		{
			if(x[0] != x[1] || y[0] != y[1])
			{
				LOWMSG("Test FAILED... aborting");
				out << tag << " Test failed: Coordinates differ. " << x[0] << ", " << y[0] << " != " << x[1] << ", " << y[1];
				Tcl_AppendResult(_pTcl, out.str().c_str(), 0);
				return TCL_ERROR;
			}

			ss0 >> v[0] >> x[0] >> y[0];
			ss1 >> v[1] >> x[1] >> y[1];
			float rel_err = (v[0] != 0? fabs((v[1]-v[0])/v[0]) : 1.);
			if(v[0] == 0 && v[1] == 0)
				rel_err = 0;
			LOWMSG(tag << " " << x[0] << " " << y[0] << " " << v[0] << " -> (" << v[1] << " " << rel_err*100 << "% )");
			if(rel_err > max_err)
				max_err = rel_err;
		}
		WARNINGMSG(tag << ": Requested error is " << (bTrue? " < ": " > ") <<
				relerr*100 << "% " << " maximum error is " << max_err*100 << "%");
		if( (max_err > relerr) == bTrue)
		{
			LOWMSG("Test FAILED... aborting");
			out << tag << " Test failed: " << max_err << (bTrue? " > ": " < ") << relerr;
			Tcl_AppendResult(_pTcl, out.str().c_str(), 0);
			return TCL_ERROR;
		}
		LOWMSG("Test PASSED... continuing");
		return TCL_OK;
	}

	if(specified("array.3D"))
	{
		std::stringstream ss;
		ss << getString("array.3D");
		float value  = getFloat("value");
		float relerr = getFloat("error");
		Kernel::Coordinates init   = getCoordinates("init");
		Kernel::Coordinates end    = getCoordinates("end");

		float x, y, z, v, max_err=0;
		ss >> x >> y >> z;
		while(!ss.eof())
		{
			float oldx = x;
			float oldy = y;
			float oldz = z;
			ss >> v >> x >> y >> z;
			if(		oldx < init._x || oldx >= end._x ||
					oldy < init._y || oldy >= end._y ||
					oldz < init._z || oldz >= end._z)
				continue;
			total += v;
			nFields++;
			float rel_err = (value != 0? fabs((v-value)/value) : 1.);
			if(value == 0 && v == 0)
				rel_err = 0;
			LOWMSG(tag << " " << oldx << " " << oldy << " " << oldz << " "
					<< v << " -> (" << value << " " << rel_err*100 << "% )");
			if(rel_err > max_err)
				max_err = rel_err;
		}
		WARNINGMSG(tag << ": Requested error is " << (bTrue? " < ": " > ") <<
				relerr*100 << "% " << " maximum error is " << max_err*100 << "%");
		std::stringstream out;
		if( (max_err > relerr) == bTrue)
		{
			LOWMSG("Test FAILED... aborting");
			out << tag << " Test failed: " << max_err << (bTrue? " > ": " < ") << relerr;
			Tcl_AppendResult(_pTcl, out.str().c_str(), 0);
			return TCL_ERROR;
		}
		out << total/nFields;
		Tcl_AppendResult(_pTcl, out.str().c_str(), 0);
		LOWMSG("Test PASSED... continuing");
		return TCL_OK;
	}
	else if(specified("float"))
	{
		float value  = getFloat("value");
		float relerr = getFloat("error");
		float orig   = getFloat("float");
		float rel_err = (value != 0? fabs((orig-value)/value) : 1.);
		if(value == 0 && orig == 0)
			rel_err = 0;
		WARNINGMSG(tag << ": Requested error is " << (bTrue? " < ": " > ")
				<< relerr*100 << "% maximum error is " << rel_err*100 << "%");
		std::stringstream out;
		if( (rel_err > relerr) == bTrue)
		{
			out << tag << " Test failed: " << rel_err << (bTrue? " > ": " < ") << relerr;
			Tcl_AppendResult(_pTcl, out.str().c_str(), 0);
			LOWMSG("Test FAILED... aborting");
			return TCL_ERROR;
		}
		LOWMSG("Test PASSED... continuing");
		return TCL_OK;
	}
	else if(specified("equal"))
	{
		std::string txt1 = getString("one");
		std::string txt2 = getString("two");
		WARNINGMSG(tag << ": Requested test " << reformatText(txt1) << " == " << reformatText(txt2));
		if((txt1 == txt2) == bTrue)
			LOWMSG("Test PASSED... continuing");
		else
			ERRORMSG("Test FAILED... aborting");

		return TCL_OK;
	}
	else if(specified("interval"))
	{
		float begin = getFloat("begin");
		float end = getFloat("end");
		float value = getFloat("value");
		WARNINGMSG(tag << ": Requested value " << value << " is expected to be within the [" << begin
				<< " ," << end << "] interval.");
		if(begin <= value && value <= end)
			LOWMSG("Test PASSED... continuing");
		else
			ERRORMSG("Test FAILED... aborting");
		return TCL_OK;
	}

	Tcl_AppendResult(_pTcl, "options allowed are 'array', 'float', 'equal' and 'interval'", 0);
	return TCL_ERROR;
}

std::string TestCmd::reformatText(const std::string &orig)
{
	std::string toReturn;

	for(size_t sz = 0; sz < orig.size(); sz++)
	{
		switch(orig[sz])
		{
		case '\n':
		case '\t':
			toReturn += " ";
			break;
		default:
			toReturn += orig[sz];
			break;
		}
		if(sz > 10 && orig.size() > 15)
		{
			toReturn += "(...)";
			break;
		}
	}
	return toReturn;
}

}
