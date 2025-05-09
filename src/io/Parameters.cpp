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

#include "Parameters.h"
#include <tcl.h>
#include <cstdlib>

using namespace IO;
using std::map;
using std::string;
using std::pair;
using std::vector;
using std::stringstream;
using Domains::Global;

Parameters::~Parameters()
{
}

Parameters::Parameters(std::istream &is)
{
	unsigned size, charSize;
	is >> size;
	string type, value;
	bool bVal;
	for(unsigned i=0; i<size; ++i)
	{
		//Read key
		is >> charSize;
		char *pBuffer = new char[charSize+2];
		is.read(pBuffer, charSize+1); //because there is a remaining space.
		pBuffer[charSize+1] = 0; //finish the char
		string key(pBuffer+1);
		delete [] pBuffer;
		//read type and value
		is >> type >> charSize;
		pBuffer = new char[charSize+2];
		is.read(pBuffer, charSize+1);
		pBuffer[charSize+1] = 0; //finish the char
		_map[key] = std::make_pair(type, string(pBuffer+1));
		delete [] pBuffer;
	}
	is >> type;
	if(type != "Used")
		ERRORMSG("'Used' != '" << type << "'. Corrupted restart file!");
	is >> size;
	for(unsigned i=0; i<size; ++i)
	{
		is >> type >> bVal;
		_used[type] = bVal;
	}
}

void Parameters::restart(std::ostream &os) const
{
	os << _map.size() << "\n";
	for(mmap::const_iterator it=_map.begin(); it!=_map.end(); ++it)
		os << it->first.size() << " " << it->first << " " << it->second.first << " " << it->second.second.size() << " " << it->second.second << std::endl;
	os << "Used " << _used.size() << " ";
	for(map<string, bool>::const_iterator it=_used.begin(); it!=_used.end(); ++it)
		os << it->first << " " << it->second << " ";
}

float Parameters::getFloat(const string &key) const
{
	stringstream ss;
	ss << get(key, "float");
	float ret;
	ss >> ret;
	if(!ss)
		ERRORMSG("'" << get(key, "float") << "' does not seem to be a proper float");
 	return ret;
}

std::vector<float> Parameters::getFloats(const string &key) const
{
    stringstream text;
    text << get(key, "array<float>");

    std::vector<float> result;
    while(!text.eof())
    {
        float tmp;
        text >> tmp;
        if(!text) {
            ERRORMSG("'" << get(key, "array<float>") << "' does not seem to be a proper array<float>");
        }
        else {} // nothing to do
        result.push_back(tmp);
    }
    return result;
}

int Parameters::getInt(const string &key) const
{
	stringstream ss;
	ss << get(key, "int");
	float ret;
	ss >> ret;
	if(!ss)
		ERRORMSG("'" << get(key, "int") << "' does not seem to be a proper int");
	return int(ret);
}

//how many tells how many divisions per string. i.e., a string can have "substrings"
vector<string> Parameters::getStrings(const string &key, unsigned howMany) const
{	
	stringstream isText;
	isText << get(key, "array<string>");

	vector<string> result;
	string temp;
	while(!isText.eof())
	{
		string add;
		isText >> temp;
		add += temp;
		if(isText.fail())
			break;
		for(unsigned i=1; i<howMany; ++i)
		{
			isText >> temp;
			add += " ";
			add += temp;
		}
		if(!isText)
			ERRORMSG("'" << get(key, "array<string>") << "' does not seem to be a proper array<string>");
		result.push_back(add);
	}
	return result;
}

string Parameters::getString(const string &key) const
{
	return get(key, "string");
}

map<string, bool> Parameters::getBoolMap(const string &key) const
{
	map<string, bool> theMap;
	stringstream isText;
	isText << get(key, "map<string,bool>"); //fill the stream with text
	string k,v;
	while(!isText.eof())
	{
		isText >> k;
		if(isText.eof())
			break;
		isText >> v;
		if(!isText)
			ERRORMSG(key << ": map<string,bool> syntax error in " << isText.str());
		if(v != "true" && v != "false")
			ERRORMSG(key << ": map<string,bool> value should be 'true' or 'false', but not " << v);
		theMap[k] = (v == "true"? true: false);
	}
	return theMap;
}

map<string, int> Parameters::getIntMap(const string &key) const
{
   map<string, int> theMap;
   stringstream isText;
   isText << get(key, "map<string,int>"); //fill the stream with text
   string k;
   int v;
   while(!isText.eof())
   {
	 isText >> k;
	 if(isText.eof())
			 break;
	 isText >> v;
	 if(!isText)
	   ERRORMSG(key << ": map<string,int> Syntax error in \"" << isText.str() << "\"");
	 theMap[k] = v;
   }
   return theMap;
}

map<string, IO::ArrheniusAlloys> Parameters::getArrheniusAlloysMap(const string &key) const
{
	map<string, IO::ArrheniusAlloys> theMap;
	stringstream isText;
	isText << get(key, "map<string,arrhenius>"); //fill the stream with text
	if (isText.str().find("%") != isText.str().npos)
	{
		string k;
		isText >> k;
		vector<float> vpref, vener, vcomp;
		while (isText)
		{
			string s_percentage;        // Stores composition for alloy with "%" symbol
			float p, e, per;            // Temporally stores each element for later append to vector
			isText >> s_percentage >> p >> e;
			if (s_percentage.size () > 0)  s_percentage.resize (s_percentage.size () - 1);  // Check if non-empty string and remove the "%" symbol
			stringstream ss_per(s_percentage);
			ss_per >> per;
			vcomp.push_back(per);
			vpref.push_back(p);
			vener.push_back(e);
		};
		vpref.pop_back();   // Pop back last element which is duplicated
		vener.pop_back();
		vcomp.pop_back();

		ArrheniusAlloys arr(vpref, vener, vcomp);
		theMap[k] = arr;
	}
	else
	{
		string k;
		float pref, ener;
		while(!isText.eof())
		{
			isText >> k;
			if(isText.eof())
				break;
			isText >> pref >> ener;
			if(!isText)
				ERRORMSG(key << ": map<string,arrhenius> syntax error in " << isText.str());
			theMap[k] = IO::ArrheniusAlloys(pref, ener);
		}
	}
	return theMap;
}

map<string, string> Parameters::getStringMap(const string &key) const
{
	map<string, string> theMap;
	stringstream isText;
	isText << get(key, "map<string,string>"); //fill the stream with text
	string k,v;
	while(!isText.eof())
	{
	  isText >> k;
	  if(isText.eof())
		  break;
	  isText >> v;
	  if(!isText)
	    ERRORMSG(key << ": map<string,string> syntax error in " << isText.str());
	  theMap[k] = v;
	}
	return theMap;
}
    
map<string, float> Parameters::getFloatMap(const string &key) const
{
	map<string, float> theMap;
	stringstream isText;
	isText << get(key, "map<string,float>"); //fill the stream with text
	string k;
	float v;
	while(!isText.eof())
	{
	  isText >> k;
	  if(isText.eof())
		  break;
	  isText >> v;
	  if(!isText)
	    ERRORMSG(key << ": map<string,float> Syntax error in \"" << isText.str() << "\"");
	  theMap[k] = v;
	}
	return theMap;
}

//sets only 1 key.
array<string, string> Parameters::getArray(const string &key) const
{
	array<string, string> theArray;
	stringstream isText;
	isText << get(key, "array<string,string>"); //fill the stream with text
	string k,v;
	while(!isText.eof())
	{
		isText >> k;
		if(isText.eof())
		  break;
		isText >> v;
		if(!isText)
			ERRORMSG(key << ": array<string,string> Syntax error in \"" << isText.str() << "\"");
		theArray.push_back(pair<string, string>(k,v));
	}
	return theArray;
}

bool Parameters::getBool(const string &key) const
{
	string value = get(key, "bool");
 
	bool b=false;
	switch(value[0])
	{
	case '0':case 'f':case 'F':case 'n':case 'N':
	  b = false;
	  break;
	case '1':case 't':case 'T':case 'y':case 'Y':
	  b = true;
	  break;
	default:
	  ERRORMSG("bool: syntax error in " << value << " Valid arguments are 0,1, true, false, yes, not");
	  break;
	}
	return b;
}

Kernel::Coordinates Parameters::getCoordinates(const string &key) const
{
	string value = get(key, "coordinates");
	stringstream ss;
	ss << value;
	float x,y,z;
	ss >> x >> y >> z;
	if(!ss)
		ERRORMSG("'" << get(key, "coordinates") << "' does not seem to be a proper coordinate (3D vector)");
	return Kernel::Coordinates(x, y, z);
}

Polynomial Parameters::getPolynomial(const string &key) const
{
	stringstream isText;
	isText << get(key, "map<string,float>");
	vector <vector <double> > vPoly;
	vector <double> poly, bounds, center;
	unsigned currentPol = 1;
	while(!isText.eof())
	{
		string pol, textNum;
		unsigned polNum;
		double val, xo;
		isText >> pol;
		if(pol == "xi")
		{
			double bound;
			isText >> bound;
			if(!isText)
				ERRORMSG("'" << get(key, "map<string,float>") << "' does not seem to be a proper polynomial");
			bounds.push_back(bound);
			continue;
		}
		if(pol == "xo")
		{
			isText >> xo;
			if(!isText)
				ERRORMSG("'" << get(key, "map<string,float>") << "' does not seem to be a proper polynomial");
			center.push_back(xo);
			continue;
		}
		isText >> val;
		textNum = pol[4];
		polNum = std::atoi(textNum.c_str());
		if(polNum == currentPol)
		{
			poly.push_back(val);
		}
		else
		{
			vPoly.push_back(poly);
			poly.clear();
			poly.push_back(val);
			currentPol = polNum;
		}
	}
	if(poly.size())
		vPoly.push_back(poly);
	return Polynomial(vPoly, bounds, center);
}

ArrheniusAlloys Parameters::getArrheniusAlloys(const string &key) const
{
	return toArrheniusAlloys(key, get(key, "arrhenius"), AA_FULL);
}

ArrheniusAlloys Parameters::toArrheniusAlloys(const string &key, const string &value, AATYPE type)
{
	stringstream isText;
	isText << value;
	size_t found = value.find("%");
	if (found != string::npos)
	{
		vector<float> vcomp, vpref, vener;
		while (!isText.eof())
		{
			string s_percentage;        // Stores composition for alloy with "%" symbol
			isText >> s_percentage;
			if(isText.eof())
				break;
			if(!isText)
				ERRORMSG(key << ".'" << value << "' does not seem to be a proper arrheniusAlloy.");
			float p=1, e=0, per;            // Temporally stores each element for later append to vector
			if(type == AA_FULL)
				isText >> p >> e;
			else if(type == AA_ENERGY)
				isText >> e;
			else
				isText >> p;
			if(!isText)
				ERRORMSG(key << ".'"  << value << "' does not seem to be a proper arrheniusAlloy.");
			if (s_percentage.size () > 0)  s_percentage.resize (s_percentage.size () - 1);  // Check if non-empty string and remove the "%" symbol
			stringstream ss_per(s_percentage);
			ss_per >> per;
			vcomp.push_back(per / 100.); // From now it will be treated as alloy fraction
			vpref.push_back(p);
			vener.push_back(e);
		};
		return ArrheniusAlloys(vpref, vener, vcomp);
	}
	else                         // if not alloy only one pair of data is given and composition is set to 0%
	{
		float p = 1, e = 0;
		if(type == AA_FULL)
			isText >> p >> e;
		else if(type == AA_ENERGY)
			isText >> e;
		else
			isText >> p;
		if(!isText)
			ERRORMSG(key << ".'"  << value << "' does not seem to be a proper arrheniusAlloy");
		return ArrheniusAlloys(p, e);
	}
}

const string & Parameters::get(const string &key, const string &type) const
{
	mmap::const_iterator iBegin = _map.lower_bound(key), iEnd = _map.upper_bound(key);
	map<string, bool>::iterator it2 = _used.find(key);
	for(mmap::const_iterator it = iBegin; it!=iEnd; ++it)
	{
		if(it->second.first == type || it->second.first.empty())
		{
			if(!it->second.second.size())
				ERRORMSG("Parameter '" << key << "' does not seem to have a value (it exists, but is empty)");
			assert(it2 != _used.end());
			it2->second = true;
			return it->second.second;
		}
	}
	if(_map.find(key) != _map.end())
		ERRORMSG(type << " " << key << " not found");
	else
		ERRORMSG(key << " not found");
	return iBegin->first;
}

bool Parameters::specified(const string &key, const string &type) const
{
	mmap::const_iterator iBegin = _map.lower_bound(key);
	mmap::const_iterator iEnd   = _map.upper_bound(key);
	for(mmap::const_iterator it=iBegin; it!=iEnd; ++it)
		if(it->second.first == type)
		{
			_used[key] = true;
			return true;
		}
	return false;
}

bool Parameters::specified(const string &key) const
{
	if(_map.find(key) != _map.end())
	{
		_used[key] = true;
		return true;
	}
	return false;
}

//takes a string associated with the key, and sources it into TCL as a procedure.
void Parameters::loadProcedure(Tcl_Interp *pTcl, const string &key, unsigned argc) const
{
	stringstream cmd;
	cmd << "proc " << key << " { ";
	for(unsigned i=0; i<argc; ++i)
		cmd << "arg" << i << " ";
	cmd << " } { " << get(key, "proc") << " }";
	if(Tcl_EvalEx(pTcl,  cmd.str().c_str(), -1, 0) != TCL_OK)
		ERRORMSG("Procedure '" << cmd.str() << "' failed: " << Tcl_GetStringResult(pTcl));
	else
		MEDMSG("Procedure '" << cmd.str() << "' sourced");
}


std::map<std::string, Arrhenius> Parameters::getArrheniusProc(Tcl_Interp *pTcl, const std::string &key) const
{
	map<string, Arrhenius> temp;
	std::stringstream buffer;
	try
	{
		std::stringstream cmd;
		cmd << key;
		if(Tcl_EvalEx(pTcl,  cmd.str().c_str(), -1, 0) != TCL_OK)
			ERRORMSG("(ArrheniusProc) Execution of the script '" << cmd.str() << "' failed: " << Tcl_GetStringResult(pTcl));
		buffer << Tcl_GetStringResult(pTcl);
	}
	catch(...)
	{
		ERRORMSG("Exception caught while executing " << key << ": " << Tcl_GetStringResult(pTcl));
	}
	Tcl_ResetResult(pTcl);

	string k;
	float v1, v2;
	while(!buffer.eof())
	{
		buffer >> k;
		if(buffer.eof())
			break;
		buffer >> v1 >> v2;
		if(!buffer)
			ERRORMSG(key << " getArrheniusProc Syntax error in \"" << buffer.str() << "\"");
		temp[k] = Arrhenius(v1, v2);
	}
	return temp;
}

map<string,ArrheniusAlloys> Parameters::getArrheniusAlloysProc(Tcl_Interp *pTcl, const string &key, AATYPE type) const
{
	map<string,ArrheniusAlloys> temp;
	std::stringstream buffer;

	try
	{
		std::stringstream cmd;
		cmd << key;
		if(Tcl_EvalEx(pTcl,  cmd.str().c_str(), -1, 0) != TCL_OK)
			ERRORMSG("(MapProc) Execution of the script '" << cmd.str() << "' failed: " << Tcl_GetStringResult(pTcl));
		buffer << Tcl_GetStringResult(pTcl);
	}
	catch(...)
	{
		ERRORMSG("Exception caught while executing " << key << ": " << Tcl_GetStringResult(pTcl));
	}
	Tcl_ResetResult(pTcl);

	string k;
	while(!buffer.eof())
	{
		buffer >> k;
		if(buffer.eof())
			break;
		string txt;
		try
		{
			txt = getTextInBrackets(buffer);
		}
		catch(const char *txt)
		{
			ERRORMSG(txt << " in " << key << " reading " << buffer.str());
		}
		if(buffer.fail())
			ERRORMSG(key << " getArrheniusProc Syntax error in \"" << buffer.str() << "\"");
		temp[k] = toArrheniusAlloys(key, txt, type);
	}

	return temp;
}

//input / output from script
map<string,float> Parameters::getFloatProc(Tcl_Interp *pTcl, const string &key) const
{
	std::stringstream buffer;
	map<string,float> temp;

	try
	{
		if(Tcl_EvalEx(pTcl,  key.c_str(), -1, 0) != TCL_OK)
			ERRORMSG("(MapProc) Execution of the script '" << key << "' failed: " << Tcl_GetStringResult(pTcl));
		buffer << Tcl_GetStringResult(pTcl);
	}
	catch(...)
	{
		ERRORMSG("Exception caught while executing " << key);
	}
	Tcl_ResetResult(pTcl);

	string k;
	float v;
	while(!buffer.eof())
	{
	  buffer >> k;
	  if(buffer.eof())
		  break;
	  buffer >> v;
	  if(!buffer)
		ERRORMSG(key << " getFloatProc Syntax error in \"" << buffer.str() << "\"");
	  temp[k] = v;
	}

	return temp;
}

float Parameters::getFloatProc(Tcl_Interp *pTcl, const string &key, const string &idx) const
{
	std::stringstream buffer;
	try
	{
		std::stringstream cmd;
		cmd << key << " " << idx;
		if(Tcl_EvalEx(pTcl,  cmd.str().c_str(), -1, 0) != TCL_OK)
			ERRORMSG("(Float) Execution of the script '" << cmd.str() << "' failed: "  << Tcl_GetStringResult(pTcl));
		buffer << Tcl_GetStringResult(pTcl);
	}
	catch(...)
	{
		ERRORMSG("Exception caught while executing " << key << " " << idx);
	}
	float res;
	buffer >> res;
	Tcl_ResetResult(pTcl);

	return res;
}

//args is a variable number of arguments. It is just passed to the function.
float Parameters::getFloatProcArgs(Tcl_Interp *pTcl, const string &key, const string &args) const
{
	std::stringstream buffer;
	try
	{
		std::stringstream cmd;
		cmd << key << " " << args;
		if(Tcl_EvalEx(pTcl,  cmd.str().c_str(), -1, 0) != TCL_OK)
			ERRORMSG("(FloatArgs) Execution of the script '" << cmd.str() << "' failed: "  << Tcl_GetStringResult(pTcl));
		buffer << Tcl_GetStringResult(pTcl);
	}
	catch(...)
	{
		ERRORMSG("Exception caught while executing " << key << " " << args);
	}
	float res;
	buffer >> res;
	Tcl_ResetResult(pTcl);

	return res;
}

 //input / output from script
void Parameters::setCommand(const string &type, const string &key, const string &value, const string &index)
{
	if(!specified(key, type)) //new parameter
	{
		if (index.size())
			ERRORMSG("Cannot use index when setting a new parameter.");
		insert(key, type, value);
		_used[key] = false;
	}
	else
	{
		// _map: key -> (type, value)
		mmap::iterator it=_map.end(), iBegin=_map.lower_bound(key), iEnd=_map.upper_bound(key);
		for(it=iBegin; it != iEnd; ++it)
			if(it->second.first == type)
			{
				if(index.empty())
				{
					it->second.second = value;
					if(type == "proc" && value.find('"') == string::npos)
						WARNINGMSG("Quotes \" should be escaped in procedure " << key << " " << value);
				}
				else
				{
					bool bAdd = true;
					if (type == "arrhenius")
					{
						stringstream isText;
						isText << it->second.second;
						LOWMSG("setParam isText: " << isText.str());
						it->second.second = "";
						string k, v, w;

						isText >> k >> v >> w;
						if(k == index)
						{
							it->second.second += index + " " + value + " ";
							bAdd = false;
						}
						else
							it->second.second += k + " " + v + " " + w + " ";
					}
					else if(type == "map<string,string>" || type == "map<string,float>" ||
					   type == "map<string,bool>"   || type == "map<string,int>")
					{
						stringstream isText;
						isText << it->second.second;
						it->second.second = "";
						string k,v;

						while(!isText.eof())
						{
							isText >> k;
							if(isText.eof())
							  break;
							isText >> v;
							if(isText.fail())
								ERRORMSG(key << ": setCommand Syntax error in \"" << isText.str() << "\"");
							if(k == index)
							{
								it->second.second += k + " " + value + " ";
								bAdd = false;
							}
							else
								it->second.second += k + " " + v + " ";
						}
					}
					else if (type == "array<string,string>")
					//if the key is repeated, takes it from the middle and adds it at the end.
					{
						stringstream isText;
						isText << it->second.second;
						it->second.second = "";
						string k, v;

						while(!isText.eof())
						{
							isText >> k;
							if(isText.eof())
							  break;
							isText >> v;
							if(isText.fail())
								ERRORMSG(key << ": setCommand Syntax error in \"" << isText.str() << "\"");
							if(k != index)
							    it->second.second += k + " " + v + " ";
						}
					}
					else if (type == "array<string>")
					{
						ERRORMSG(key << ": array<string> does not accept the index argument.");
					}
					else if(type == "map<string,arrhenius>")
					{
						stringstream isText;
						isText << it->second.second;
						it->second.second = "";
						string k,v1,v2;
						if (isText.str().find("%") != isText.str().npos)
						{
							string comp;
							while(!isText.eof())
							{
								isText >> k;
								if(isText.eof())
								  break;
								isText >> comp >> v1 >> v2;
								if(isText.fail())
									ERRORMSG(key << ": setCommand Syntax error in \"" << isText.str() << "\"");
								if(k == index)
								{
									it->second.second += k + " " + value + " ";
									bAdd = false;
								}
								else
									it->second.second += k + " " + comp + " " + v1 + " " + v2 + " ";
							}
						}
						else
						{
						    while(!isText.eof())
						    {
						    	isText >> k;
						    	if(isText.eof())
						    		break;
						    	isText >> v1 >> v2;
						    	if(isText.fail())
						    		ERRORMSG(key << ": setCommand Syntax error in \"" << isText.str() << "\"");
						    	if(k == index)
						    	{
						    		it->second.second += k + " " + value + " ";
						    		bAdd = false;
						    	}
						    	else
						    		it->second.second += k + " " + v1 + " " + v2 + " ";
						    }
						}
					}
					else
						ERRORMSG("Type " << type << " does not accept parameter index");
					if(bAdd) //the key was not find, we have to add the new element
						it->second.second += index + " " + value + " ";
				}
			}
	}
}

void Parameters::addCommand(const string &type, const string &key, const string &value, const string &index)
{
	if(type == "array<string,string>")
	{
		if(index.empty())
			ERRORMSG("add needs a valid index.");

		mmap::iterator it=_map.end(), iBegin=_map.lower_bound(key), iEnd=_map.upper_bound(key);
		for(it=iBegin; it != iEnd; ++it)
			if(it->second.first == type)
			{
				it->second.second += index + " " + value + " ";
				break;
			}
	}
	else if(type == "array<string>")
	{
		mmap::iterator it=_map.end(), iBegin=_map.lower_bound(key), iEnd=_map.upper_bound(key);
		for(it=iBegin; it != iEnd; ++it)
			if(it->second.first == type)
			{
				it->second.second += value;
				break;
			}
	}
	else
		ERRORMSG("add cannot be used for types different than array");
}

void Parameters::unsetCommand(const string &type, const string &key, const string &index)
{
	// so, there is index.
	mmap::iterator it=_map.end(), iBegin=_map.lower_bound(key), iEnd=_map.upper_bound(key);
	for(it=iBegin; it != iEnd; ++it)
		if(it->second.first == type)
		{
			if(index.empty())
			{
				_map.erase(it);
				return;
			}
			else
			{
				bool bFound = false;
				if(type == "map<string,string>" || type == "map<string,float>" ||
				   type == "map<string,bool>"   || type == "map<string,int>" ||
				   type == "array<string,string>")
				{
					stringstream isText;
					isText << it->second.second;
					it->second.second = "";
					string k,v;

					while(!isText.eof())
					{
						isText >> k;
						if(isText.eof())
						  break;
						isText >> v;
						if(isText.fail())
							ERRORMSG(key << ": unsetCommand Syntax error in \"" << isText.str() << "\"");
						if(k == index)
							bFound = true;
						else
							it->second.second += k + " " + v + " ";
					}
				}
				else if(type == "map<string,arrhenius>")
				{
					stringstream isText;
					isText << it->second.second;
					it->second.second = "";
					string k,v1,v2;
					while(!isText.eof())
					{
						isText >> k;
						if(isText.eof())
						  break;
						isText >> v1 >> v2;
						if(isText.fail())
							ERRORMSG(key << ": setCommand Syntax error in \"" << isText.str() << "\"");
						if(k == index)
							bFound = true;
						else
							it->second.second += k + " " + v1 + " " + v2 + " ";
					}
				}
				else
					ERRORMSG("Type " << type << " does not accept parameter index");
				if(!bFound)
					WARNINGMSG("unset: Could not find " << key << " with type " << type << " and index " << index);
				return;
			}
		}
	ERRORMSG("unset: Could not find " << key << " with type " << type);
}

string Parameters::getCommand(const string &type, const string &key, const string &index) const
{
	stringstream out;

	if(index.size())
	{
		if(type == "array<string,string>")
		{
			bool b = false;
			array<string, string> theArray = getArray(key);
			for(array<string, string>::iterator it = theArray.begin(); it != theArray.end(); ++it)
				if(it->first == index)
				{
					out << (b? " ":"") << it->second;
					b= true;
				}
		}
		else if(type == "map<string,string>")
		{
			map<string, string> theMap = getStringMap(key);
			map<string, string>::iterator it = theMap.find(index);
			if(it != theMap.end())
			out << it->second;
		}
		else if(type == "map<string,float>")
		{
			map<string, float> theMap = getFloatMap(key);
			map<string, float>::iterator it = theMap.find(index);
			if(it != theMap.end())
				out << it->second;
		}
		else if(type == "map<string,int>")
		{
			map<string, int> theMap = getIntMap(key);
			map<string, int>::iterator it = theMap.find(index);
			if(it != theMap.end())
				out << it->second;
		}
		else if(type == "map<string,bool>")
		{
			map<string, bool> theMap = getBoolMap(key);
			map<string, bool>::iterator it = theMap.find(index);
			if(it != theMap.end())
				out << it->second;
		}
		else if(type == "map<string,arrhenius>")
		{
			map<string, IO::ArrheniusAlloys> theMap = getArrheniusAlloysMap(key);
			map<string, IO::ArrheniusAlloys>::iterator it = theMap.find(index);
			if(it != theMap.end())
				out << it->second._v[0]._pref << " " << it->second._v[0]._ener;
		}
		else
			ERRORMSG("Type " << type << " does not accept parameter index");
		return out.str().c_str();
	}
	else
		return get(key, type);
}

void Parameters::insert(const std::string &key, const std::string &type, const std::string &value)
{
	_map[key] = make_pair(type, value);
}

string Parameters::getTextInBrackets(std::stringstream &ss)
{
	unsigned howMany = 1;
	string txt;
	ss >> txt;
	if(ss.fail())
		throw "Syntax error at the beginning";
	if(txt != "<")
		throw (string("Expected '<' not found: " + txt + " found instead")).c_str();
	string ret;
	while(howMany)
	{
		ss >> txt;
		if(txt == "<")
			howMany++;
		else if(txt == ">")
			howMany--;
		else
		{
			ret += " ";
			ret += txt;
		}
		if(ss.fail())
			throw "Syntax error. End of buffer while waiting for '>'";
	}
	return ret;
}
