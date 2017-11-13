/* Finite Element Method Module
 *
 * Author: ignacio.romero@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain,
 *      and       Technical University of Madrid (UPM), Madrid, Spain
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
/*
 * usercommand.cpp
 * 
 * routines that prompt the user for commands and recognize them
 * i. romero 
 * septermber 2001, transformed to C++ in dec 2005
 */

#include "Io/usercommand.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include <cstdio>

#define COMMENTCHAR      35		// ascii for #
#define CONTINUECHAR     92     // ascii for backslash
#define UDEBUG            0

using namespace std;

static void   removeSpaces(string& str);


usercommand :: usercommand() :
	_numarg(0.0),
	_keyword(""),
	_option(""),
	_hasStringArgument(false),
	_hasNumberArgument(false)
{	
}


usercommand :: usercommand(const usercommand &uc) :
	_keyword(uc.keyword()),
	_option(uc.option()),
	_numarg(uc._numarg),
	_hasNumberArgument(uc._hasNumberArgument),
	_hasStringArgument(uc._hasStringArgument)
{
}




usercommand :: usercommand(const std::string& keyword,
                           const std::string& option,
                           const double value)
:
    _keyword(keyword),
    _option(option),
    _numarg(value),
    _hasStringArgument(false),
    _hasNumberArgument(true)
{        
}


usercommand :: usercommand(const string& inputsentence) :
	_numarg(0.0),
	_keyword(""),
	_option(""),
	_hasStringArgument(false),
	_hasNumberArgument(false)
{    
	int start=0, end=0;
	char c;
	_numarg = 0.0;
    string sentence = inputsentence;
    removeSpaces(sentence);

	if (UDEBUG == 1) cout << endl << "usercommand constructor from sentence(" << sentence.length() << "):" << sentence << endl;
    
    // advance spaces
	while ( (c=sentence[start]) == ' ') start++;
		
    // get the keyword
    end = start;
    while (isalnum(c=sentence[end]) || c == '_' || c == '.') end++;
	_keyword = sentence.substr(start, end-start);
	if (UDEBUG == 1) cout << "     keyword:" << _keyword << "\n";
	
	// advance until option
	if (end == sentence.length()) return;
	start = end;
	while ( (c=sentence[start]) == ' ' || c == '=') start++;
	if (end == sentence.length()) return;
	
	// get the option/argument. This is a bit more involved. It can be a word->the option, 
	// or a number->the argument, or a string delimited by " or '
	// we assume it is an option if the first character is not a digit
	end = start;
	char firstc = sentence[start];
	
	//the string with double quotes
	if (firstc == 34 )
	{
		end++;
		while(sentence[end] != 34) end++;
		_option = sentence.substr(start+1, end-start-1);
		_hasStringArgument = true;
		if (UDEBUG == 1) cout << "     option:" << _option << endl;
	}

	// the string with single quote
	else if (firstc == 39 )
	{
		end++;
		while(sentence[end] != 39) end++;
		_option = sentence.substr(start+1, end-start-1);
		_hasStringArgument = true;
		if (UDEBUG == 1) cout << "     option:" << _option << endl;
	}
	
	
	// option
	else if ( !( isdigit(firstc) || firstc == '-' || firstc == '+') )
	{
		while (isalnum(c=sentence[end]) || c == '.' || c == '/' || c == '_') end++;
		_option = sentence.substr(start, end-start);
		if (UDEBUG == 1) cout << "     option:" << _option << endl;
		_hasStringArgument = true;
	}
	
    // scan numerical arguments
	if (end == sentence.length()) return;
	
	istringstream rest(sentence.substr(end, sentence.length()-end));
	rest >> _numarg;
	_hasNumberArgument = true;
}



usercommand :: ~usercommand()
{
}



const std::string&   usercommand :: keyword() const
{
    return _keyword;
}



const std::string&   usercommand :: option() const
{
    return _option;
}



// this function returns the number obtained after evaluating the usercommand
// so, eg,     radius = 3.14   , will give uc.value() ---> 3.14
//             radius = "_pi"   , wile give uc.value() ---> 3.14159....
double usercommand :: value() const
{
	double r(0.0);
	
	if (_hasStringArgument)
	{
	}
	else
		r = _numarg;
		
	return r;	
}



int usercommand :: integerValue() const
{
    return static_cast<int>(value());
}


/*
usercommand& usercommand :: operator=(const usercommand &uc)
{
	_keyword = uc.keyword();
	_option  = uc.option();
	_numarg  = uc._numarg;
	
	_hasStringArgument = uc._hasStringArgument;
	_hasNumberArgument = uc._hasNumberArgument;
	
	return *this;
}
*/




void usercommand :: print(std::ostream& of) const
{
	of << "     usercommand = keyword |" << keyword() << "|, option |" << option()
		<< "|, argument : " << _numarg << "\n";
}




// an empty vector of commands
commandLine :: commandLine() 
{
}



// the commandLine contains a vector of usercommand pointers that have been allocated in the
// constructor commandLine(string), and we must free their memory
commandLine :: ~commandLine()
{
}



commandLine :: commandLine(const std::string& textline)
{
    readFromString(textline);
}




bool commandLine :: hasKeyword(const std::string& word) const
{
    bool ret = false;
    
    std::vector<usercommand>::const_iterator iter = begin();
	while( iter != end() )
	{
        if ( iter->keyword() == word)
        {
            ret = true;
            break;
        }
		++iter;
	}
    
    
    return ret;
}





std::string commandLine :: optionFor(const std::string& keyword) const
{
    std::string theOption;

    std::vector<usercommand>::const_iterator iter = begin();
	while( iter != end() )
	{
        if ( iter->keyword() == keyword)
        {
            theOption = iter->option();
            break;
        }
		++iter;
	}

    return theOption;
}



int  commandLine :: readFromStream(ifstream &in)
{
	string line;

	if (!scanCommandLine(in, line)) return 0;
	
	// fill the commandLine with new commands
	readFromString(line);
	return 1;
}




/* this function reads one line of a file and puts its content inside a commandLine
   structure. The line must be a series of commands separated by commas.
   it returns the number of commands read.
*/
size_t  commandLine :: readFromString(const string  &line)
{
	if (UDEBUG == 1) cout << "inside readFromString |" << line << "|" << endl;
	
	// empty all the commands in the commandLine
	this->clear();
	
    // check for empty line. If empty, return with no error
	if (line.empty()) return 0;
	
    // read each sentence of the command line and create a command from each
	int start = 0;
	int end   = 0;
	char c;
    while ( (c=line[end]) != '\0')
    {
        if (c == ',')
        {
			usercommand uc(line.substr(start, end-start));
			this->push_back(uc);
			start = ++end;
        }
        
        else if ( c == '"' )
        {
            end++;
            while ( line[end] != '"') end++;
            end++;
        }
        
        else end++;
    }
	usercommand uc(line.substr(start, end-start));
	this->push_back(uc);
	
	return size();
}




void commandLine :: print(std::ostream& of) const
{
	of  << "\n\n CommandLine:";
	for (int k=0; k< size(); k++)  (*this)[k].print(of);
}







/* Scans a line of the input file in the FELIKS format, that is, taking out
 * comments and appending lines with the '\' symbol at the end.
 * returns 0 if end of file reached. Otherwise 1.
 * note that the return value could be 1 but the line empty.
 */
bool commandLine :: scanCommandLine(istream &infile, string &line)
{
    // erase whatever line contained
	line.clear();
	
	// get one line from the stream
	string tmp;
	if (UDEBUG == 1) cout << "inside scanCommandLine" << endl;
    if (!getline(infile, tmp)) return false;
    if (UDEBUG == 1) cout << "     full line:" << tmp << endl;

	// transfer the contents of the line 'tmp' to a stringstream
	// count spaces at beginning
	int start=0;
	while (tmp[start] == ' ' || tmp[start] == '\t') start++;
	if (start == tmp.size()) return true;

	// count spaces at the end
	int end = tmp.size()-1;
	while(tmp[end] == ' ') end--;

	ostringstream cleanline;
    int n = start;
	char ch;
	while ( n <= end && (ch=tmp[n] ) != COMMENTCHAR)
	{
		cleanline << ch;
		n++;
	}

	// transfer the stringstream to the string line
	line = cleanline.str();

    // if last character of the line is CONTINUECHAR means command continues in the next line
	// work recursively
	if (!line.empty() && line[line.size()-1] == CONTINUECHAR) 
    {
		line.erase(line.size()-1, 1);
		if (!scanCommandLine(infile, tmp)) return false;
		line = line + tmp;
    }

	if (UDEBUG == 1) cout << "     This is the line read from the stream:" << line << endl;
	return true;
}




double commandLine :: valueFor(const std::string& keyword) const
{
    double theValue=0.0;
    
    std::vector<usercommand>::const_iterator iter = begin();
	while( iter != end() )
	{
        if ( iter->keyword() == keyword)
        {
            theValue = iter->value();
            break;
        }
		++iter;
	}
    
    return theValue;
}





// output a commandLine
ostream&  operator<<(ostream& os, const commandLine& cl)
{
	vector<usercommand>::const_iterator iter = cl.begin();
	os << "\n";
	while( iter != cl.end() )
	{
		(*iter).print(os);
		++iter;
	}
	return os;
}



// inputs a commandLine from a stream
// usage cin >> cline
//

istream&  operator>>(istream &is, commandLine &cl)
{
	string line;
    
	// The first thing to do is to empty whatever the commandLine held
	// before entering this function, for it is typical to reuse clines
	// in loops within parsers.
	cl.clear();
	bool ok = commandLine::scanCommandLine(is, line);
	
	if (ok) cl.readFromString(line);	
	return is;
}




static void   removeSpaces(string& str)
{
    string old = str;
    str.clear();
    for (int a=0; a<old.size(); a++)
        if ( old[a] != ' ') str.append(1, old[a]);
}


