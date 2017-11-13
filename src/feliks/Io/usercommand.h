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
 * usercommand.h
 * 
 * classes to deal with input files and parse their information.
 *
 *
 * A comandlist is a list of usercommands separated by commas and ended with an end-of-line.
 * A usercommand is a keyword (=) option/ (numeric argument) (numeric argument). Lines
 * that start with '#' are considered as comments
 *
 * For example this is a commandlist:
 *
 *           # element definition
 *           eltype, label = 3, type = bar3d, thickness = 0.33, crosssection = circle
 *
 * Commandlists can span more than one line if the last (non-white) character of the line
 * is the backslash '\'. For example, this is another commandlist:
 *
 *          eltype, label = 5, type = frod, material = 2, area = 4.0, i11  = 5e-3, \
 *              i22 = 3.3e-3, j = 2e-3
 *
 * Usage:
 *
 *   To use objects of this class, first one reads the file that contains the commands:
 *          ifstream in("analysis.feliks");
 *			while ( in >> cl) 
 *               {
 *                  // process the commandlist, doing for example:	 
 *					cout << cl;
 *               }
 *
 *   Once the commandlist has been read, to use each usercommand one proceeds usually as
 *
 for (int j=0; j< cl.size(); j++)
 {
    usercommand &uc = cl.command(j);
    if(     uc.keyword() == "area") area     = uc.value();
    else if(uc.keyword() == "i11")  inertia1 = uc.value();
    else if(uc.keyword() == "i22")  inertia2 = uc.value();
 ...
 }
 *
 *
 *    
 *
 */

#ifndef _usercommand_h
#define _usercommand_h

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

/*
 * data structure.
 * main command, a character option and a numerical argument
 */

class usercommand
{
	
public:
                            usercommand();
                            usercommand(const std::string& line);
                            usercommand(const usercommand &uc);
                            usercommand(const std::string& keyword,
                                        const std::string& option,
                                        const double value);
                            ~usercommand();
	
    int                     integerValue() const;
    const std::string&      keyword() const;
	const std::string&      option() const;
    double                  value() const;
    
	void                    promptUser();
	void                    print(std::ostream& of=std::cout) const;	    
    
    
private:
    std::string             _keyword;
    std::string             _option;
    double                  _numarg;
	bool                    _hasStringArgument;
    bool                    _hasNumberArgument;	
};




/* in the input file or elsewhere it is common to define a line with commands, separated by commas. All
 the commands in this line are held in a 'commandline'
 usage:
 
 commandLine cl;
 ...
 usercommand &uc = cl.userCommand(0);
 usercommand &uc2 = cl.userCommand(1);
 string uc.option();
 double uc.value();
 
 */


class commandLine : public std::vector<usercommand>
{

public:
                        commandLine();
                        commandLine(std::ifstream &file);
                        commandLine(const std::string& textline);
                        ~commandLine();
	
	friend std::ostream&  operator<<(std::ostream &os, const commandLine &cl);
	friend std::istream&  operator>>(std::istream &is, commandLine &cl);
    
    bool                hasKeyword(const std::string& word) const;
    std::string         optionFor(const std::string& word) const;
	void				print(std::ostream &of=std::cout) const;
	int					readFromStream(std::ifstream &file);
	size_t				readFromString(const std::string  &string);
    double              valueFor(const std::string& word) const;
	
	// scan a line from a stream ignoring comments and pasting continued lines
	static bool         scanCommandLine(std::istream &infile, std::string &line);

};

#endif
