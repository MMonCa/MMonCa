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
 * helper.c
 *
 * help commands for feliks
 *
 * i. romero, january 2002
 */

#include "Io/helper.h"

#include <iostream>
#include <fstream>
#include <string>

#include <cstdlib>
#include <cstring>
#include "Io/message.h"
#include "Main/feliks.h"

using namespace std;


// a unique object 
helper theHelper;


void helper :: describeFELIKS() const
{
  printf("\n\n  F EL I K S, version %4.1f, %s", FELIKS_VERSION, FELIKS_DATE); /* defined in feliks.h */
  printf("\n  A general purpose finite element code for research and development.");
  printf("\n\n  By: Ignacio Romero");
  printf("\n      E.T.S. Ingenieros Industriales");
  printf("\n      Universidad Politecnica de Madrid");
  printf("\n      ignacio.romero@upm.es\n\n");
}





/* a very basic help message, showing the major help options */
void helper ::  helpBasics(void) const
{
  printf("\n\n Basic help for FELIKS:");
  printf("\n  To get a list of all commands, type 'help commands'.");
  printf("\n  To get some help about command xxx, type 'help xxx'.");
  printf("\n  To obtain a brief description of the code, type 'help feliks'.");
  printf("\n  To exit the program, type 'quit'.\n"); 
}




void helper :: commandList(void) const
{
  printf("\n              Command List");
  printf("\n              ------------");
  printf("\n  assemble  dump  help  info  integration option plot quit show  solve"); 
}





/* This function gives more specific help about a command
   that is given in the "command" argument
*/
void helper :: describe(const std::string& command) const
{
	if (command == "")						helpBasics();
	else if ( command == "commands" )		commandList();
	else if ( command == "assemble")		printHelpFile("assemble.txt");
	else if ( command == "dump")			printHelpFile("dump.txt");
	else if ( command == "info")			printHelpFile("info.txt");
	else if ( command == "feliks")			describeFELIKS();
	else if ( command == "integration")		printHelpFile("integration.txt");
	else if ( command == "options")			printHelpFile("option.txt");
	else if ( command == "quit")			printHelpFile("quit.txt");
	else if ( command == "show")			printHelpFile("show.txt");
	else if ( command == "solve")			printHelpFile("solve.txt");
	else cout << "There is no help for command '%s'." << command << endl;
}



/* reads a file that must be in the directory Flow/HelpFiles
   and prints the information in the screen.
*/
void helper :: printHelpFile(const std::string& filename) const
{
	std::string line, path;
	
	if (filename.empty()) 
		ErrorMessage("in PrintHelpFile. Empty filename");
	
	else
    {
		path  = FELIKSPATH;
		path += "/Io/HelpFiles/" + filename;
		
		ifstream f(path.c_str());
		
		if ( !f)
			printf("Help file '%s' does not exist", path.c_str());
		else
		{
			printf("\n--------------------------------------------------------------------------");
			while ( getline(f, line) )   cout << endl << line;
			printf("\n--------------------------------------------------------------------------");
			printf("\n");
			f.close();
		}
	}
}








