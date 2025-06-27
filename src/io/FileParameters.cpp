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

#include "FileParameters.h"
#include <sstream>
#include <algorithm>
#include <dirent.h>
#include <sys/stat.h>
#include <fstream>
#include "Diagnostic.h"
 
using namespace IO;
using std::string;
using std::map;
using std::endl;
using std::pair;

FileParameters::FileParameters(const std::string &szPath)
{
	openDir(szPath, "");
}

FileParameters::~FileParameters()
{
	for(map<string, bool>::iterator it=_used.begin(); it!= _used.end(); ++it)
		if(it->second == false)
			MEDMSG("Parameter " << it->first << " not used");
}

void FileParameters::dump(std::ostream &aOut, std::string const& aFilename)
{
	if(!aFilename.empty()) {
		aOut << "written into file: " << aFilename;
		std::ofstream file(aFilename);
		for(auto const& item : _map)
		{
			file << item.first << ':' << item.second.first << ':' << item.second.second << '\n';
		}
	}
	else {
		for(auto const& item : _map)
		{
			aOut << item.first << ':' << item.second.first << ':' << item.second.second << '\n';
		}
	}
}

void FileParameters::openDir(const std::string &szPath, const std::string &key)
{
   DIR *pDir = opendir(szPath.c_str());
   dirent *pDirent;
   if(pDir == 0)
	ERRORMSG("Cannot open " << szPath);   
   MEDMSG("Reading config in " << key);
   insert(key, "dir", "");
   while ((pDirent=readdir(pDir))) 
   {
      std::string entryText(pDirent->d_name);
      if(entryText[0] == '.')
      	continue;
      std::string newKey;
      if(key.empty())
      	newKey  = entryText;
      else
        newKey = key + '/' + entryText;
      struct stat buffer;
      stat((szPath+'/'+entryText).c_str(), &buffer);
      if(S_ISREG(buffer.st_mode))
    	  openFile(szPath+'/' + entryText, newKey);
      else if(S_ISDIR(buffer.st_mode))
    	  openDir(szPath+'/' + entryText, newKey);
      else
    	  LOWMSG("Cannot understand " << entryText);

   }
   closedir(pDir);
}

void FileParameters::openFile(const std::string &szFileName, const std::string &key)
{
   std::ifstream ifs(szFileName.c_str());
   if(ifs.fail())
    ERRORMSG("Can't open " << szFileName);
   insert(key, "file", "");
   int nLine = 0;
   while(!ifs.eof())
   {
	string name, szType, value;
	ifs >> szType >> name;
	if(ifs.eof())
		break;
	std::string newKey(key+'/'+name);
	std::getline(ifs, value);
	while(std::isspace(value[0]))
		value.erase(value.begin(), value.begin()+1);
        if(value.size() == 0)
                continue;
	if(value[0] == '{') //read until another '}'
	{
		string totalValue("");
		string buffer = value;
		string::iterator it=buffer.begin()+1;
		unsigned count = 1;
		while(count)
		{
			while(it == buffer.end())
			{
				std::getline(ifs, buffer);
				buffer += '\n';
				it=buffer.begin();
				if(ifs.eof())
					ERRORMSG("End of file " << szFileName << " while looking for '}'");
			}
			if(*it == '{')
				count++;
			if(*it == '}')
				count--;
			if(count)
				totalValue += *it++;
		}
		value = totalValue;
	}

	if(value.size())
        while(value[value.length()-1] == '\\')
	{
		value[value.length()-1] = ' ';
		string buffer;
		std::getline(ifs, buffer);
		if(ifs.eof())
			ERRORMSG("End of file " << szFileName << " while looking for the line after '\\'");
		value += buffer;
	}

	MEDMSG(nLine++ << " " << newKey << " = (" << szType << ") '" << value << "'");
	if(szType == "bool"      || szType == "int"                || szType == "float" ||
	   szType == "string"    || szType == "map<string,int>" ||
	   szType == "map<string,string>" || szType == "map<string,float>" ||
	   szType == "arrhenius" || szType == "map<string,bool>"   ||
	   szType == "map<string,arrhenius>" ||
	   szType == "array<string,string>" || szType == "array<string>" ||
	   szType == "proc" || szType == "coordinates")
	{
		insert(newKey, szType, value);
		_used[newKey] = false;
	}
	else if(szType == "comment" || szType == "//" || szType == "#")
		continue;//do nothing
	else
	     ERRORMSG(szType << " : syntax error in " << szFileName);
    }
}
