/*
 * CascadeEvent.cpp
 *
 *  Created on: May 6, 2015
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

#include "CascadeEvent.h"
#include "domains/MCClient.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include "kernel/Mesh.h"
#include "kernel/SubDomain.h"
#include "kernel/Analyzer.h"
#include <sstream>

using std::string;
using std::vector;
using std::stringstream;
using std::map;

namespace Cascade {

CascadeEvent::CascadeEvent(Kernel::Domain *pD, const std::vector<std::string> &f,
		bool bReact, bool periodic, bool voluminic, bool bDisplace, bool bCorrectX,
		const std::string &fileName, float rate)  :
				Event(pD), _format(f), _bReact(bReact), _periodic(periodic), _voluminic(voluminic),
				_bDisplace(bDisplace), _bCorrectX(bCorrectX), _rate(rate)
{
	_pist.open(fileName.c_str());
	if(_pist.fail())
		ERRORMSG("cascade: file " << fileName << " not found");

	//Algorithm to store the random positions within the file. It assumes that the first line of every cascade starts with a # New cascade
	//It works if there are some headers. The only condition is that every line which has only one number, n,
	//(no matter if it is surrounded by spaces or tabs) and the rest lines of that cascade are OK.

	bool bFormat = false;
	unsigned position = 0;
	string theLine;
	getline(_pist, theLine);
	if(theLine == "# New cascade")
	{
		_startPos.push_back(position);
		do
		{
			position = _pist.tellg();
			getline(_pist, theLine);
			if(theLine == "# New cascade")
				_startPos.push_back(position);
		} while(!_pist.eof());
	}
	else
		do
		{
			_startPos.push_back(position);
			stringstream ss;
			unsigned numberOf;
			ss << theLine;
			ss >> numberOf;
			while(numberOf--)
			{
				getline(_pist,theLine);
				if(_pist.eof())
					ERRORMSG("Error in cascade file. Format 'Number of particles' seems to be wrong. Unexpected end of file");
			}
			bFormat = true;
			position = _pist.tellg();
			getline(_pist, theLine);
		} while(!_pist.eof());

	LOWMSG(_startPos.size() << " cascades detected. Format detected is " << (bFormat? "Number of particles":"'# New cascade' tag."));
	if(_startPos.empty())
		ERRORMSG("File " << fileName << " contains no cascades. Are the new cascades being started with '# New cascade' or the number of particles?");
	_cPos = _startPos;
}

CascadeEvent::~CascadeEvent()
{
	_pist.close();
}

void CascadeEvent::perform(Kernel::SubDomain *p, unsigned)
{
	p->_evLog.performed(0, Kernel::Event::CASCADE, 0, 0, 0, 0);
	string firstLine;
	stringstream firstBuffer;

	Kernel::Mesh *pMesh = p->_pDomain->_pMesh;
	Kernel::Coordinates m,M;
	pMesh->getDomain(m,M);
	float x=0, y=0, z=0;
	unsigned randomPos = p->_rng.rand()*_startPos.size();
	assert(randomPos < _startPos.size());
	if(_pist.eof())
		_pist.clear();

	_pist.seekg(_startPos[randomPos],_pist.beg);
	getline(_pist,firstLine);
	firstBuffer << firstLine;
	unsigned numberOfParticles = 0;
	if(firstLine != "# New cascade")	//format number of particles
	{
		stringstream ss;
		ss << firstLine;
		ss >> numberOfParticles;
	}

	if(_voluminic)
		x = m._x + p->_rng.rand()*(M._x - m._x);
	if(_bDisplace)
	{
		y = m._y + p->_rng.rand()*(M._y - m._y);
		z = m._z + p->_rng.rand()*(M._z - m._z);
	}
	if(_bCorrectX)
	{
		Kernel::Coordinates cell, Cell;
		double minX = M._x;
		for(Kernel::Mesh::iterator mi =pMesh->begin(); mi != pMesh->end(); ++mi)
		{
			mi->getCorners(cell, Cell);
			if(cell._y < y && Cell._y > y && cell._z < z && Cell._z > z)
				if(!Domains::global()->PM()->isGasLike(mi->getMaterial()) && cell._x < minX)
					minX = cell._x;
		}
		x = minX;
	}

	//for each particle in the cascade
	unsigned howMany = 0;
	while(true)
	{
		string line;
		stringstream buffer;
		vector<string> fieldsInFile;
		getline(_pist, line);
		if(_pist.eof() || line == "# New cascade" || (numberOfParticles && numberOfParticles == howMany++))
			break;
		buffer << line;
		do
		{
			string token;
			buffer >> token;
			fieldsInFile.push_back(token);
		} while(!buffer.eof());
		if(buffer.fail())
			fieldsInFile.pop_back();
		create(pMesh, x, y, z, _format, fieldsInFile, m, M, _bReact, _periodic, _voluminic);
	}

	//removing the cascade position, so it is not repeated.
	if(_startPos.size() > 1)
	{
		_startPos[randomPos] = _startPos.back();
		_startPos.pop_back();
	}
	else
		_startPos = _cPos;
}

/* Format are 4 strings, each of them containing the format of:
 * 	1-> Point Defect
 * 	2, 3, 4 -> x,y,z
 * fileLine is an array of strings in one file line
 *
 * The function has to transform the first format in the particle type, generate variables for each column
 * in filename, and evaluate the 2nd, 3rd and 4th format
 */
void CascadeEvent::create(Kernel::Mesh *pMesh,
		float impactx, float impacty, float impactz,
		const vector<string> &format, const vector<string> &fileLine,
		const Kernel::Coordinates &m, const Kernel::Coordinates &M,
		bool bReact, bool periodic, bool voluminic)
{
	unsigned typeColumn = format[0][0] - 'A';
	//create variable map
	map<string, double> variableMap;
	for(unsigned column = 0; column < fileLine.size(); ++column)
	{
		string varName = "A";
		varName[0] += column;
		variableMap[varName] = atof(fileLine[column].c_str());
		HIGHMSG(varName << " = " << atof(fileLine[column].c_str()) << "; ");
	}

	Kernel::Analyzer<double> xAnalyzer(format[1]);
	double x = xAnalyzer(variableMap);
	Kernel::Analyzer<double> yAnalyzer(format[2]);
	Kernel::Analyzer<double> zAnalyzer(format[3]);
	double y = yAnalyzer(variableMap);
	double z = zAnalyzer(variableMap);

	x += impactx;
	y += impacty;
	z += impactz;

	if(periodic)
	{
		float sizeY = M._y - m._y;
		float sizeZ = M._z - m._z;
		while(y >= M._y) y -= sizeY;
		while(y <  m._y) y += sizeY;
		while(z >= M._z) z -= sizeZ;
		while(z <  m._z) z += sizeZ;
		if(voluminic)
		{
			float sizeX = M._x - m._x;
			while(x >= M._x) x -= sizeX;
			while(x <  m._x) x += sizeX;
		}
	}

	Kernel::Coordinates c(x, y, z);
	if(c.isInto(m, M))
	{
		unsigned index = pMesh->getIndexFromCoordinates(c);
		Kernel::M_TYPE mt = pMesh->getElement(index)->getMaterial();
		if(Domains::global()->PM()->getMaterial(mt)._bAmorphous )
			if(fileLine[typeColumn] == "I" || fileLine[typeColumn] == "V" )
				pMesh->getElement(index)->incrAmorphParts();

		const string name = fileLine[typeColumn];
		Kernel::Event::E_TYPE et = Domains::global()->PM()->getDefectType(name, mt);
		switch(et)
		{
		case Kernel::Event::MOBILEPARTICLE:
			Domains::global()->client()->createMP(name, "", c, bReact);
			break;
		default:
		  {
			vector<string> tokens;
			IO::ParameterManager::getTokens(name, ':', tokens);
			if(tokens.size() != 2 || Domains::global()->PM()->getDefectType(tokens[1], mt) != Kernel::Event::CLUSTER)
			{
				WARNINGMSG("Cannot recognize " << name << " as a cluster. ':' expected!");
				return;
			}
			Domains::global()->client()->createMC(tokens[1], tokens[0], c);
		  }
		}
	}
	else //be sure it is out!
		Domains::global()->client()->createMP(fileLine[typeColumn],"", Kernel::Coordinates(1e19,1e19,1e19), bReact);
}


} /* namespace Cascade */
