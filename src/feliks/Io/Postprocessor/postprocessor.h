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
/* postprocessor.h
 *
 *
 *
*/

#ifndef _postprocessor_h
#define _postprocessor_h


#include <string>
#include <iostream>
#include <memory>
#include "Math/matrix.h"
#include <vector>
#include "Io/Postprocessor/resultdata.h"


class commandLine;
class model;
class smoother;

class postprocessor{

protected:
	std::string name;
	bool    compress;
	bool    activated;
    bool    evprint;
	double  logfreq;
	double  lastlog;
	std::vector<std::string> results;
	smoother *theSmoother;
		
	
public:

	postprocessor();
	postprocessor(const commandLine &cl);
	virtual ~postprocessor(){};
	
	
	inline std::string getName()	{return name;}

	//void  log(double time, model &m)=0;
	virtual void	print(std::ostream &of=std::cout);
	
	virtual bool	dumpMesh(const model &m) = 0;
    virtual bool	dumpModes(const model &m, matrix evectors)=0;
	bool			dumpSolution(const model &m, const double time);
	virtual bool	dumpSolution2(const model &m, const double time)=0;

	double			nextLogTime() const;
	static  void	scan(const commandLine &cl, std::auto_ptr<postprocessor> *p);
};
	

#endif

