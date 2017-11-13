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
 *  postprocessor.c
 *  feliks
 *
 *  Created by Ignacio Romero on Fri Dec 03 2004.
 *  Copyright (c) 2004 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include <memory>
#include <string>


#include "Io/Postprocessor/postprocessor.h"
#include "Io/logger.h"
#include "Io/usercommand.h"
#include "Io/message.h"
#include "Model/model.h"
#include "Io/Postprocessor/paraview.h"
#include "Io/Postprocessor/resultdata.h"
#include "smoother.h"



postprocessor :: postprocessor() :
	name(""),
	compress(false),
    evprint(false),
	activated(false),
	logfreq(0.0),
	lastlog(0.0),
	theSmoother(0)
{
}



postprocessor :: postprocessor(const commandLine &cl) :
	name(""),
	compress(false),
    evprint(false),
	activated(false),
	logfreq(0.0),
	lastlog(0.0),
	results(),
	theSmoother(0)
{
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if      (uc.keyword() == "compress"  && uc.option() == "no") compress = false;
        else if (uc.keyword() == "evalspots"  && uc.option() == "yes")evprint = true;
		else if (uc.keyword() == "frequency")		logfreq = uc.value();
		else if (uc.keyword() == "result")			results.push_back( uc.option()  );
		else if (uc.keyword() == "smoother"  && uc.option() == "spr") theSmoother = new SPRsmoother();
		else if (uc.keyword() == "smoother"  && uc.option() == "lump") theSmoother = new lumpsmoother();
	}		
	
	// default method for stress smoothing
	if (theSmoother == 0) theSmoother = new SPRsmoother();
}



bool postprocessor ::	dumpSolution(const model &m, const double time)
{
	logger::POSTtime.start();
	
	// if the frequency is 0.0, log always. Otherwise only if it's time
	extern double global_macheps;
	if (logfreq > 0.0 && time != 0.0 && nextLogTime() - time > global_macheps) return true;

	// call the specific postprocessor to do its thing...
	bool ret = dumpSolution2(m, time);
	
	// update the counter that keeps track of log outputs
	if (logfreq > 0.0) 
		while (time >= lastlog + logfreq)  lastlog += logfreq;
	logger::POSTtime.stop();

	return ret;
}



double postprocessor ::	nextLogTime() const
{
	return (lastlog + logfreq);
}



void postprocessor :: print(ostream &of)
{
	of << "\n\n";

	if (name == "")
	{
		of << " No postprocessor for output results " << "\n\n";
		return;
	}
	
	of	<< "\n\n";
	printCentered(of, " P o s t p r o c e s s o r");
	of 	<< endl << endl
		<< "\n Type          : " << name;
	if (logfreq > 0.0) 	of << "\n Log frequency : " << logfreq;
	else                of << "\n Log frequency : every time step";
	if (compress)       of << "\n Compress files: yes";
	else				of << "\n Compress files: no";
	of << "\n Smoothing     : " << theSmoother->getName(); 
	of << "\n Output vars.  :";
	if ( results.empty() ) of << " ----- " << endl;
	else
	{
		for (int a=0; a<results.size(); a++) of << "\t\t " << results[a] << "\n";
	}
	of << "\n";
}



void postprocessor :: scan(const commandLine &cl, std::auto_ptr<postprocessor> *p)
{
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if      (uc.keyword() == "type" && uc.option() == "paraview")	*p = std::auto_ptr<paraview>(new paraview(cl));
	}	
}
