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
 * logger.cpp
 *
 * i. romero, september 2001.
 * January 2004, redo the whole logging structure
 *
 * logs activity in feliks analysis
 *
 */

#include "Io/logger.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include "boost/foreach.hpp"

#include "Main/feliks.h"
#include "Io/usercommand.h"
#include "Model/Node/node.h"
#include "Model/Parts/body.h"
#include "Model/model.h"

#include <algorithm>
#include <functional>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <cctype>
#include <ctime>
#include <typeinfo>

#include "Io/message.h"
#include "Postprocessor/postprocessor.h"
#include "Analysis/analysis.h"
#include "Analysis/assembler.h"
#include "Analysis/energies.h"
#include "Model/Sets/elset.h"
#include "Model/Sets/nodeset.h"
#include "General/idata.h"
#include "General/timer.h"
#include "Math/vector.h"

#ifdef WITHMPI
#include "mpi.h"
#endif


using namespace std;


// we must define the static part of the class so that it exists
// even if no logger has been declared
vector<logger*>	logger :: allLoggers(0);
ofstream		logger :: sumf;
ofstream        logger :: mainlog;
ofstream        logger :: warnings;
timer			logger :: CPUtime;
timer           logger :: PREtime;
timer           logger :: ANAtime;
timer			logger :: PARSEtime;
timer			logger :: POSTtime;
time_t			logger :: walltime = 0;
string			logger :: summaryFilename = "";
string			logger :: mainlogFilename = "";
string          logger :: warningFilename = "";


// a couple of variables to keep track of the history of the solution. We want to make
// them invisible to the user
static int stepCount=0;
static int iterCount=0;


// length of a text line
const int linelength = 70;

static char* getHost();
static char* getOS();




/* the new way of defining logging variables. The syntax can be:
logger, variable = ..., [varnumber = ...], settype = ..., (set)label = [all]..., 
	name = ..., [frequency = ...], 
	[filename = ...], [xaxis = ...], [yaxis = ...], [precision = ...],
	[xlabel = ...] [, ylabel = ...]
*/
logger :: logger(const commandLine &clist) :
	active(true),
	vartype(LOGVAR_NONE) ,
	settype(LOGSET_NONE) ,
	setName(""),
	frequency(0.0),
	lasttime(0.0),
	logfilename(""),
	varnum(0),
	varx(0),
	vary(0),
	xlabel("Time"),
	ylabel(),
	precision(4),
    droptol(1e-8)
{
    usercommand uc;
    char         remove[30]="";
	
    /* first scan to see what type of variable we want to log
        It does not need to be the first data, but it must be there
        Start at k=1 because the first command is "output" */
    int k = 0;
    while (++k < clist.size() )
    {
        uc = clist[k];
        if (uc.keyword() == "variable")
        {
			if		(uc.option() ==  "bodies")   vartype = LOGVAR_BODIES;
            else if (uc.option() ==  "stress")   vartype = LOGVAR_GRADIENT;
            else if (uc.option() ==  "gradient") vartype = LOGVAR_GRADIENT;
            else if (uc.option() ==  "rate")     vartype = LOGVAR_RATE;
            else if (uc.option() ==  "velo")     vartype = LOGVAR_RATE;
            else if (uc.option() ==  "acce")     vartype = LOGVAR_ACCEL;
			else if (uc.option() ==  "connectivity") vartype = LOGVAR_CONNECTIVITY;
            break;
        }
    }
    
    
	// We proceed to scan the rest of the options
    k = 0;
    while (++k < clist.size())
    {
        uc = clist[k];
		
        if      (uc.keyword() == "varnumber")   varnum = (int) uc.value();
		else if (uc.keyword() == "filename")    logfilename = uc.option();
		else if (uc.keyword() == "frequency")   frequency = std::max<double>(0.0, uc.value());
		else if (uc.keyword() == "settype")
		{
			if		(uc.option() == "element")  settype = LOGSET_ELEMENT;
			else if (uc.option() == "elset")    settype = LOGSET_ELSET;
			else if (uc.option() == "nodeset")  settype = LOGSET_NODESET;
			else if (uc.option() == "node")		settype = LOGSET_NODE;
		}
		
		else if (uc.keyword() == "precision") precision = (int) uc.value();
		
		/*  define the (set)label of the set over which the log acts. It can be a number (referring
			to the corresponding node, nodeset, element, elset) or 'all' (if it refers to node or element */
		else if (uc.keyword() == "setlabel" || uc.keyword() == "label") setlabel = (int) uc.value();            
		else if (uc.keyword() == "name") 		setName = uc.option();
		else if (uc.keyword() == "xaxis") 	varx = (int) uc.value();
		else if (uc.keyword() == "yaxis") 	vary = (int) uc.value();		
		else if (uc.keyword() == "xlabel") 	xlabel=  uc.option();
		else if (uc.keyword() == "ylabel") 	ylabel=  uc.option();		
	}
	
	
	// some final checks. This should go in the specific constructor. Right now,
	// these are the remainders of the old C structs
	if (logfilename.empty())
	{
		ostringstream strm;
		if      (vartype == LOGVAR_BODIES)   strm << "feliks_bodies.vlog";
		else if (vartype == LOGVAR_STRESS)   strm << "feliks_str"  << allLoggers.size() << ".vlog";
		else if (vartype == LOGVAR_RATE)     strm << "feliks_rate" << allLoggers.size() << ".vlog";
		else if (vartype == LOGVAR_ACCEL)    strm << "feliks_acc"  << allLoggers.size() << ".vlog";
		
		logfilename = strm.str();
	}
	
	string removename = logfilename;
	modifyNameForMPI(removename);
	
	sprintf(remove, "rm -f %s", removename.c_str() );
	//cout << "\n name to erase " << removename << flush;
	system(remove); // erase old versions
}




// this function checks if the set on which the logger reports exists. If it does not, the logger
// is deactivated
void logger :: check(const model& theModel)
{
	if  (settype == LOGSET_NODE)
    {
		if ( &(theModel.getNode(setlabel)) == 0) active = false;
	}
    
	else if (settype == LOGSET_NODESET)
    {
		if ( &(theModel.getNodeset(setName)) == 0 ) active = false;
    }
}




void logger :: checkLoggers(const model& theModel)
{
	BOOST_FOREACH( logger* l, allLoggers) 	l->check(theModel);
}




void logger :: close()
{
	FILE *tmp=NULL;
	char tmpname[]="feliks.gnuplot";
	
	// if there is a graphic to be plotted
	if (vary > 0)
	{
		// open a temporary file where the gnuplot commands are written
		tmp = fopen(tmpname, "wr");
		
		// write gnuplot file
		fprintf(tmp, "set terminal pdf");
		fprintf(tmp, "\nset output \"%s.pdf\"", logfilename.c_str());
		fprintf(tmp, "\nset grid");
		fprintf(tmp, "\nset xlabel '%s';", xlabel.c_str());
		fprintf(tmp, "\nset ylabel '%s';", ylabel.c_str());
		fprintf(tmp, "\nplot \"%s\" using %d:%d with linespoints title '' ", logfilename.c_str(), varx+1, vary+1);
		fclose(tmp);
		
		// execute the gnuplot command
		system ("gnuplot feliks.gnuplot");
	}
}






void logger :: closeLogFiles()
{
    // current time
    time_t lt = time(NULL);
	CPUtime.stop();
    
    // close general log file
    if (mainlog.good())
    {
		mainlog.precision(3);
		mainlog << "\n\n Total elapsed time  = " << fixed << difftime(lt, walltime) << " seconds (wall time)."
				<< "\n Total CPU time      = " << CPUtime.seconds()      << " seconds."
				<< "\n Preprocess CPU time = " << PREtime.seconds()      << " seconds."
				<< "\n Parsing time        = " << PARSEtime.seconds()    << " seconds."
				<< "\n Postprocessing time = " << POSTtime.seconds()     << " seconds."
				<< "\n Solution CPU time   = " << ANAtime.seconds()      << " seconds."
		<< "\n\n Analysis finished on " << ctime(&lt);
		mainlog.close();
    }
    
    // close summary file
    if (sumf.good())
    {
        sumf << "\n" << "\n";
		sumf.precision(3);
		sumf << endl << endl
			<< "Time information" << "\n"
			<< "----------------" << "\n"
			<< " Analysis finished on " << ctime(&lt) << "\n"
			<< " Total elapsed time  = " << fixed << difftime(lt, walltime) << " seconds (wall time)." << "\n"
			<< " Total CPU time      = " << fixed << CPUtime.seconds()      << " seconds." << "\n"
			<< " Preprocess CPU time = " << fixed << PREtime.seconds()      << " seconds." << "\n"
			<< " Solution CPU time   = " << fixed << ANAtime.seconds()      << " seconds." << "\n";
		sumf.close();
    }
    
	
	// close warnings file
    if (warnings.good()) warnings.close();

	
	// close summary file
	string conflogFilename("feliks.conv") ;
#ifdef WITHMPI
	{	
        stringstream ss ; ss << MPI::COMM_WORLD.Get_rank() ; 
        conflogFilename.append(ss.str()) ; 	
	}
#endif
 	ofstream conf(conflogFilename.c_str(), ofstream::app);
	if (conf.good())
    {
        conf << "\n---------------------------------------------------------------------------" << "\n";
		conf << "Total number of time steps : " << stepCount << "\n";
		conf << "Total number of iterations : " << iterCount << "\n";
		conf.close();
    }
	
	
	// close all variable logs
	for_each( allLoggers.begin(), allLoggers.end(), mem_fun(&logger::close) );
	
    // Deallocate all the loggers
    freeLoggersList();
}








/* this is the main function that performs the logging.
* Evaluates a logger, and does a different thing depending on the variable
* logged, the set, etc...
*/
void logger ::  evaluate(const double time, const model &m, const assembler& theAssembler, const energies &energies)
{
    FILE    *fp=NULL;
    extern   double global_macheps, global_dt;
	double   gap, tolerance;
	
    // checks for quick return
    if (vartype == LOGVAR_NONE || !active) return;
    
    // if the frequency is 0.0, log always. Otherwise only if it's time
	// the tolerance depends (as stated) on the number of sums and current time
	gap  = time - lasttime - frequency;
	
	// some integrators do not have global_dt, so we use a reference unit time
	double rt ( global_dt > 0 ? global_dt : 1.0);
	tolerance = time*global_macheps * time/rt;
    if (frequency > 0.0 && time != 0.0 && gap < -tolerance) return;
    
	
    // open the corresponding log text file
    fp = fopen(logfilename.c_str(), "at");
	ofstream of(logfilename.c_str(), ofstream::app);
	of.precision(precision);
    
    vector<int> countCells(50);	
    switch(vartype)
    {
        case(LOGVAR_DOF):
            if      (settype == LOGSET_NODE)
            {
                if (time > 0.0) of << endl;
				of << scientific << time;
                const node &nd = m.getNode(setlabel);
                nd.printDOFs(of);
            }
			
            else if (settype == LOGSET_NODESET)
            {
                if (time > 0.0) of << endl;
				of << scientific << time;
				const nodeset &nds = m.getNodeset(setName);
					
                nds.printDofs(of);
            }
			
            break;
                  
	case(LOGVAR_BODIES):
            if (time > 0.0) of << endl;
                cout << "\n LOGVAR_BODIES disabled" << endl;
            /* 
             of << scientific << time << " " << setw(4) << m.theBodies.size()-1;

        BOOST_FOREACH(body* bd, m.theBodies )
        {
        	if (typeid( *bd ) == typeid( biocell ) )
			{
				biocell* c=static_cast<biocell*>( bd );
				int g = c->getGeneration();
				countCells[(g-1)]++;
			}	
		}
            
		for (int j=0; j<50; j++) of << " " << setw(4) << countCells[j];
             */
		break;

        case(LOGVAR_FORCE):
            break;
                        
        case(LOGVAR_RATE):
            if (time > 0.0) of << endl;
			of << scientific << time;
            if      (settype == LOGSET_NODE)
            {
                const node &nd = m.getNode(setlabel);
                nd.printRates(3, of);
            }
				
            else if (settype == LOGSET_NODESET)
            {
				const nodeset &nds = m.getNodeset(setName);
				nds.printRates(1, of);
            }
				
			break;
            
        case(LOGVAR_ACCEL):
            if (time > 0.0) of << endl;
			of << scientific << time;
            if      (settype == LOGSET_NODE)
            {
                const node &nd = m.getNode(setlabel);
				nd.printRates(2, of);
            }
				
            else if (settype == LOGSET_NODESET)
            {
				const nodeset &nds = m.getNodeset(setName);
				nds.printRates(2, of);
            }
								
            break;            
            
        default:
            break;
    }
    
    // make note that a evaluation was made
    lasttime = time;
    
    fclose(fp);
	of.close();
}





// this belongs to the new way of logging
/* deallocates all the memory of the loggers */
void logger :: freeLoggersList()
{
	for (int k=0; k< allLoggers.size(); k++)
		delete allLoggers[k];
	
	allLoggers.clear();
}





void logger :: listLoggers(ostream &of)
{
    if (allLoggers.size() > 0)
    {
        of  << endl << endl << endl;
		printCentered(of, "V a r i a b l e    L o g s");
        
		for (int k=0; k<allLoggers.size(); k++)
			allLoggers[k]->print(of);
	}
}




void logger :: logVariables(const double time, const model &theModel, const energies &energies, const assembler &theAssembler)
{
	vector<logger*>::iterator iter=allLoggers.begin();
	
	while( iter != allLoggers.end())
	{
		(*iter)->evaluate(time, theModel, theAssembler, energies);
		++iter;
	}
}

	

// this belongs to the new way of logging
/* create the default loggers, the ones that do the main log file, the convergence file, the
 summary file
 */
void logger :: openLogFiles(void)
{	
	startTimers();
    
	// remove old files. All the files that start with feliks
	string command("rm -f feliks[._-]*");
	modifyNameForMPI(command);
    system(command.c_str());	
	
	// main log file
	string mainlogFilename("feliks.log") ;
	modifyNameForMPI(mainlogFilename);
	mainlog.open(mainlogFilename.c_str());
	
	if (!mainlog) 
		ErrorMessage("in OpenLogFile. Cannot open file.");
	else
	{
		logger::printHeader(mainlog);
		
		// the call to ctime allocates some memory space that can not be freed ~ 900bytes
		// but it is only allocated once, even if it is called more times */
//        char host[256], login[256];
//        if (gethostname( host,  256) == 0) mainlog << "\n Hostname : " << host;
//		if (getlogin_r ( login, 256) == 0) mainlog << "\n Username : " << login;
    }
    
	// summary file
	string infoFilename("feliks.info") ;
	modifyNameForMPI(infoFilename);
	sumf.open(infoFilename.c_str());
	if (!sumf) 
		ErrorMessage("in OpenLogFiles. Cannot open summary file.");
	else
	{
        char host[256];
		logger::printHeader(sumf);
		sumf << "F EL I K S   summary file";
		sumf << "\n----------------------";
        sumf << "\nUsername             : " << getlogin();
        if (gethostname( host,  256) == 0) sumf << "\n Hostname : " << host;
		sumf << "\n\nAnalysis started on " << ctime(&walltime) << "\n\n";
	}
    
	// warnings file
	string warningsFilename("feliks.warnings");
	modifyNameForMPI(warningsFilename);
	warnings.open(warningsFilename.c_str());
	if (!warnings) 
		ErrorMessage("in OpenLogFiles. Cannot open summary file.");
	else
	{
		logger::printHeader(warnings);
		warnings << "F EL I K S  warnings file" << endl;
		warnings << "----------------------" << endl
		<< "Username             : " << getlogin() << "\n"
		<< "Hostname             : " << getHost() << "\n"
		<< "OS                   : " << getOS() << "\n";
		warnings << "Analysis started on " << ctime(&walltime) << "\n\n";
	}
	
    // file with convergence results
	ofstream conf;
	string conflogFilename("feliks.conv") ;
	modifyNameForMPI(conflogFilename);
	conf.open(conflogFilename.c_str());
	if (!conf)
        ErrorMessage("in OpenLogFiles. Cannot open convergence file.");
    else
    {
        conf <<  "F EL I K S   summary file of convergence results";
        conf <<  "\nAnalysis started on " << ctime(&walltime);
        conf << "\n  Step     Time          D-time    #iter  #trial   Init_resid    End_resid";
        conf << "\n---------------------------------------------------------------------------" << "\n";
		conf.close();
    }	
}



void logger :: print(ostream &of)
{
    of  << "\n\n Variable logger information";

	of  << "\n   Var type      : ";
    if      (vartype == LOGVAR_DOF)      of  << "Degree of freedom";
    else if (vartype == LOGVAR_GRADIENT) of  << "Gradient";
    else if (vartype == LOGVAR_FORCE)    of  << "Force";
    else if (vartype == LOGVAR_STRESS)   of  << "Stress";
    else if (vartype == LOGVAR_REACTION) of  << "Reaction";
    else if (vartype == LOGVAR_RATE)     of  << "Rate";
    else if (vartype == LOGVAR_ACCEL)    of  << "Acceleration";
    else if (vartype == LOGVAR_ENERGY)   of  << "Energy";
	else if (vartype == LOGVAR_CONNECTIVITY) of<< "Connectivity";
	
	of  << "\n   In set        : ";
    if (settype == LOGSET_NODE)     of  << "node " << setlabel;
    else if (settype == LOGSET_NODESET)  of  << "nodeset " << setlabel;
    else if (settype == LOGSET_ELEMENT)  of  << "element " << setlabel;
    else if (settype == LOGSET_ELSET)    of  << "elset " << setlabel;
    else if (settype == LOGSET_MESH)     of  << "the whole model";
    of  << "\n   Filename      : " << logfilename;
    
    
    of  << "\n   Var label     : ";
    if      (varnum  == -1)              of  << "all";
    else                                     of  << varnum+1;
    
    of  << "\n   Log frequency : ";
    if (frequency > 0.0)
        of  << frequency;
    else
        of  << "every time step";
    of  << "\n   Last log time : " << lasttime;
	of  << "\n   Precision     : " << precision;
    of  << "\n   Drop tolerance: " << droptol;

	of  << "\n   Graphic output: ";
	
	if (vary > 0)  
	{
		of << "\n\t   Abscissa : " << vary;
		of << "\n\t   Ordinate : " << varx;
		of << "\n\t   Abscissa label : " << xlabel;
		of << "\n\t   Ordinate label : " << ylabel;
	}
	else
		of << "none.";
	
	if ( !active) 
		of  << "\n   Warning       : logger inactive because set is not in model.";
}





void logger :: printHeader(ostream &of)
{
	of 
	<< "=============================================================================\n"
	<< "|                                                                           |\n"
	<< "|                              F    E    IL   K    S                        |\n"
	<< "|                                                                           |\n"
	<< "|                                                                           |\n"
	<< "|                                                                           |\n"
	<< "|             A general purpose finite element analysis program             |\n"
	<< "|                                                                           |\n"
	<< "|              Copyright Ignacio Romero, ignacio.romero@upm.es              |\n"
	<< "|                      Universidad Politecnica de Madrid                    |\n"
	<< "| version  " << FELIKS_VERSION << "  (" << FELIKS_DATE << ")" << setw(46) << "|" << endl
	<< "|                                                                           |\n"
	<< "=============================================================================\n\n"
	<< endl;
}



void logger :: scan(const commandLine &cl)
{
	logger *l(0);
	usercommand uc;
	string       option;
	
    /* first scan to see what type of variable we want to log
	 It does not need to be the first data, but it must be there
	 Start at k=1 because the first command is "output" */
    int k = 0;
    while (++k < cl.size() )
    {
        uc = cl[k];
        if (uc.keyword() == "variable")
        {
			if      (uc.option() ==  "energy")      l = new energylogger(cl);
         	else if (uc.option() ==  "gradient")	l = new gradientlogger(cl);
			else if (uc.option() ==  "dof")         l = new doflogger(cl);
			else if (uc.option() ==  "mass")		l = new masslogger(cl);
			else if (uc.option() ==  "reaction" || uc.option() == "reactions")	l = new reactionlogger(cl);
			//else if (uc.option() ==  "connectivity") l = new connectivityLogger(cl);
            break;
        }
    }
	
	// we need to subclass all logger types. Right only two missing
	if		( l == 0) l = new logger(cl);
	
	allLoggers.push_back(l);
}



void logger :: startAnalysisTimer()
{
	PREtime.stop();
	ANAtime.start();
}



void logger :: startTimers()
{
	// start the time that account for the CPU time spent in FELIKS and obtain current time
	CPUtime.start();
	PREtime.start();
    walltime = time(NULL);
}	




/* this function fills up a line of data in the convergence file. It contains:
the load step number, the (pseudo)time, the # of iterations, the initial residual,
the final residual */
void logger :: toConvFile(int step, double time, double dt, int niter, int ntrials, double res0, double resf)
{
	string conflogFilename("feliks.conv") ;
	modifyNameForMPI(conflogFilename);
 	ofstream conf(conflogFilename.c_str(), ofstream::app);
	if (!conf) 
		ErrorMessage("in OpenLogFiles. Cannot open summary file.");
	else
	{		
		conf << setw(5) << step << "  "
			<< setw(10) << scientific << time << "  " << dt << "  " 
			<<  setw(3)  << niter << "    " << ntrials << "     "
			<<  setw(11) << res0  << "  " << resf << "\n";
		conf.close();
	}
	stepCount++;
	iterCount += niter;
}





/*
 * this function just redirects the output of a printf command to th
 * log file that must have been already opened. Uses the standard function vfprintf
 * that works like fprintf but for a va_list argument. See K&R pg 174
 */
void logger :: toLogFile(char *fmt, ...)
{
    va_list args;
    FILE *logsf=NULL;
    
    logsf = fopen(mainlogFilename.c_str(), "at");
    va_start(args,fmt);
    vfprintf(logsf, fmt , args);
    fclose(logsf);
    
    va_end(args);
}





/* this is a modified ToLogFile function that write a line of the form:
*            first         :  rest
*  and is useful to get a unified set of messages from the analysis file
*/
void logger :: toLogFile2(const char *first, const char *fmt, ...)
{
    va_list args;
    FILE *logsf=NULL;
    
    logsf = fopen(mainlogFilename.c_str(), "at");
    
    fprintf(logsf, "\n  %-25s : ", first);
    va_start(args,fmt);
    vfprintf(logsf, fmt , args);
    va_end(args);
    fclose(logsf);
}



/* 
* this function just redirects the output of a printf command to th
 * log file that must have been already opened. Uses the standard function vfprintf
 * that works like fprintf but for a va_list argument. See K&R pg 174
 */
void logger :: vToLogFile(char *fmt, va_list args)
{
    FILE *logsf=NULL;
    
    logsf = fopen(mainlogFilename.c_str(), "at");
    vfprintf(logsf, fmt , args);
    fclose(logsf);
}




/*
 * this function just redirects the output of a printf command to th
 * summary file that must have been already opened. Uses the standard function vfprintf
 * that works like fprintf but for a va_list argument. See K&R pg 174
 */
void logger :: toSummaryFile(char *fmt, ...)
{
}




/* 
* this function just redirects the output of a printf command to th
 * sum file that must have been already opened. Uses the standard function vfprintf
 * that works like fprintf but for a va_list argument. See K&R pg 174
 */
void logger :: vToSummaryFile(char *fmt, va_list args)
{
}





doflogger :: doflogger(const commandLine& cl) :
	logger(cl),
	P2nodeset(0),
	P2node(0),
    bodyname(),
    dimension(-1),
    manifoldLabel(-1)
{
	vartype = LOGVAR_DOF;
	
    
    if (settype == LOGSET_NONE)
    {
        int k = 0;
        while (++k < cl.size() )
        {
            usercommand uc = cl[k];
            if      (uc.keyword() == "body")         bodyname      = uc.option();
            else if (uc.keyword() == "dimension")    dimension     = uc.integerValue();
            else if (uc.keyword() == "label")        manifoldLabel = uc.integerValue();
        }
        
        if (!bodyname.empty() && dimension >= 0 && manifoldLabel >= 0)
            settype = LOGSET_MANIFOLD;
    }
    
    
	if (logfilename.empty() ) 
	{
		ostringstream strm;		
		strm << "feliks_dof" << allLoggers.size() << ".vlog";
		logfilename = strm.str();
	}
	
	modifyNameForMPI(logfilename);	
}




void doflogger ::  evaluate(const double time, const model &m, 
                            const assembler& theAssembler, const energies &energies)
{
	if (!active) return;
	
    extern   double global_macheps, global_dt;
	double   gap, tolerance;
        
    // if the frequency is 0.0, log always. Otherwise only if it's time
	// the tolerance depends (as stated) on the number of sums and current time
	gap = time - lasttime - frequency;
	
	// some integrators do not have global_dt, so we use a reference unit time
	double rt ( global_dt > 0 ? global_dt : 1.0);
	tolerance = time*global_macheps * time/rt;
    if (frequency > 0.0 && time != 0.0 && gap < -tolerance) return;
    
    // open the corresponding log text file
	ofstream of(logfilename.c_str(), ofstream::app);
	of.precision(precision);
	if (time > 0.0) of << "\n";
    
	if   (settype == LOGSET_NODE)
	{
		if (P2node == 0) P2node = &( m.getNode(setlabel) );
		of << scientific << time;
		P2node->printDOFs(of);
	}
	
	else if (settype == LOGSET_NODESET)
	{
		if (P2nodeset == 0) P2nodeset = &( m.getNodeset(setName) );
		of << scientific << time;
		P2nodeset->printDofs(of);
	}
    
				
    // make note that a evaluation was made
    lasttime = time;
    
	of.close();
}



void doflogger :: print(ostream &of)
{
    of  << "\n\n Variable logger information";
	
	of  << "\n   Var type      : Degrees of freedom.";
    of  << "\n   Filename      : " << logfilename;
    of  << "\n   Log frequency : ";
    if (frequency > 0.0)
        of  << frequency;
    else
        of  << "every time step";
	
	of  << "\n   In set        : ";
    if      (settype == LOGSET_NODE)     of  << "node " << setlabel;
    else if (settype == LOGSET_NODESET)  of  << "nodeset " << setName;
    else if (settype == LOGSET_MANIFOLD)
    {
        of << "body " << bodyname << ", manifold dimension: " << dimension << ", label: " << manifoldLabel;
    }

	of  << "\n   Precision     : " << precision;
	of  << "\n   Graphic output: ";
	if (vary > 0)  
	{
		of << "\n\t   Abscissa : " << vary;
		of << "\n\t   Ordinate : " << varx;
		of << "\n\t   Abscissa label : " << xlabel;
		of << "\n\t   Ordinate label : " << ylabel;
	}
	else
		of << "none.";
	
	of << flush;
}




energylogger :: energylogger(const commandLine& cl) :
	logger(cl)
{
	vartype = LOGVAR_ENERGY;
	if (logfilename.empty() ) logfilename = "feliks_energy.vlog";
	
	modifyNameForMPI(logfilename);
}



void energylogger ::  evaluate(const double time, const model &m, const assembler& theAssembler, const energies &theEnergies)
{
    extern   double global_macheps, global_dt;
	double   gap, tolerance;
    
    // if the frequency is 0.0, log always. Otherwise only if it's time
	// the tolerance depends (as stated) on the number of sums and current time
	gap = time - lasttime - frequency;
	
	// some integrators do not have global_dt, so we use a reference unit time
	double rt ( global_dt > 0 ? global_dt : 1.0);
	tolerance = time*global_macheps * time/rt;
    if (frequency > 0.0 && time != 0.0 && gap < -tolerance) return;
    
    // open the corresponding log text file
	ofstream of(logfilename.c_str(), ofstream::app);
	of.precision(precision);
    
	// the first line is treated differently the first time, to allow header in gnuplot format
	if (time == 0.0) 
	{
		of << "#  time  ";
		theEnergies.printNames(of);
	}
	of << endl << scientific << setw(10) << time;
    
    energies cleanEnergies = theEnergies;
    cleanEnergies.chop(droptol);
    
	cleanEnergies.print(of);
	of.close();	
    
    // make note that a evaluation was made
    lasttime = time;    
}



void energylogger :: print(ostream &of)
{
    of  << "\n\n Variable logger information";
	
	of  << "\n   Var type      : Energy, entropy, and momenta.";
	of  << "\n   In set        : the whole domain.";
    of  << "\n   Filename      : " << logfilename;
    of  << "\n   Log frequency : ";
    if (frequency > 0.0)
        of  << frequency;
    else
        of  << "every time step";
    of  << "\n   Last log time : " << lasttime;
	of  << "\n   Precision     : " << precision;	
	of  << "\n   Graphic output: ";
	if (vary > 0)  
	{
		of << "\n\t   Abscissa : " << vary;
		of << "\n\t   Ordinate : " << varx;
		of << "\n\t   Abscissa label : " << xlabel;
		of << "\n\t   Ordinate label : " << ylabel;
	}
	else
		of << "none.";
}





gradientlogger :: gradientlogger(const commandLine& cl) :
	logger(cl)
{
	vartype = LOGVAR_GRADIENT;

	// some final checks
	if (logfilename.empty())
	{
		ostringstream strm;
		strm << "feliks_grad" << allLoggers.size() << ".vlog";
		logfilename = strm.str();
	}
	
	modifyNameForMPI(logfilename);
}




/* this is the main function that performs the logging.
 * Evaluates a logger, and does a different thing depending on the variable
 * logged, the set, etc...
 */
void gradientlogger ::  evaluate(const double time, const model &m, const assembler& theAssembler, const energies &energies)
{
    extern   double global_macheps, global_dt;
	double   gap, tolerance;
    
    // if the frequency is 0.0, log always. Otherwise only if it's time
	// the tolerance depends (as stated) on the number of sums and current time
	gap  = time - lasttime - frequency;
	
	// some integrators do not have global_dt, so we use a reference unit time
	double rt ( global_dt > 0 ? global_dt : 1.0);
	tolerance = time*global_macheps * time/rt;
    if (frequency > 0.0 && time != 0.0 && gap < -tolerance) return;
    
    // open the corresponding log text file
	ofstream of(logfilename.c_str(), ofstream::app);
	of.precision(precision);

    
	if (time > 0.0) of << "\n";
	of << scientific << setw(10) << time;

	if (settype == LOGSET_ELSET)
	{
		elset &es(m.getElset(setName));
		es.printGradients(of);
	}
	
	else if (settype == LOGSET_ELEMENT)
	{
		element &el(m.getElement(setlabel));
		el.printGradients(of);
	}
	        
    
    // make note that a evaluation was made
    lasttime = time;
    
	of.close();	
}


void gradientlogger :: print(ostream &of)
{
	logger::print(of);
}





masslogger :: masslogger(const commandLine& cl) :
logger(cl)
{
	vartype = LOGVAR_MASS;
	if (logfilename.empty() ) logfilename = "feliks_mass.vlog";
	
	modifyNameForMPI(logfilename);
}



void masslogger ::  evaluate(const double time, const model &m, const assembler& theAssembler, const energies &theEnergies)
{
    extern   double global_macheps, global_dt;
	double   gap, tolerance;
    
    // if the frequency is 0.0, log always. Otherwise only if it's time
	// the tolerance depends (as stated) on the number of sums and current time
	gap = time - lasttime - frequency;
	
	// some integrators do not have global_dt, so we use a reference unit time
	double rt ( global_dt > 0 ? global_dt : 1.0);
	tolerance = time*global_macheps * time/rt;
    if (frequency > 0.0 && time != 0.0 && gap < -tolerance) return;
    
    // open the corresponding log text file
	ofstream of(logfilename.c_str(), ofstream::app);
	of.precision(precision);
    
	// the first line is treated differently the first time, to allow header in gnuplot format
	if (time == 0.0) of << "#     time       mass    volume";
	of << "\n" << scientific << setw(12) << time << setw(12) << theEnergies.mass << setw(12) << theEnergies.volume;
	of.close();	
    
    // make note that a evaluation was made
    lasttime = time;    
}



void masslogger :: print(ostream &of)
{
    of  << "\n\n Variable logger information";
	
	of  << "\n   Var type      : Mass";
	of  << "\n   In set        : the whole domain.";
    of  << "\n   Filename      : " << logfilename;
    of  << "\n   Log frequency : ";
    if (frequency > 0.0)
        of  << frequency;
    else
        of  << "every time step";
    of  << "\n   Last log time : " << lasttime;
	of  << "\n   Precision     : " << precision;	
	of  << "\n   Graphic output: ";
	if (vary > 0)  
	{
		of << "\n\t   Abscissa : " << vary;
		of << "\n\t   Ordinate : " << varx;
		of << "\n\t   Abscissa label : " << xlabel;
		of << "\n\t   Ordinate label : " << ylabel;
	}
	else
		of << "none.";
}






reactionlogger :: reactionlogger(const commandLine& cl) :
	logger(cl),
	P2nodeset(0),
	P2node(0)
{
	vartype = LOGVAR_REACTION;
	
    if (settype == LOGSET_NONE)
    {
        int k = 0;
        while (++k < cl.size() )
        {
            usercommand uc = cl[k];
            if      (uc.keyword() == "body")         bodyname      = uc.option();
            else if (uc.keyword() == "dimension")    dimension     = uc.integerValue();
            else if (uc.keyword() == "label")        manifoldLabel = uc.integerValue();
        }
        
        if (!bodyname.empty() && dimension >= 0 && manifoldLabel >= 0)
            settype = LOGSET_MANIFOLD;
    }
    
	// some final checks
	if (logfilename.empty())
	{
		ostringstream strm;
		strm << "feliks_reac" << allLoggers.size() << ".vlog";
		logfilename = strm.str();
	}
	
	modifyNameForMPI(logfilename);
}



void reactionlogger ::  evaluate(const double time, const model &m, const assembler& theAssembler, const energies &energies)
{
	extern   double global_macheps, global_dt;
	double   gap, tolerance;
    
    // if the frequency is 0.0, log always. Otherwise only if it's time
	// the tolerance depends (as stated) on the number of sums and current time
	gap  = time - lasttime - frequency;
	
	// some integrators do not have global_dt, so we use a reference unit time
	double rt ( global_dt > 0 ? global_dt : 1.0);
	tolerance = time*global_macheps * time/rt;
    if (frequency > 0.0 && time != 0.0 && gap < -tolerance) return;
    
    // open the corresponding log text file
	ofstream of(logfilename.c_str(), ofstream::app);
	of.precision(precision);
	
    
	if (time > 0.0) of << "\n";
	of << scientific << setw(10) << time;
	
	model *mm;
	mm = const_cast<model*>(&m);
	
	theAssembler.assembleReactions(*mm, time);
	if(settype == LOGSET_NODE)
	{
		if (P2node == 0) P2node = &( m.getNode(setlabel) );
		of << scientific << std::showpos << right << setw(of.precision()+8)
            << "    "
			<< const_cast<node*>(P2node)->force()(0) << "    " 
			<< const_cast<node*>(P2node)->force()(1) << "    "
			<< const_cast<node*>(P2node)->force()(2);
	}
	
	else if (settype == LOGSET_NODESET)
	{
		if (P2nodeset == 0) P2nodeset = &( m.getNodeset(setName) );
        
        if (P2nodeset->size() > 0)
        {            
            int dofsPerNode = (*(P2nodeset->begin()))->getNDofs();
            double *ftot = new double[dofsPerNode];
            for (int i=0; i<dofsPerNode; i++) ftot[i] = 0.0;
            
            extern double global_macheps;	
            
            std::set<node*>::const_iterator iter = P2nodeset->begin();	
            while (iter != P2nodeset->end() )  
            {
                for (int i=0; i<dofsPerNode; i++)
                    ftot[i] += (*iter)->force()(i);
                ++iter;
            }
            
            for (int i=0; i<dofsPerNode; i++)
            {
                double f = fabs( ftot[i] ) > 1000.0*global_macheps ? ftot[i] : 0.0;
                of << " " << scientific << std::showpos << right << setw(of.precision()+8) << f;
            }
            
            delete [] ftot;
        }
	}
    
    
    // make note that a evaluation was made
    lasttime = time;
	of.close();	
}


/*
connectivityLogger :: connectivityLogger (const commandLine& cl) :
	logger(cl),
	hasbeenCalled(false)
{
	vartype = LOGVAR_CONNECTIVITY;
    settype = LOGSET_MESH;
	
	if (logfilename.empty()) logfilename = "feliks.connect";
}



void connectivityLogger ::  evaluate(const double time, 
                                     const model &m, 
                                     const assembler& theAssembler,
                                     const energies &energies)
{
	if (hasbeenCalled) return;
	hasbeenCalled = true;

	ofstream of;
	of.open(logfilename.c_str());

	model &mm = const_cast<model&>(m);

	BOOST_FOREACH( const body *b, mm.theBodies ) printBodyNodeConnectivity(*b, of);
	
	of.close();
}



void connectivityLogger :: print(ostream& of)
{
    of  << "\n\n Variable logger information";
	
	of  << "\n   Var type      : Connectivity";	
	of  << "\n   In set        : Whole model";
    of  << "\n   Filename      : " << logfilename;
}
*/



static char* getHost()
{
	static char nohost[] = "UNKNOWN HOST";
	char *host = getenv("HOST");
	
	if (host == NULL) host = getenv("HOSTNAME");
	if (host == NULL) host = nohost;
	return host;
}



static char* getOS()
{
	static char noOS[] = "UNKNOWN OS";
	char *os = getenv("OSTYPE");
	
	if (os == NULL) os = getenv("OS");
	if (os == NULL) os = noOS;
	return os;
}


void modifyNameForMPI(string& logfilename)
{
#ifdef WITHMPI
	stringstream ss ; ss << MPI::COMM_WORLD.Get_rank();
	logfilename.append(ss.str()) ; 	
#endif
}


std::ostream& printCentered(std::ostream& of, const std::string &message)
{
	int len = message.length();
	int white(linelength-len);
	
	string gap(white/2,' ');
	string centered(gap);
	centered += message;
	centered += gap;
	of << centered;
	return of;
}



