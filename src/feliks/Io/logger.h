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
 * logger.h
 *
 * i. romero, september 2001
 *
 * logs activity in feliks analysis
 *
 */

#ifndef _logger_h
#define _logger_h


#include <string>
#include <vector>
#include <ostream>

#include <cstdarg>
#include <cstdio>


#include "Analysis/analysis.h"
#include "Analysis/energies.h"
#include "General/timer.h"
#include "Model/Sets/nodeset.h"



class assembler;
class commandLine;
class doflogger;
class model;


/* type of set where the logger acts on */
enum logset{
    LOGSET_NONE,
    LOGSET_NODE,
    LOGSET_NODESET,
    LOGSET_ELEMENT,
    LOGSET_ELSET,
    LOGSET_MESH,
    LOGSET_MANIFOLD
};




/* quantity that is logged */
enum logvar{
    LOGVAR_NONE,
	LOGVAR_BODIES,	// number of bodies in the model
    LOGVAR_DOF,
    LOGVAR_GRADIENT,
    LOGVAR_FORCE,
    LOGVAR_STRESS,
    LOGVAR_REACTION,
    LOGVAR_RATE,
    LOGVAR_ACCEL,
    LOGVAR_ENERGY,
    LOGVAR_MASS,
	LOGVAR_CONNECTIVITY,
};



#define ERROR(mesg) do{fprintf(stderr, "[%s:%d] ", __FILE__, __LINE__);       \
perror(mesg); } while(ERROR_spin)
#define PERROR_MALLOC() { /* extern enum switch_t   writer_thread_switch;         \
extern pthread_mutex_t writer_thread_switch_lock;    \
extern pthread_cond_t  writer_thread_switch_cond; */   \
/* pthread_mutex_lock( &writer_thread_switch_lock );    \
writer_thread_switch=ON;                             \
pthread_mutex_unlock( &writer_thread_switch_lock );  \
kill( 0, SIGUSR1 );                                  \
sched_yield();                                       \
pthread_mutex_lock( &writer_thread_switch_lock );    \
while( writer_thread_switch == ON )                  \
pthread_cond_wait( &writer_thread_switch_cond,     \
&writer_thread_switch_lock );   \
pthread_mutex_unlock( &writer_thread_switch_lock ); */  \
(void) fprintf(stderr, "[%s:%d] ", __FILE__, __LINE__);\
perror( "malloc()" );  


// one logger for each logging command
class logger{
	
protected:
	bool	 active;
    double   frequency; 	// how often data must be logged
    double   lasttime;		// last time data was output to file
    double   droptol;       // below this, a zero is output
    string   logfilename;
    FILE    *logfile;		// file where this particular data is logged
    logset   settype;		// the type of the set that is logged
    int      setlabel;
	string   setName;
    void    *set;			// the actual set
    logvar   vartype;		// the type of the variable that is logged
    int      varnum;		// -1:all, >=0 the labeled
	int      varx;			// variable for gnuplot 
	int      vary;
	string   xlabel, ylabel;
	int      precision;
	
	// some private variables that are unique
	// they contain information about the existing loggers
	static  std::vector<logger*> allLoggers; 
	static  std::string summaryFilename;
	static  std::string mainlogFilename;
	static  std::string warningFilename;
	
	
	// private functions
	static  void	printHeader(ostream &os);
	
	// these are the main tasks of a logger
	virtual void	check(const model& theModel);
	virtual void	close();
	virtual void	evaluate(const double time, const model &m, const assembler& theAssembler, const energies &energies);
	virtual void	print(ostream &of=cout);
	
	
	friend class analysis;
	
public:

	// the following are not really loggers, but rather streams. It is useful to include them in the
	// logger class 
	static  ofstream sumf;
	static  ofstream mainlog;
	static  ofstream warnings;
	

	// all the clocks of the analysis
	static  timer CPUtime;		// total CPU time
	static  timer PREtime;		// preprocess time
	static  timer ANAtime;		// analysis time
	static  timer PARSEtime;	// parsing time
	static  timer POSTtime;		// time spent dumping post files
	static  time_t walltime;
	
	
	logger(const commandLine &cl);
	virtual ~logger(){}
	static void		scan(const commandLine &cl);

	// open default logger streams (not really loggers)
	static void		openLogFiles();
	static void		closeLogFiles();
	
	
	static void		checkLoggers(const model& theModel);
	static void		freeLoggersList();
	static void		listLoggers(ostream &of=cout);
	static void		logVariables(const double time, const model &m, const energies &energy, const assembler& theAssembler);

	

	// functions to open and close log file
	void OpenLogFiles(void); 
	void CloseLogFiles(void);

	// print to log file
	static void toConvFile(int step, double time, double dt, int niter, int ntrial, double res0, double resf);
	static void toLogFile(char *fmt, ...);
	static void vToLogFile(char *fmt , va_list args);
	static void toLogFile2(const char *first, const char *fmt, ...);
	static void toSummaryFile(char *fmt, ...);
	static void vToSummaryFile(char *fmt , va_list args);
	
	//time logging functions
	static void startTimers();
	static void startAnalysisTimer();	
};




class doflogger : public logger
{

public:
                    doflogger(const commandLine& cl);
                    ~doflogger(){}
	
	virtual void	evaluate(const double time, const model &m, const assembler& theAssembler, 
                             const energies &energies);
	virtual void	print(ostream &of=cout);

    
    
private:
    nodeset         privateNodeset;
	const nodeset*  P2nodeset;
	const node*     P2node;
    std::string     bodyname;
    int             dimension, manifoldLabel;
};




class energylogger : public logger
{
	
public:
	energylogger(const commandLine& cl);
	~energylogger(){}
	
	virtual void	check(const model& theModel){}
	virtual void	evaluate(const double time, const model &m, const assembler& theAssembler, const energies &energies);
	virtual void	print(ostream &of=cout);
};





class gradientlogger : public logger{
	
public:
	gradientlogger(const commandLine& cl);
	~gradientlogger(){}

	virtual void evaluate(const double time, const model &m, const assembler& theAssembler, const energies &energies);
	virtual void print(ostream &of=cout);
};



class masslogger : public logger
{
	
public:
	masslogger(const commandLine& cl);
	~masslogger(){}
	
	virtual void	check(const model& theModel){}
	virtual void	evaluate(const double time, const model &m, const assembler& theAssembler, const energies &energies);
	virtual void	print(ostream &of=cout);
};



class reactionlogger : public logger
{
	
public:
	reactionlogger(const commandLine& cl);
	~reactionlogger(){}
    
	virtual void evaluate(const double time, const model &m, 
                          const assembler& theAssembler, const energies &energies);
	
private:
    nodeset         privateNodeset;
	const nodeset*	P2nodeset;
	const node*		P2node;
    std::string     bodyname;
    int             dimension, manifoldLabel;
};



/*
class connectivityLogger : public logger
{

public:
	connectivityLogger(const commandLine& cl);
	~connectivityLogger(){}
	virtual void evaluate(const double time, const model &m, 
                          const assembler& theAssembler, const energies &energies);
	virtual void	print(ostream &of=cout);
	

private:
	bool hasbeenCalled;
};
*/


void			modifyNameForMPI(string& filename);
std::ostream&	printCentered(std::ostream&, const std::string &message);




#endif

