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
 *  tscontrol.h
 *
 *  Luis Lacoma        	March-2005
 *	  				Revised June-2005
 *
 */

#ifndef _tscontrol_h
#define _tscontrol_h

#include <string>

#include "Model/model.h"
#include "Io/usercommand.h"
#include "Analysis/energies.h"
#include "Io/usercommand.h"


class tscontrol
{     

public:
	double   dt;
	double   mindt;							// dt below which auto stepping stops
	double   maxdt;							// dt above which auto stepping stops
		
	tscontrol(commandLine &cl);
	tscontrol();
	virtual ~tscontrol(){};

	static void		scan(commandLine &cl, tscontrol **t);
	void			setLinks(energies &ener, model &m);
	
	
	double			adaptiveForwardDt();
	void			advance();
	void			advanceDtVariables();
	bool			check(std::ostream& of=std::cout);
	double			computeIntegrationTolerance(model &m);
	double			decreaseDt();

	inline void		setMaxDt(double mdt)		{mindt = mdt;}
	inline void		setMinDt(double mdt)		{maxdt = mdt;}
	inline void		setDt(double dtc)			{dt = dtc;}

	
	virtual double  chooseNewDt(double dt, double error, double tol, double *dtP, double *errorP, double *tolP, double *p)=0;
	inline double   getError(){return error;}
	inline  double  getPreviousDt(){return dtHistory[0];}
	inline  std::string  getName(){return name;}
	virtual void    info(std::ostream &of=std::cout);
	void            setParameters();
	void            setData(double tol, int nriter);


protected:
	std::string     name;
	double          error;					// error employed for time step selection
	double          tol;
	double          dtHistory[2];			// vector with time steps at t_{n+1}, t_n and t_{n-1}
	double          errorHistory[2];		// vector with errors at t_{n+1}, t_n and t_{n-1}
	double          tolHistory[2];			// vector with tolerances at t_{n+1}, t_n and t_{n-1}
	double          parameter[3];			// vector with parameters of time step controller
	int             nSteps;					// number of accepted steps
	double          orderMethod;			// order of convergence of the integration method. Usually, it is 2
	int             rejectedSteps;			// number of steps which have been rejected
	double          userTolerance;			// tolerance factor chosen by user
	energies        *energy;
	model           *linkedModel;
    
	
    

};



class fixedTSC : public tscontrol
{
	public:
	fixedTSC(commandLine &cl);
	fixedTSC(double dt=1.0);

	double chooseNewDt(double dt, double error, double tol, double *dtP, double *errorP, double *tolP, double *p);
	void info(std::ostream &of=std::cout);
};



class iterationTSC : public tscontrol
{
	private:
		int     targetIterations;

	
	public:
	iterationTSC(commandLine &cl);
	
	double chooseNewDt(double dt, double error, double tol, double *dtP, double *errorP, double *tolP, double *p);
	void info(std::ostream &of=std::cout);
};



class iTSC : public tscontrol
{
	private:
		double p0;
	
	public:
		iTSC(commandLine &cl);
	
		double chooseNewDt(double dt, double error, double tol, double *dtP, double *errorP, double *tolP, double *p);
		void   info(std::ostream &of=std::cout);
};



class cflTSC : public tscontrol
{
private:
	double	maxev;
	int		computeEvery;
	int		counter;
	double  scaletmin;
	
public:
	cflTSC(commandLine &cl);
	cflTSC();

	
	double chooseNewDt(double dt, double error, double tol, double *dtP, double *errorP, double *tolP, double *p);
	void   info(std::ostream &of=std::cout);
};


#endif



