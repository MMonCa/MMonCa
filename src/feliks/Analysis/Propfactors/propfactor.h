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
 * propfactor.h
 *
 * FELIKS, a general purspose finite element program
 * c I. Romero, Universidad Politecnica de Madrid.
 *
 * i.romero, march 2002, december 2003
 *
 * proportional factor functions. Each "force" or imposed displacement
 * defined in the model is multiplied by a scalar function during the analysis 
 * so that loads and imposed displacements can be constant, change linearly,
 * exponentially, ... or with any linear combination of these.
 * 
 * This scalar function is made of the sum of several scalar factors. We
 * refer to this sum as the proportional factor. Each of these prfactors can
 * be made of as many individual functions as the user defines and it
 * has to be associated with a label that identifies it uniquely.
 *
 * In the analysis definition file, the user says something similar to:
 *
 *         scaling, label = 7, type = linear, start = 3.1 , end = 10.6 , slope = 3.14
 *         scaling, label = 7, type = sine  , start = 0.0 , ampl = 7.0 , period = 1.0
 *         scaling, label = 7, type = constant, start = 3.0, end = 12.0, value = -23.9
 *
 * Then, in the model file, for each imposed load or boundary condition the user
 * defines
 *
 *         bc, nodeset = ... , scaling = 7, ux = 3.0
 *
 * which means that the 1st dof of all the nodes in the nodeset have an imposed value of 14.3 x prop_3(t),
 *
 *
 * Every load and bc defined in the model must have an associated scaling factor.
 * If no proportional value is defined, the default value 0 is assigned. If the
 * proportional factor 0 is not defined by the user, it is implicitely defined
 * as type=linear, start=0, end=infty, slope = 1.0
 *
 * 
 */



#ifndef _propfactor_h
#define _propfactor_h


#include <list>
#include <map>
#include <iostream>

#include "Math/tensor.h"

#ifdef WITHTBB
#include "tbb/spin_mutex.h"
#endif

class commandLine;



class propfactor
{

public:
                        propfactor();
                        propfactor(const commandLine& cl);
	virtual             ~propfactor(){};
	
	virtual double      eval(const double t, const blue::ivector& p) const =0;
	virtual void        print(std::ostream &of=std::cout) const =0;
	
    
    
protected:
	bool                printme;		// display message of evaluation
    double              start;			// start time for the function f(t)
    double              end;			// same, but end
    double              repeat;			// for repeating loads
};




// propfactor combination
class factorcombo
{
	
private:
	int                     label;		// identifying numeric label
    std::list<propfactor*>  factors;	// must be pointers because they can be constant, linear, sine ...
    
    
	// map with the factorcombos associated with label
	static std::map<int,factorcombo*> combos;
	
	
public:
	factorcombo(int labelc);
	factorcombo(int labelc, commandLine *clist, int nfactorsc);
	~factorcombo();

	static void freePropfactorCombinations();
	static void scan(const commandLine &cl);

	//This should be called when the analysis is initialized
	static void initializePropfactorCombinations();

	// this is used to add more propfactor combinations to the linked list
	// commented iro 31/8/2006 static void add(factorcombo &fc);
	static void add(int combolabel, propfactor &p);

	/* this are the main evaluation function. When we want to obtain the value of all
		the propfactor combinations at one instant t, we first call EvalAllPropfactorCombinations(),
		and then, one by one with PropfactorCombinationValue().
		Alternatively, one could evaluate one by one each PropfactorCombination with the functions
		described next. The preferred strategy is this first one, for it saves time to evaluate
		only once each propfactorcombination */
	static double eval(const int label, const double t, const blue::ivector& p);

	double PropfactorCombinationValue(factorcombo fc);
	double eval(const double t, const blue::ivector& p) const;

	// return a proportional factor combination given its label
	static factorcombo& getPropfactorCombinationFromLabel(int label);
	
	
	// functions to print the factors
	void			print(std::ostream &of=std::cout) const;
	static void		printAllValues(const double time, std::ostream &of=std::cout);
	static void		listPropfactorCombinations(std::ostream &of=std::cout);
	int				getLabel() const {return label;}
};



class constantprf : public propfactor
{
	private:
		double c;
		
	public:
		constantprf(const commandLine& cl);
	
		double eval(const double t, const blue::ivector& p) const;
		void   print(std::ostream &of=std::cout) const;
};




class linearprf : public propfactor
{
private:
	double init, slope;
	
public:
	linearprf(const commandLine& cl);
	linearprf(double startc, double endc, double initc, double slopec);
    
	double eval(const double t, const blue::ivector& p) const;
	void   print(std::ostream &of=std::cout) const;
};



class triangleprf : public propfactor
{
private:
	double peak, height;
	
public:
	triangleprf(const commandLine& cl);
    
	double eval(const double t, const blue::ivector& p) const;
	void   print(std::ostream &of=std::cout) const;
};



class sineprf : public propfactor
{
private:
	double amplitude, frequency, period, phase;
	
public:
	sineprf(const commandLine& cl);
    
	double eval(const double t, const blue::ivector& p) const;
	void   print(std::ostream &of=std::cout) const;
};




#endif
