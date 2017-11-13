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
 * propfactor.c
 *
 * i. romero, december 2003
 *
 * as explained in propfactor.h the proportional factors f_i(t) (propfactors) for the
 * loading and imposed displacements are grouped in propfactor combination. Every
 * load or imposed displacement is multiplied by the value of the corresponding
 * propfactor combination, defined as
 *
 *                  F(t) = Sum_{i=0}^Npropfactors  f_i(t)
 *
 * where Npropfactors is the number of proportional factors in the combination. All the
 * combinations defined are stored in a table called prflist. Each of the entries
 * in prflist is a propfactor combination, i.e., an array of propfactors. The array
 * prflistindex keeps track of how many propfactors are stored in each combination.
 */


#include <list>
#include <map>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>
#include <sstream>
#include "boost/foreach.hpp"

#include "Analysis/Propfactors/propfactor.h"
#include "Io/message.h"
#include "Io/usercommand.h"
#include "Math/tensor.h"
#include "Io/logger.h"


using namespace blue;


// initialization of static variables
map<int,factorcombo*>	factorcombo :: combos;


factorcombo :: factorcombo(int labelc) :
label(labelc),
factors(0)
{
}


/*
 factorcombo :: factorcombo(int labelc, commandLine *clist, int nfactorsc) :
 value(0.0),
 label(labelc),
 lastEvalTime(-1.0),
 factors(0),
 print(false)
 {
 cout << "is this function ever used??" << endl;
 
 // loop over all the proportional factors sent
 propfactor *prf=NULL;
 for (int i=0; i< nfactorsc; i++)
 {
 // read proportional factor in this combination
 const usercommand &uc = clist[i][0];
 
 if ( uc.keyword() != "type")
 Message("\n Error scanning propfactor %d. 'type' must be the first keyword", label);
 else
 {
 if      (uc.option() == "constant")  prf = new constantprf(clist[i]);
 else if (uc.option() == "linear")    prf = new linearprf(clist[i]);
 else if (uc.option() == "triangle")  prf = new triangleprf(clist[i]);
 else if (uc.option() == "sine")      prf = new sineprf(clist[i]);
 
 if (prf != NULL)
 factors.push_back(prf);
 else
 Message("\bn Error scanning propfactor %d. Type %s is not defined", label, (uc.option).c_str());
 }
 }
 }
 */




/* commented by iro 31/08/2006
 // add a propfactor combination to the static linked list
 void factorcombo :: add(factorcombo &fc)
 {
 static bool replaced0 = false;
 
 // if the combination labeled 0 has not been defined, it must
 if ( combos.find(0) == combos.end() )  initializePropfactorCombinations();
 
 // if fc has label '0' this means it is a replacement for
 // the default proportional factor combination, the one
 // which sits in the first slot of the linked list.
 // In this case, do not add, but replace, but only if it is the first
 // replacement. Further additions only add more scaling factors
 if (fc.label == 0 || !replaced0)
 {
 replace0 = true;
 combos.erase(0);
 combos.insert( make_pair(fc.label, &fc) );
 }
 
 // no propfactor combination can have negative label
 else if (fc.label < 0)
 {
 Message("\n Warning: Proportional factors can not have negative labels.");
 delete &fc;
 }
 
 // in the rest of cases, add the propfactor combo to the linked list
 else
 combos.insert( make_pair(fc.label, &fc) );
 }
 */


void factorcombo :: add(int combolabel, propfactor &p)
{
	static bool replaced0=false;
	
	// if the combination labeled 0 has not been defined, it must
    if ( combos.find(0) == combos.end() )  initializePropfactorCombinations();
    
	// if the combination to which this factor must be added has never been
	// used before, create it
	factorcombo *combo = NULL;
	if ( combos.find(combolabel) == combos.end() ) combo = new factorcombo(combolabel);
	else  		                                   combo = (combos.find(combolabel))->second;
    
	// if fc has label '0' this means it is a replacement for
    // the default proportional factor combination, the one
    // which sits in the first slot of the linked list.
    // In this case, do not add, but replace, but only if it is the first
	// replacement. Further additions only add more scaling factors
    if (combolabel == 0 && !replaced0)
    {
		replaced0 = true;
		combos.erase(0);
		combo = new factorcombo(0);
		combos.insert( make_pair(combolabel, combo) );
	}
	
	// add the propfactor to the (newly created) combination
	combo->factors.push_back(&p);
	combos.insert( make_pair(combolabel, combo) );
}








/* evaluates, at time t, the propfactor combination number ncombo, as stored in
 the global list. It also stores in its variable 'value' the evaluation result,
 for reuse.
 */
double factorcombo :: eval(const double t, const ivector& p) const
{
    double value = 0.0;
	list<propfactor*>::const_iterator iter = factors.begin();
	while (iter != factors.end())
	{
		value += (*iter)->eval(t, p);
		++iter;
	}
    
	return value;
}




double factorcombo :: eval(const int label, const double t, const ivector& p)
{
    factorcombo &fc = getPropfactorCombinationFromLabel(label);
	return fc.eval(t, p);
}




factorcombo :: ~factorcombo()
{
	// before erasing the list itself we must empty the memory it contains
	list<propfactor*>:: iterator iter = factors.begin();
	while(iter != factors.end() )
	{
		delete *iter;
		++iter;
	}
}




/* this function deallocates all the memory used for all combinations.
 * It should be called once, at the very end.
 */
void factorcombo ::  freePropfactorCombinations()
{
	map<int,factorcombo*>::iterator iter = combos.begin();
	
    while (iter != combos.end() )
    {
        delete iter->second;
		++iter;
    }
}



/* traverse the linked list to see if one of the stored factors has
 * a label that coincides with the asked one.
 * Many times the same label is asked repeatedly, so it pays off to
 * save the last requested label and result, just in case the same one
 * is asked for.
 */
factorcombo& factorcombo ::  getPropfactorCombinationFromLabel(int label)
{
	map<int, factorcombo*>:: iterator  iter = combos.find(label);
    
	// propfactor with label l not found
    if ( iter == combos.end() )
        ErrorMessage("\n Warning: scaling factor %d is not defined.", label);
    
	return *(iter->second);
}




/* initialized the linked list of proportional factor combinations, so that
 in the first one there is a linear factor
 */
void factorcombo :: initializePropfactorCombinations()
{
	if ( combos.find(0) == combos.end() )
    {
		propfactor  *prf        = new linearprf(0.0, HUGE, 0, 1.0);
        factorcombo *firstcombo = new factorcombo(0);
		firstcombo->factors.push_front(prf);
		combos.insert( make_pair(0, firstcombo) );
    }
}




void factorcombo :: print(ostream &of) const
{
    of  << endl << endl
    << " Scaling function number " << label
    << " (" << factors.size() << " components)"
    << flush;
	
	list<propfactor*>::const_iterator iter = factors.begin();
	while(iter != factors.end() )
	{
		(*iter)->print(of);
		++iter;
	}
}




void factorcombo ::	printAllValues(const double time, ostream &of)
{
	ivector dummy;
	map<int,factorcombo*>::iterator iter = combos.begin();
	while (iter != combos.end() )
	{
		of  << "\n  Scaling factor " << setw(3) << iter->first
        << " : " << scientific << setw(12) << iter->second->eval(time, dummy);
		++iter;
	}
}




void factorcombo :: listPropfactorCombinations(ostream &of)
{
    of  << "\n\n\n\n";
	stringstream title;
	title << "P r o p o r t i o n a l   f a c t o r s  (" << combos.size() << ")";
	string str(title.str());
	printCentered(of, str);
    
	map<int,factorcombo*>::iterator iter = combos.begin();
	while (iter != combos.end() )
	{
		iter->second->print(of);
		++iter;
	}
    of  << endl;
}



void factorcombo :: scan(const commandLine &cl)
{
	propfactor  *p=NULL;
	int          combolabel(0);
	
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if      (uc.keyword() == "combination")  combolabel = (int) uc.value();
		else if (uc.keyword() == "label")		combolabel = (int) uc.value();
		else if (uc.keyword() == "type" && uc.option() == "constant") p = new constantprf(cl);
		else if (uc.keyword() == "type" && uc.option() == "linear")   p = new linearprf(cl);
		else if (uc.keyword() == "type" && uc.option() == "triangle") p = new triangleprf(cl);
		else if (uc.keyword() == "type" && uc.option() == "sine")     p = new sineprf(cl);
	}
	if (p != NULL) add(combolabel, *p);
}




// a constructor of all propfactors, just to set up default values
propfactor :: propfactor() :
start(0.0),
end(HUGE),
repeat(0.0),
printme(false)
{
}




// a constructor of all propfactors, just to set up default values
propfactor :: propfactor(const commandLine& cl) :
start(0.0),
end(HUGE),
repeat(0.0),
printme(false)
{
    
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
        
		// common options for all propfactors
		if      (uc.keyword() == "start")  start  = uc.value();
		else if (uc.keyword() == "end")    end    = uc.value();
		else if (uc.keyword() == "print" && uc.option() == "yes" ) printme = true;
	}
}




/* create a proportional factor of the type CONSTANT. The function is
 * f(t) =    value    ,   start <= time <= end
 *           0        ,   otherwise
 */
constantprf :: constantprf(const commandLine& cl) :
propfactor(cl),
c(1.0)
{
   for (int i=1; i< cl.size(); i++)
    {
        // scan for the word "type"
        const usercommand &uc = cl[i];
		
        if      (uc.keyword() == "value")  c = uc.value();
    }
}




double constantprf :: eval(const double t, const ivector& p) const
{
	if (start < t && t <= end)  return c;
	else                        return 0.0;
}



void constantprf :: print(ostream &of) const
{
	of << "\n  +constant pr. factor active in (" << start << ", ";
	if   (end > 0.9*HUGE) of << "infinite]";
	else                  of << end << "]";
	of << "\n   value = " << c << flush;
}




linearprf :: linearprf(const commandLine& cl) :
propfactor(cl),
init(0.0),
slope(1.0)
{
    for (int i=1; i< cl.size(); i++)
    {
        // scan for the word "type"
        const usercommand &uc = cl[i];
		
        if (uc.keyword() == "initial")	init   = uc.value();
        else if (uc.keyword() == "slope")	slope  = uc.value();
    }
}




linearprf :: linearprf(double startc, double endc, double initc, double slopec) :
init(initc),
slope(slopec)
{
	start = startc;
	end   = endc;
}




double linearprf :: eval(const double t, const ivector& p) const
{
	if (start < t && t <= end)  return init + (t-start)*slope;
	else                        return 0.0;
}




void linearprf :: print(ostream &of) const
{
	of  << "\n  +linear pr. factor active in (" << start << ", ";
    if   (end > 0.9*HUGE) of << "infinite]";
    else                  of << end << "]";
	of << "\n   initial value = " << init << ",  slope = " << slope;
}





triangleprf :: triangleprf(const commandLine& cl) :
propfactor(cl),
peak(-1.0),
height(1.0)
{
    for (int i=1; i< cl.size(); i++)
    {
        // scan for the word "type"
        const usercommand &uc = cl[i];
		
        if      (uc.keyword() == "peak")   peak   = uc.value();
        else if (uc.keyword() == "height") height = uc.value();
        else if (uc.keyword() == "repeat") repeat = uc.value();
    }
    
	if (peak == -1.0) peak = 0.5*(start+end);
}




double triangleprf :: eval(const double t, const ivector& p) const
{
	if      (start < t && t <= peak) return (t-start)*height/(peak-start);
	else if (peak  < t && t <= end)  return (t-end)  *height/(peak-end);
	else							 return  0.0;
}




void triangleprf :: print(ostream &of) const
{
	of  << "\n  +triangular pr. factor active in (" << start << ", ";
	if   (end > 0.9*HUGE) of << "infinite]";
	else                  of << end << "]";
	of << "\n   peak time = " << peak << ",  height = " << height << flush;
}



sineprf :: sineprf(const commandLine& cl) :
propfactor(cl),
amplitude(0.0),
frequency(0.0),
period(0.0),
phase(0.0)
{
	for (int i=1; i< cl.size(); i++)
    {
        // scan for the word "type"
        const usercommand &uc = cl[i];
		
		if      (uc.keyword() == "amplitude")	amplitude   = uc.value();
        else if (uc.keyword() == "frequency")	frequency	= uc.value();
        else if (uc.keyword() == "period")		period		= uc.value();
        else if (uc.keyword() == "phase")		phase		= uc.value();
	}
	
	if (frequency == 0 && period == 0)
	{
        Message("\n Error in the definition of propfactor of type sine.");
        Message("\n Either the frequency or the period must be defined.");
		frequency = 1.0;
	}
	
	// either the period or the frequency is > 0 (or both)
	if (period > 0.0) frequency = 1.0/period;
	else              period    = 1.0/frequency;
	
}




double sineprf :: eval(const double t, const ivector& p) const
{
	double w = 2.0*M_PI*frequency;
	
	if (start < t && t<= end) return amplitude*sin(w*(t-start)+phase);
	else					  return 0.0;
}




void sineprf :: print(ostream &of) const
{
	of  << "\n  +sine pr. factor active in (" << start << ", ";
	if   (end > 0.9*HUGE) of << "infinite]";
	else                  of << end << "]";
    
	of << "\n   amplitude = " << amplitude << ", frequency = " << frequency << " Hz," << flush;
	of << "period = " << period << ", phase shift = " << phase << flush;
}



