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
 *  dofnumberer.h
 *  feliks
 *
 *  Created by Ignacio Romero on 25/6/09.
 *
 *  The objects in this class are reponsible for filling up the IDmaps in
 *  bodies, collections of nodes, single nodes, etc... In the most general
 *  case, all the information on the connectivity of the elements in 
 *  a body/model can be used to optimize the numbering
 */


#ifndef _dofnumberer_h
#define _dofnumberer_h

#include <list>

class element;
class model;
class modelpart;
class IDmap;


class DOFnumberer
	{
	
	public:
                        DOFnumberer();
		virtual         ~DOFnumberer(){}
		
		virtual int		getLastDOFAssigned() const {return _nextDOF-1;}
		virtual void	setNextDOFtoAssign(const int n)	{_nextDOF = n;}
		virtual void	resetDOFNumbers(modelpart& b);
		virtual void	resetDOFNumbers(const std::list<modelpart*>::iterator start, 
                                        const std::list<modelpart*>::iterator end);

		virtual void	assignDOFNumbers(const std::list<modelpart*>::iterator start, 
                                         const std::list<modelpart*>::iterator end)=0;
		virtual void	assignDOFNumbers(modelpart& b)=0;
        virtual void    assignDOFNumbers(element& e) = 0;
		virtual void	fillIDs(IDmap &m);
	
	private:
		int	_nextDOF;
	};




/* The simplest node numberer just transverses the list of nodes and
 * sets the IDs in increasing order
 */
class plainDOFnumberer : public DOFnumberer
	{

	public:
		virtual void	assignDOFNumbers(const std::list<modelpart*>::iterator start, 
                                         const std::list<modelpart*>::iterator end);
		virtual void	assignDOFNumbers(modelpart& b);
        virtual void    assignDOFNumbers(element& e);
		
	protected:
		int _nextDOF;
	};




/* The dof numberer for coupled problems solved in a staggered way,
 * so that each dof group is numbered independently
 */
class splitDOFnumberer : public DOFnumberer
	{
	
	public:
		splitDOFnumberer();
		
		virtual void	assignDOFNumbers(const std::list<modelpart*>::iterator start, 
                                         const std::list<modelpart*>::iterator end);
		virtual void	assignDOFNumbers(modelpart& b);
        virtual void    assignDOFNumbers(element& e);

		void			assignUDOFNumbers(const std::list<modelpart*>::iterator start, 
                                          const std::list<modelpart*>::iterator end);
		void			assignUDOFNumbers(modelpart& b);
		
		virtual void	assignPDOFNumbers(const std::list<modelpart*>::iterator start, 
                                          const std::list<modelpart*>::iterator end);
		void			assignPDOFNumbers(modelpart& b);
	
		void			setNextUDOFtoAssign(const int n)	{_UnextDOF = n;}
		void			setNextPDOFtoAssign(const int n)	{_PnextDOF = n;}

		
	private:

		int _UnextDOF;
		int _PnextDOF;
	};


#endif

