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
 *  loading.h
 *  feliks
 *
 *  Created by Ignacio Romero on 4/5/10.
 *  Copyright 2010 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#ifndef _loading_h
#define _loading_h


#include <iostream>
#include "Math/tensor.h"


class assembler;
class commandLine;
class energies;
class factorcombo;
class longvector;
class node;
class model;


class loading
{


public:
	enum loadingT
    {
		LOAD_POINT,
		LOAD_SURFACE,
		LOAD_BODY,
		LOAD_BLAST,
        LOAD_NODESET
	};	


                    loading(const int propfactorlabel);
                    loading(const commandLine& cl);
	virtual         ~loading(){}
	
	virtual void	accumulateEnergy(energies& ener, const double time)=0;
	virtual void	assembleLoading(assembler& asb, const double time, const double f)=0;
	virtual void	transferLoadsToNodes(const double time)=0;
	void			replaceNode(node& deleted, node& survivor){}

    factorcombo&	getScalingFactor();
	int				getPropfactorLabel() const {return factorlabel;}
	loadingT		getType() const {return theType;}
	virtual void	initialize(model& theModel);
	virtual void	print(std::ostream& of=std::cout) const;
	

	static void		scan(const commandLine &cl, model& theModel);

	
protected:
	factorcombo		*factor;			// a pointer to the proportional factor
	loadingT		theType;
	int				factorlabel;		// the label of the proportional factor
	blue::ivector			fU;					// force on the U dof
	blue::ivector			fR;					// force on the R dof
	double			fP;					// force on the P dof
	blue::ivector			fD;					// force on the D dof
	
	bool			uloaded[3], rloaded[3], dloaded[2], ploaded;
	
private:
	loading();	
};



#endif

