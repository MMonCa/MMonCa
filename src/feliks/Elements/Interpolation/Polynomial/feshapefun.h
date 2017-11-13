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
 * shapef.h
 *
 * shape functions routines
 *
 * i. romero , feb 2002
 */


#ifndef _shapef_h
#define _shapef_h
	
	
#include "Geometry/quadpoint.h"
#include "Elements/Interpolation/shapefun.h"
#include "Math/tensor.h"
#include "Main/feliks.h"
#include <vector>
#include <iostream>

#ifdef WITHTBB
#include "tbb/scalable_allocator.h"
#endif


#define PUSH      +1
#define NOPUSH     0
#define SECONDD   +1
#define NOSECONDD  0


class element;
class elmtface;
class node;
class cell;




class FEshapefunbuilder
{
		
public:
		
    // the only way to initialize a FEshapefunbuilder is to associate it to an element
                        FEshapefunbuilder(const element& e);
                        FEshapefunbuilder(const element& e, const int nshp);
	virtual             ~FEshapefunbuilder(){}
	FEshapefunbuilder&	operator=(const FEshapefunbuilder &sh);
	
    size_t              size() const;
    const shapefun&     operator[](const size_t k) const;
    shapefun&           operator[](const size_t k);
    vector<shapefun>&   getTheShapefunctions(){return shp;}
    

	
	// evaluates all the shape functions and their derivatives at a gauss point, 
	// with respect to the current coordinates. Also
	//	the jacobian of the transformation. 
	void evaluate(const quadpoint& qp, double& j, const bool dd=false, const int push=PUSH);
		
	// the same as above but the derivatives are w/r to the isoparametric coordinates
	virtual void evaluateRef(const quadpoint& gp);
	
	// the same as above but also the second derivatives of the shape functions
	virtual void evaluateDD(const quadpoint& gp, double& jacobian);



		
    
protected:
	vector<shapefun>    shp;
	void            eval1d(		const quadpoint& gp, double& j, const bool dd, const int push);
	void            eval2d(		const quadpoint& gp, double& j, const bool dd, const int push);
	void            eval3d(		const quadpoint& gp, double& j, const bool dd, const int push);
	void            evalCST(	const quadpoint& gp, double& j, const bool dd, const int push);
	void            evalQUAD(	const quadpoint& gp, double& j, const bool dd, const int push);
	void            evalTRIQUAD(const quadpoint& gp, double& j, const bool dd, const int push);
	void            evalBRICK(	const quadpoint& gp, double& j, const bool dd, const int push);
	void            evalTET(	const quadpoint& gp, double& j, const bool dd, const int push);
	void            evalTET10(	const quadpoint& gp, double& j, const bool dd, const int push);
    FEshapefunbuilder();
};


inline size_t               FEshapefunbuilder::size() const {return shp.size();}
inline const shapefun&      FEshapefunbuilder::operator[](const size_t k) const {return shp[k];}
inline shapefun&            FEshapefunbuilder::operator[](const size_t k)		{return shp[k];}
inline void                 FEshapefunbuilder::evaluate(const quadpoint& qp, double& j, const bool dd, const int push)
                                {eval3d(qp, j, dd, push);}




// These are shape functions defined over a surface in a 3D problem. They are not shape functions
// for plane problem.
class surfaceShapefuns : public FEshapefunbuilder
{
	
public:
	surfaceShapefuns(const element& e);
	surfaceShapefuns(const element& e, const int nshp);
	virtual ~surfaceShapefuns(){}

	void evaluate(		const quadpoint& qp, double& j);
	void evaluateRef(	const quadpoint& gp);
	void evaluateOrthogonal(const blue::ivector& a0, const blue::ivector& a1, surfaceShapefuns& shpx, double& sj);
	
private:
	surfaceShapefuns();
	
};




class faceShapefuns : public FEshapefunbuilder
{

public:
	faceShapefuns(const elmtface& aFace);
	virtual ~faceShapefuns(){}
	
	void evaluate(const quadpoint& qp, double& j);
	void evaluateRef(const quadpoint& gp);
	void evaluateOrthogonal(const blue::ivector& a0, const blue::ivector& a1, faceShapefuns& shpx, double& sj);
	
private:
	faceShapefuns();	
};




// These are shape functions defined over a curve in a 1, 2, or 3Dproblem. 
class curveShapefuns : public FEshapefunbuilder{
	
public:
	curveShapefuns(const element& e);
	virtual ~curveShapefuns(){}
	
	void evaluate(		const quadpoint& gp, double& jacobian);
	void evaluateRef(	const quadpoint& gp);
	
private:
	curveShapefuns();
};



#endif


