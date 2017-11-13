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
 * element.h
 *
 * ignacio romero
 * may 2000, revised july 2003
 * rewritten dec 2011, iro to incorporate evaluationEntities
 *
 *
 * this structure and set of functions hold the functionality of the elements,
 * understood as part of the analysis. Each of them will be associated with
 * a geometrical element. This structure will be in charge of obtaining the
 * element residual, stiffness, mass, etc.
 */


#pragma once
#ifndef _feliks_element_h
#define _feliks_element_h

#include <vector>
#include <iostream>
#include <set>


#include "Elements/eltype.h"
#include "Elements/elmttasks.h"
#include "General/idata.h"
#include "Analysis/energies.h"

#ifdef WITHTBB
#include "tbb/scalable_allocator.h"
#endif



class assembler;
class evalspot;
class FEshapefunbuilder;
class integrator;
class material;
class matrix;
class modelpart;
class node;
class quadpoint;
class resultdata;

namespace feliks
{
    namespace topology
    {
        class cell;
    }
}






class element
{
	
public:
                        element(const int label_);
	virtual             ~element();

	// initializations
	void                addConnectivity(element &e1);
	virtual bool        initialize();
	void                markConstrained();
	virtual void        reverse(){};
	void                replaceNode(node& deleted, node& survivor);
    void                setColor(const int c);
	void                setParentSubmesh(modelpart &b);
    bool                testImplementation() const;

		 
	// virtual functions to update info from time step to time step
    virtual void        commitCurrentState()=0;
    virtual void        resetCurrentState()=0;
	virtual void        updateCurrentState(const dofset::evaluation_time when)=0;

    

	// functions to retrieve members of the element
	void                getNodesInFace(const int face, std::vector<node*> &nodes);
	node&               getNodeInFace(const int face, const int local) const;
	void                getNodeLabelsInFace(const int face, idata &nlabels) const;
	element*            getNeighbor(const int n) const;
	node&               getNode(const int local) const;
	const material&     getMaterial() const;

    virtual std::set<evalspot*> getEvalspots() const;

    
    
	// functions to obtain information about a Element
	double              characteristicDim();
    const feliks::topology::cell* getCellP() const;
    int                 getColor() const;
	void                getDofNumbers(vector<int>& theDofs);
    const eltype&       getEltype() const;
	geometryT           getGeometryType() const;
    int                 getLabel() const;
	modelpart&          getParentSubmesh();
	int                 getPolynomialDegree();
	virtual double      maxEigenvalue() {return 0;}	// rather max w, sqrt of max evalue = 2 / dt_elem
	size_t              getNDofs() const;
	size_t              getNDofsPerNode() const;
	size_t              nFaces() const;
	size_t              nNeighbors()	const;
	size_t              getNNodes() const;
	int                 getNNodesPerFace() const;
	void                getNodeLabels(vector<int>& labels) const;
	int                 nQuadraturePoints() const;
    double              getTime() const;
	const eltype&       getType() const;

	
	// mass matrix manipulation
	bool                hasLumpedMass();
	virtual bool        lumpMassToNodes();
	virtual	bool        lumpDampingToNodes();
	bool                getDampingMatrix(matrix& damping);


	// element status
    void                activate();
    void                deactivate();
	bool                hasConstrainedDofs() const;
    bool                isActive() const;
	bool                isVisible() const;

	
	// Functions to manipulate element data
	virtual bool        check();
	int                 defineFaces(size_t *nfaces, idata **faces);
	void                setTime(const double t);
	
    
	// functions to print element data
	virtual void        print();
	void                printGradients(ostream &of=cout) const;
	void                printTgEigenv();


    // functions integration over the element
	bool                computeNumericDDEnergy(matrix &tangent); 
	void                contributeToResidual(assembler& as);
	void                contributeToResidualTangent(assembler& as);
	virtual bool        contributeToGenericIntegral(assembler& theAssembler);

    virtual bool        integrateDampingMatrix(assembler& as);
    virtual bool        integrateEnergy(energies &ener, const dofset::evaluation_time& when) const =0;
    virtual bool        integrateDEnergy()=0;
	virtual bool        integrateMassMatrix(assembler& as);
	
	void                approximatedEigenvalues() const;
	virtual  bool       explicitForcesToNodes();
	virtual  bool       specialIntegral(assembler &as) {return false;}
	virtual  bool       sharedFunction(assembler& as) {return false;}

    virtual bool        gradient(const size_t ip, resultdata& gr, blue::ivector& gpcoor, double& dvol) const {return false;}


	// Functions to interpolate data inside element
	double interpolatePressure	(const FEshapefunbuilder& shp, const dofset::evaluation_time when=dofset::tna) const;
	double interpolatePressureRate		(const FEshapefunbuilder& shp, const dofset::evaluation_time when=dofset::tna);
	double interpolatePressureRate2		(const FEshapefunbuilder& shp, const dofset::evaluation_time when=dofset::tna);

	void interpolateAcceleration		(const FEshapefunbuilder& shp, blue::ivector &acc, const dofset::evaluation_time when=dofset::tna);
	void interpolateAngularVelocity		(const FEshapefunbuilder& shp, blue::ivector& omega);
	void interpolateBodyAlpha			(const FEshapefunbuilder& shp, blue::ivector& alpha);
	void interpolateBodyOmega			(const FEshapefunbuilder& shp, blue::ivector& avelo);
	void interpolateBodyRotationVector	(const FEshapefunbuilder& shp, blue::ivector& Theta);
	void interpolatePosition			(const FEshapefunbuilder& shp, 
                                         blue::ivector &pos, 
                                         const dofset::evaluation_time when=dofset::tna) const;
	void interpolateReferencePosition	(const FEshapefunbuilder& shp, blue::ivector &pos);
	void interpolateRotationIncrementVector(const FEshapefunbuilder& shp, blue::ivector& dtheta, blue::ivector& dthetapr);
	void interpolateRotationVector		(const FEshapefunbuilder& shp, blue::ivector& theta);
	void interpolateVelocity			(const FEshapefunbuilder& shp, 
                                         blue::ivector &vel, 
                                         const dofset::evaluation_time when=dofset::tna) const;


	
	// this functions defines an order relation by comparing labels
	inline friend bool operator<(const element &el1, const element &el2){return el1.label < el2.label;}
	inline friend bool operator==(const element &el1, const element &el2){return el1.label == el2.label;}
	static bool epLess(const element *ep1, const element *ep2);
	
	
protected:

    element(const int lbl_, const int type_,     const feliks::topology::cell* c, const int nv, node **ndlist);
	element(const int lbl_, const eltype& type_, const feliks::topology::cell* c, const int nv, node **ndlist);

    
    const eltype*                   etype;
    material*                       theMaterialDescription;
	modelpart*                      parentSubmesh;			
	energies                        theEnergies;
    const feliks::topology::cell*     theCell;
    
    size_t                          nfaces;             // hypersurfaces (sides)
    bool                            hasConstrDofs;      // true if element has a node with constrained dofs
	bool                            lumpedMass;			// true if the element mass is assumed to be lumped
	
    vector<node*>                   nodes;              // the node variables of the element
	vector<element*>                neighbors;
	
    double                          mass;               // element mass
    double                          eigmax;             // upper bound for maximum eigenvalue
	double                          storedStiffness;
    
    
private:
    
                                    element();
    bool                            active;
	int                             color;				// element partitioned color    
    int                             label;				// label assigned to the element by user
	double                          time;
    int                             typelabel;          // element type number, defined by user
};



// convenience functions
bool                staticTangent(element& el, matrix &tangent);











// inlines
inline void                 element :: activate()                   {active = true;}
inline void                 element :: deactivate()                 {active = false;}
inline int                  element :: getColor() const             {return color;}
inline const eltype&        element :: getEltype() const            {return *etype;}
inline geometryT            element :: getGeometryType() const      {return etype->geometry();}
inline int                  element :: getLabel() const             {return label;}
inline const material&      element :: getMaterial() const          {return *theMaterialDescription;}
inline node&                element :: getNode(const int l) const  {return *(nodes[l]);}
inline size_t               element :: nNeighbors()	const		{return neighbors.size();}
inline const feliks::topology::cell* element :: getCellP() const {return theCell;}
inline size_t               element :: getNNodes() const			{return nodes.size();}
inline modelpart&           element :: getParentSubmesh()			{return *parentSubmesh;}
inline double               element :: getTime() const             {return time;}
inline const eltype&        element :: getType()  const            {return *etype;}
inline bool                 element :: hasLumpedMass()				{return lumpedMass;}
inline bool                 element :: hasConstrainedDofs() const	{return hasConstrDofs;}
inline bool                 element :: isActive() const			{return active;}
inline bool                 element :: isVisible() const			{return etype->isVisible();}
inline void                 element :: setColor(const int c)       {color = c;}
inline void                 element :: setTime(const double t)		{time = t;}
inline void                 element :: setParentSubmesh(modelpart &b){parentSubmesh = &b;}


#endif

