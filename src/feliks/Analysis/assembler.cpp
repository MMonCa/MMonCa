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
 * assembler.cpp
 *
 * i. romero
 * july 2000, last modification december 2003
 * march 2006, translated to c++
 *
 * these are the functions responsible for assembling the element contributions
 * into the global residual vector and tangent matrix
 */


#include "Analysis/assembler.h"

#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include "boost/foreach.hpp"

#include "Analysis/analysis.h"
#include "Analysis/Dofs/dofset.h"
#include "Analysis/Loading/pointload.h"
#include "Analysis/Propfactors/propfactor.h"
#include "Elements/evalspot.h"
#include "Elements/element.h"

#include "General/idata.h"
#include "General/feliksutil.h"

#include "Io/message.h"
#include "Main/feliks.h"
#include "Math/matrix.h"
#include "Model/Parts/modelpart.h"
#include "Model/model.h"
#include "Model/Node/node.h"
#include "Math/Sparse/sparse.h"
#include "Math/Sparse/petscmatrix.h"
#include "Model/Parts/modelpart.h"

#ifdef WITHTBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#endif


using namespace blue;


assembler :: assembler() :
	_v(0), _m(0) ,
	_first_operator_split(false),
	_second_operator_split(false),
	_time_derivatives(0)
{
	_tasks.reset();
}

assembler :: assembler(longvector& v, sparseMatrix& m) :
	_v(&v), _m(&m)  ,
	_first_operator_split(false),
	_second_operator_split(false),
	_time_derivatives(0)
{
	_tasks.reset();
}


void assembler :: assemble	(const node& nd, const int dof, const double f)
{
#ifdef WITHTBB
	tbb::spin_mutex::scoped_lock my_lock;
	my_lock.acquire(const_cast<vectorDofset&>(nd.getUDS()).assemble_mutex);
#endif
	
	if ( !nd.isDofConstrained(dof) ) _v->data[ nd.getID(dof)] += f;
	
#ifdef WITHTBB
	my_lock.release();
#endif
	
}


void assembler :: assemble	(const ivector& f, const vectorDofset& uds)
{
#ifdef WITHTBB
	tbb::spin_mutex::scoped_lock my_lock;
	my_lock.acquire(const_cast<vectorDofset&>(uds).assemble_mutex);
#endif
	if ( ! uds.isConstrained(0) ) _v->data[ uds.getID(0) ] += f[0];
	if ( ! uds.isConstrained(1) ) _v->data[ uds.getID(1) ] += f[1];
	if ( ! uds.isConstrained(2) ) _v->data[ uds.getID(2) ] += f[2];	
	
#ifdef WITHTBB
	my_lock.release();
#endif
}


void assembler :: assemble	(const ivector& f, const directorDofset& uds)
{
#ifdef WITHTBB
	tbb::spin_mutex::scoped_lock my_lock;
	my_lock.acquire(const_cast<directorDofset&>(uds).assemble_mutex);
#endif
	if ( ! uds.isConstrained(0) ) _v->data[ uds.getID(0) ] += f[0];
	if ( ! uds.isConstrained(1) ) _v->data[ uds.getID(1) ] += f[1];	
#ifdef WITHTBB
	my_lock.release();
#endif
}



void assembler :: assemble	(const ivector& f, const rotationDofset& uds)
{
#ifdef WITHTBB
	tbb::spin_mutex::scoped_lock my_lock;
	my_lock.acquire(const_cast<rotationDofset&>(uds).assemble_mutex);
#endif
	
	if ( ! uds.isConstrained(0) ) _v->data[ uds.getID(0) ] += f[0];
	if ( ! uds.isConstrained(1) ) _v->data[ uds.getID(1) ] += f[1];	
	if ( ! uds.isConstrained(2) ) _v->data[ uds.getID(2) ] += f[2];	

#ifdef WITHTBB
	my_lock.release();
#endif
}



void assembler :: assemble	(const itensor& t, vectorDofset& uda, vectorDofset& udb)
{
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			int IDi = uda.getID(i);
			for (size_t j=0; j<3; j++)
				if ( ! udb.isConstrained(j) ) _m->addTo( t(i,j), IDi, udb.getID(j) );
		}
	}
}


void assembler :: assemble	(const double f, vectorDofset& uda, vectorDofset& udb)
{
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) && ! udb.isConstrained(i) ) 
			_m->addTo( f, uda.getID(i), udb.getID(i) );
	}
}


void assembler :: assemble	(const itensor& t, vectorDofset& uda, rotationDofset& udb)
{
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( ! udb.isConstrained(j) ) _m->addTo( t(i,j), uda.getID(i), udb.getID(j) );
		}
	}
}

void assembler :: assemble	(const itensor& t, rotationDofset& uda, vectorDofset& udb)
{
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( ! udb.isConstrained(j) ) _m->addTo( t(i,j), uda.getID(i), udb.getID(j) );
		}
	}
}

void assembler :: assemble	(const itensor& t, rotationDofset& uda, rotationDofset& udb)
{
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( ! udb.isConstrained(j) ) _m->addTo( t(i,j), uda.getID(i), udb.getID(j) );
		}
	}
}



void assembler :: assemble	(const itensor& t, vectorDofset& uda, directorDofset& udb)
{
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<2; j++)
				if ( ! udb.isConstrained(j) ) _m->addTo( t(i,j), uda.getID(i), udb.getID(j) );
		}
	}
}



void assembler :: assemble	(const itensor& t, directorDofset& uda, vectorDofset& udb)
{
	for (size_t i=0; i<2; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( ! udb.isConstrained(j) ) _m->addTo( t(i,j), uda.getID(i), udb.getID(j) );
		}
	}
}



void assembler :: assemble	(const itensor& t, directorDofset& uda, directorDofset& udb)
{
	for (size_t i=0; i<2; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<2; j++)
				if ( ! udb.isConstrained(j) ) _m->addTo( t(i,j), uda.getID(i), udb.getID(j) );
		}
	}
}



void assembler :: assemble	(const itensor& t, Unode& nda, Unode& ndb)
{
	vectorDofset& uda = nda.getUDS();
	vectorDofset& udb = ndb.getUDS();
	
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( ! udb.isConstrained(j) ) _m->addTo( t(i,j), uda.getID(i), udb.getID(j) );
		}
	}
}



void assembler :: assemble	(const ivector& v, vectorDofset& dsa, scalarDofset& dsb)
{
	for (size_t i=0; i<3; i++)
		if (! dsa.isConstrained(i) ) 
			if ( ! dsb.isConstrained() ) 
				_m->addTo( v[i], dsa.getID(i), dsb.getID() );
}


void assembler :: assemble	(const ivector& v, scalarDofset& dsa, vectorDofset& dsb)
{
	if (! dsa.isConstrained() ) 
		for (size_t j=0; j<3; j++)
			if ( ! dsb.isConstrained(j) ) 
				_m->addTo( v[j], dsa.getID(0), dsb.getID(j) );
}


void assembler :: assemble	(const double d, scalarDofset& dsa, scalarDofset& dsb)
{
	if (! dsa.isConstrained() ) 
		if ( ! dsb.isConstrained() ) 
			_m->addTo( d , dsa.getID(), dsb.getID() );
}



void assembler :: assemble			(const double f, const scalarDofset& sds)
{
	if (! sds.isConstrained() ) 
		_v->data[ sds.getID(0) ] += f;
}



void assembler ::  assemble(const longvector& v, const vector<int>& doflabels)
{
	_v->assembleInVector(v, doflabels);
}



void assembler ::  assemble(const matrix& m, const vector<int>& doflabels)
{
	_m->assemble(m, doflabels);
}




/* Loops over all the elements of the analysis model, requesting their element
 damping matrices and adding their contribution to the global matrix. 
 */
void assembler :: assembleDampingMatrix(model &m, sparseMatrix &damping)
{
    // Reset global mass matrix
    damping.setZero();
	
	// Define assembler tasks
	_tasks.reset();
	_tasks.mass_matrix = true;
	
	// loop through 
    BOOST_FOREACH(modelpart* p , m.theParts )
        p->integrateDampingMatrix(*this);
	

    _m->finalizeAssembly();
    DebugMessage("Consistent damping assembled.");
}





/* Add external forces to global residual R = fact*Fext - Fint */
void assembler :: assembleExternalForces(const std::list<loading*>& loads, double time, double fact)
{
    if (time < 0.0 || fact <= 0.0) return;	
	
    BOOST_FOREACH(loading* ld, loads)
        ld->assembleLoading(*this, time, fact);

    DebugMessage("External forces incorporated to residual.");	
}




void assembler :: assembleGenericContribution(modelpart& sm)
{
	if ( !sm.isActive() ) return;
	
    BOOST_FOREACH(evalspot* ev, sm.evalspots )
		ev->genericIntegral(*this);
    
    BOOST_FOREACH(element* el, sm.elements)
        el->contributeToGenericIntegral(*this);

    _m->finalizeAssembly();
    DebugMessage("Generic integration performed.");
}	



void assembler :: assembleL2ProjectionMatrix(const model &m)
{
	// Define assembler tasks
	_tasks.reset();
	_tasks.L2_projection_matrix = true;
	
    // Reset global tangent
    _m->setZero();
	
    // indicate that, if the tangent has the LU factors already computed, they are obsolete
    _m->setNotFactorized();
    
    // loop through bodies
    BOOST_FOREACH(modelpart* p, m.theParts )
        assembleGenericContribution( *p );

	
	_m->finalizeAssembly();
	DebugMessage("L2 projection matrix assembled.");
}



/* Loops over all the elements of the analysis model, requesting their element
 mass matrices and adding their contribution to the global matrix. 
 */
void assembler :: assembleMassMatrix(const model &m)
{
    // Reset global mass matrix
    _m->setZero();
	
	// loop through bodies
    BOOST_FOREACH(modelpart* p, m.theParts )
        p->integrateMassMatrix(*this);
	
	_m->finalizeAssembly();
    DebugMessage("Consistent mass assembled.");
}





#ifdef WITHTBB

class tbbIntegrateNodalForces
{
public:
    tbbIntegrateNodalForces( vector<evalspot*>& evs) : my_ev(&evs){}
    
    void operator()( const tbb::blocked_range<size_t>& r) const
    {
        vector<evalspot*> &ev( *my_ev );
        for (size_t a=r.begin(); a!=r.end(); a++)
            ev[a]->integrateDEnergy();
    }
    
private:
    vector<evalspot*> *const my_ev;
};

#endif

/*
bool modelpart :: integrateNodalForces()
{
	if (! isActive() ) return false;
	
#ifdef WITHTBBBBBBB
	// tbb threaded version with tbb
    tbb::parallel_for(tbb::blocked_range<size_t>(0,evalspots.size()), 
                      tbbIntegrateNodalForces(evalspots),
                      tbb::auto_partitioner() );
#else
	//serial version
    
#endif
	return true;
}
*/



void assembler :: assembleNodalForces(model &m, const double time, const bool bodies, 
                                      const bool contacts, const bool external) const
{
    BOOST_FOREACH(modelpart* b, m.theParts )
        b->integrateDEnergy();
    
    
    BOOST_FOREACH(modelpart* b, m.theParts )
        b->changeSignOfNodalForces();
    

	// loop over the point loads and transfer their loads to the nodal forces
	if (external)
    {
        BOOST_FOREACH(loading* ld, m.theLoads)
            ld->transferLoadsToNodes(time);
    }
    
    DebugMessage("Nodal forces assembled.");
}




/* the nodal reactions is the sum of all the internal forces plus all the
 * distributed body forces. The internal forces include the inertial ones,
 * in a dynamic simulation.
 */
void assembler :: assembleReactions(model &m, const double time) const
{
	m.resetNodalForces();
	
	// the reactions are just the sum of all the computable forces: inertial + internal - fexternal
	// the classical reactions are just the ones at the nodes with no dofs
	assembleNodalForces(m, time, true, true, true);
    
    BOOST_FOREACH(modelpart* b, m.theParts )
        b->changeSignOfNodalForces();

	DebugMessage("Reactions assembled.");
}




/* Loops over all the elements obtaining the residual contribution and adding it
 to the global residual vector
 It can assemble, according to argument 'part':
 part = 'S'  : only the static part of the residual
 part = 'D'  : only the dynamic part of the residual
 part = 'A'  : the whole residual
 
 *Warning: the eqtime is not the time for which we are solving (normally referred to tn+1) but
 the time at which the equilibrium is enforced, which might be between tn and tn+1
 
 *The imposed boundary conditions are enforced not by incrementing the corresponding dofs,
 but by setting the corresponding increments, and modifying the residual.
 
 */
void assembler :: assembleResidual(const model &m, double eqtime, char part)
{	
	// Define assembler tasks
	_tasks.reset();
	
	char lpart = tolower(part);
	if (lpart == 's' || lpart == 'a') _tasks.static_residual  = true;
	else                              _tasks.static_residual  = false;
	if (lpart == 'd' || lpart == 'a') _tasks.dynamic_residual = true;
	else                              _tasks.dynamic_residual = false;
	
	extern analysisTypeT global_analysistype;
	if (global_analysistype != FEANALYSIS_TRANSIENT && global_analysistype != FEANALYSIS_STAGGERED)
		_tasks.dynamic_residual = false;
	
	
	// Reset global residual
	_v->setZero();
	
    BOOST_FOREACH(modelpart* bd, m.theParts )                 
        assembleGenericContribution(*bd);
    	
	assembleExternalForces(m.theLoads, eqtime, 1.0);
	DebugMessage("Global residual and tangent assembled.");
}




void assembler :: assembleResidualContribution(modelpart& theBody)
{
	if ( !theBody.isActive() ) return;
	
    
    std::vector<element*>::iterator iter = theBody.elements.begin();
    while (iter != theBody.elements.end())
    {
        (*iter)->contributeToResidual(*this);
        ++iter;
    }    
}




/* Loops over all the elements obtaining the residual AND tangent contribution and adding it
 to the global vector/matrix. This is done in such a way that is faster that computing
 the residual plus computing the tangent.
 
 Warning: the eqtime is not the time for which we are solving (normally referred to tn+1) but
 the time at which the equilibrium is enforced, which might be between tn and tn+1
 */
void assembler :: assembleResidualTangent(const model &m, double eqtime)
{		
	// Define assembler tasks
	_tasks.reset();
	_tasks.static_residual  = true;
	_tasks.static_tangent   = true;

	if (_time_derivatives > 0)
	{
		_tasks.dynamic_residual = true;
		_tasks.dynamic_tangent  = true;
	}
	
	// Reset global residual and tangent
	_v->setZero();
	_m->setZero();
	
	// indicate that, if the tangent has the LU factors already computed, they are obsolete
	_m->setNotFactorized();
	
    
	// loop through parts
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        if (bd->isActive())
        {
            BOOST_FOREACH(evalspot* ev, bd->evalspots )
                ev->integrateResidualTangent(*this);
            
            BOOST_FOREACH(element* el, bd->elements)
                el->contributeToResidualTangent(*this);
        }
    }


	assembleExternalForces(m.theLoads, eqtime, 1.0);
	
	DebugMessage("Global residual and tangent assembled.");
	_m->finalizeAssembly();
}



#ifdef WITHTBBBB
// data structure for tbb threaded version
class tbbResidualTangentContribution : public assembler
{	
    
private:
    assembler               *my_assembler;
    vector<evalspot*>       *my_ev_entities;
    
public:
	
	tbbResidualTangentContribution(vector<evalspot*> ev, assembler& a) :
        my_ev_entities(&ev), my_assembler(&a) {}
	tbbResidualTangentContribution(const tbbResidualTangentContribution &a) : 
        my_ev_entities(a.my_ev_entities), my_assembler(a.my_assembler) {}	
	
    void operator()( const tbb::blocked_range<size_t>& r) const
    {
        for (size_t a=r.begin(); a!=r.end(); a++)	
            (*my_ev_entities)[a]->integrateResidualTangent(*my_assembler);
    }
};
#endif



void assembler :: assembleResidualTangentContribution(modelpart& sm)
{
    if (!sm.isActive()) return;
        
#ifdef WITHTBBB
	// tbb threaded version with tbb
	// we need to put a lock in the assembly
	tbb::parallel_for(tbb::blocked_range<size_t>(0,elements.size()), slicedAssembleRT(elements, residual, tangent));
	
#else

    BOOST_FOREACH(evalspot* ev, sm.evalspots )
        ev->integrateResidualTangent(*this);

#endif
    _m->finalizeAssembly();
}




void assembler :: assembleSpecialIntegral(const model &m)
{
	// reset global residual vector
	_v->setZero();
    
	// loop through bodies
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(element* el, bd->elements)
		{
			// Go to the element function and compute integral
			el->specialIntegral(*this);
		}
	}		
    _m->finalizeAssembly();
}




void assembler :: assembleStaticTangent(const model &m)
{
	// Define assembler tasks
	_tasks.reset();
	_tasks.static_tangent = true;
	
    // Reset global tangent
    _m->setZero();
    _m->setNotFactorized();
    
    // loop through bodies
    BOOST_FOREACH(modelpart* bd, m.theParts )
		assembleGenericContribution( *bd );
		
	_m->finalizeAssembly();
    DebugMessage("Static tangent assembled.");
}





/* Loops over all the elements of the analysis model, requesting their element
 tangent matrices and adding their contribution to the global tangent. 
 */
void assembler :: assembleTangent(const model &m)
{
	// Define assembler tasks
	_tasks.reset();
	_tasks.static_tangent = true;
	extern double  ftgstiff;
	ftgstiff = 1.0;
	
	if (_time_derivatives > 0)
		_tasks.dynamic_tangent  = true;
	
    // Reset global tangent
    _m->setZero();
	
    // indicate that, if the tangent has the LU factors already computed, they are obsolete
    _m->setNotFactorized();
    
    // loop through bodies
    BOOST_FOREACH(modelpart* bd, m.theParts )
		assembleGenericContribution( *bd );
		
	_m->finalizeAssembly();
	DebugMessage("Total Tangent assembled.");
}



BCassembler :: BCassembler(longvector& v, sparseMatrix& m) :
	assembler(v, m)
{
	// Define assembler tasks
	_tasks.reset();
	_tasks.static_tangent   = true;
	if (_time_derivatives > 0)
		_tasks.dynamic_tangent  = true;
	
}


BCassembler :: BCassembler(const assembler &a)
{
	_v = a._v;
	_m = a._m;
	
	// Define assembler tasks
	_tasks.reset();
	_tasks.static_tangent   = true;
	if (_time_derivatives > 0)
		_tasks.dynamic_tangent  = true;
}



void BCassembler :: assemble(const itensor& m, vectorDofset& dsa, vectorDofset& dsb)
{
	for (size_t i=0; i<3; ++i)
	{
		if (! dsa.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( dsb.isConstrained(j) )
					_v->data[ dsa.getID(i) ] -= m(i,j) *  dsb.increment(j);
		}
	}	
}


void BCassembler :: assemble(const double f, vectorDofset& dsa, vectorDofset& dsb)
{
	for (size_t i=0; i<3; ++i)
	{
		if (! dsa.isConstrained(i) ) 
		{
				if ( dsb.isConstrained(i) )
					_v->data[ dsa.getID(i) ] -= f *  dsb.increment(i);
		}
	}	
}



void BCassembler :: assemble(const itensor& m, vectorDofset& dsa, rotationDofset& dsb)
{
	for (size_t i=0; i<3; ++i)
	{
		if (! dsa.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( dsb.isConstrained(j) )
					_v->data[ dsa.getID(i) ] -= m(i,j) *  dsb.increment(j);
		}
	}	
}

void BCassembler :: assemble(const itensor& m, rotationDofset& dsa, vectorDofset& dsb)
{
	for (size_t i=0; i<3; ++i)
	{
		if (! dsa.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( dsb.isConstrained(j) )
					_v->data[ dsa.getID(i) ] -= m(i,j) *  dsb.increment(j);
		}
	}	
}

void BCassembler :: assemble(const itensor& m, rotationDofset& dsa, rotationDofset& dsb)
{
	for (size_t i=0; i<3; ++i)
	{
		if (! dsa.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( dsb.isConstrained(j) )
					_v->data[ dsa.getID(i) ] -= m(i,j) *  dsb.increment(j);
		}
	}	
}



void BCassembler :: assemble	(const ivector& v, vectorDofset& dsa, scalarDofset& dsb)
{
	for (size_t i=0; i<3; i++)
		if (! dsa.isConstrained(i) ) 
			if ( dsb.isConstrained() ) 
				_v->data[ dsa.getID(i) ] -= v[i] * dsb.increment();
}


void BCassembler :: assemble	(const ivector& v, scalarDofset& dsa, vectorDofset& dsb)
{
	if (! dsa.isConstrained() ) 
		for (size_t j=0; j<3; j++)
			if ( dsb.isConstrained(j) ) 
				_v->data[ dsa.getID() ] -= v[j] * dsb.increment(j);
}


void BCassembler :: assemble	(const double f, scalarDofset& dsa, scalarDofset& dsb)
{
	if (! dsa.isConstrained() ) 
		if ( dsb.isConstrained() ) 
			_v->data[ dsa.getID() ] -= f * dsb.increment();
}





#ifdef WITHMPI

/************* Must be reviewed. What will happen when vector *_v changes its length? *************/
MPIImplicitAssembler :: MPIImplicitAssembler() :
assembler() 
{	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
}


MPIImplicitAssembler :: MPIImplicitAssembler(longvector& v, sparseMatrix& m) :
assembler(v, m) 
{
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE ,v.length, &Vecb);
	
	VecScatterCreateToAll(Vecb,&scatterCtx,&VecbCopy);
	
	double *array;
	int localFirstCol, localLastCol;
	
	VecGetOwnershipRange(Vecb, &localFirstCol, &localLastCol);
	
	VecGetArray(Vecb, &array);
	Dcopy(v.data + localFirstCol, array, localLastCol-localFirstCol);	
	VecRestoreArray(Vecb, &array);
	
	VecAssemblyBegin(Vecb);
	VecAssemblyEnd(Vecb);
	
	VecScatterBegin(scatterCtx, Vecb, VecbCopy,INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd  (scatterCtx, Vecb, VecbCopy,INSERT_VALUES, SCATTER_FORWARD);
}


MPIImplicitAssembler :: ~MPIImplicitAssembler()
{
	VecDestroy(Vecb);
	VecDestroy(VecbCopy);	
	VecScatterDestroy(scatterCtx);
}


void MPIImplicitAssembler :: assembleMassMatrix(const model &m)
{

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	
    // Reset global mass matrix
    _m->setZero();
	
	if(rank == 0){
		
		// Define assembler tasks
		_tasks.reset();
		_tasks.mass_matrix = true;
		
		
		// loop through bodies. Contactpairs are assumed to provide no mass at all
		int nbodies( m.getNBodies() );
		for (size_t b=0; b<nbodies; b++)
		{
			body& aBody( m.getBody(b) );
			
			// loop through elements
			elementIterator iter= aBody.elements.begin();
			
			while ( iter != aBody.elements.end() )
			{
				(*iter)->integrateMassMatrix(*this);			
				++iter;
			}
		}
	}
	
	_m->finalizeAssembly();
    DebugMessage("Consistent mass assembled.");

}




void MPIImplicitAssembler :: assembleResidual(const model &m, double eqtime, char part)
{

	if(rank == 0){
		// Define assembler tasks
		_tasks.reset();
		
		char lpart = tolower(part);
		if (lpart == 's' || lpart == 'a') _tasks.static_residual  = true;
		else                              _tasks.static_residual  = false;
		if (lpart == 'd' || lpart == 'a') _tasks.dynamic_residual = true;
		else                              _tasks.dynamic_residual = false;
		
		extern analysisTypeT global_analysistype;
		if (global_analysistype != FEANALYSIS_TRANSIENT && global_analysistype != FEANALYSIS_STAGGERED)
			_tasks.dynamic_residual = false;
		
		
		
		// Reset global residual
		_v->setZero();
		
		// loop through bodies
		int nbodies( m.theBodies.size() );
		for (size_t b=0; b<nbodies; b++)
			assembleResidualContribution( m.getBody(b) );
		
		// loop through contactpairs
		for (size_t b=0; b<m.theContactpairs.size() ; b++)
			if (m.getContactpair(b).isActive()) 
				assembleResidualContribution( m.getContactpair(b) );
		
		assembleExternalForces(m, eqtime, 1.0);
	}
	
	int *indecesRows = new int[_v->length];
	for (size_t i=0; i< _v->length; i++) indecesRows[i] = i;
	
	if(rank == 0){
		VecSetValues(Vecb, _v->length, indecesRows, _v->data, INSERT_VALUES);
	}
	VecAssemblyBegin(Vecb);
	VecAssemblyEnd(Vecb);
	double *array;
	
	
	VecScatterBegin(scatterCtx, Vecb, VecbCopy,INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd  (scatterCtx, Vecb, VecbCopy,INSERT_VALUES, SCATTER_FORWARD);
	
	if(rank!=0){
		VecGetArray(VecbCopy, &array);
		Dcopy(array, _v->data, _v->length);	
		VecRestoreArray(VecbCopy, &array);
	}
	
	DebugMessage("Global residual and tangent assembled.");
	_m->finalizeAssembly();
}



void MPIImplicitAssembler :: assembleResidualTangent(const model &m, double eqtime)
{			

		// Define assembler tasks
		_tasks.reset();
		_tasks.static_residual  = true;
		_tasks.static_tangent   = true;
		
		if (_time_derivatives > 0)
		{
			_tasks.dynamic_residual = true;
			_tasks.dynamic_tangent  = true;
		}
		
		
		// Reset global residual and tangent
		_v->setZero();
		_m->setZero();
	
	if(rank == 0){
		// indicate that, if the tangent has the LU factors already computed, they are obsolete
		_m->setNotFactorized();
		
		// loop through bodies
		int nbodies( m.getNBodies() );
		for (size_t b=0; b<nbodies; b++)
			assembleResidualTangentContribution( m.getBody(b) );
			
		// loop through contactpairs
		int ncp( m.theContactpairs.size() );
		for (size_t b=0; b<ncp; b++)
			if (m.getContactpair(b).isActive()) 
				assembleResidualTangentContribution( m.getContactpair(b) );
	
				assembleExternalForces(m, eqtime, 1.0);
	}
	
	DebugMessage("Global residual and tangent assembled.");
	
	int *indecesRows = new int[_v->length];
	for (size_t i=0; i< _v->length; i++) indecesRows[i] = i;
	
	if(rank == 0){
		VecSetValues(Vecb, _v->length, indecesRows, _v->data, INSERT_VALUES);
	}
	VecAssemblyBegin(Vecb);
	VecAssemblyEnd(Vecb);
	
	double *array;
		
	VecScatterBegin(scatterCtx, Vecb, VecbCopy,INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd  (scatterCtx, Vecb, VecbCopy,INSERT_VALUES, SCATTER_FORWARD);
	
	if(rank!=0){
		VecGetArray(VecbCopy, &array);
		Dcopy(array, _v->data, _v->length);	
		VecRestoreArray(VecbCopy, &array);
	}
	
	_m->finalizeAssembly();


}






void MPIExplicitAssembler :: initialize() 
{
	list< pair<node*, vector<int> > >::iterator itcon = linkedModel->boundaryConnectivity.begin() ;
	nconnections = 0 ;
	while(itcon!=linkedModel->boundaryConnectivity.end()) {
		nconnections += (*itcon).second.size();
	    ++itcon ;
	}
	
	transferTable = new mpinodebuffer [nconnections] ;
 	requestTableRec = new MPI_Request [nconnections] ;
	requestTableSend = new MPI_Request [nconnections] ;
    itcon = linkedModel->boundaryConnectivity.begin() ;
	int iconnection = 0 ;
	while(itcon!=linkedModel->boundaryConnectivity.end()) {
		int ipart = 0 ;
		while(ipart<(*itcon).second.size()) 
		{
			transferTable[iconnection].theNode   = (*itcon).first;
			transferTable[iconnection].partition = (*itcon).second[ipart] ;	
			transferTable[iconnection].label  = (*itcon).first->getLabel() ;	
			
			++ipart ; 
			++iconnection ;
		}
		++itcon ;	
	}
	
	//  	requestBoundaryMasses();
	//	sendBoundaryMasses();
	//	requestTableRec->Waitall(nconnections, requestTableRec);
	//	assembleMasses();
	
	
	// inicializacion para comunicacion sin bloqueos
	// primero num_proc 
	mpirank = MPI::COMM_WORLD.Get_rank();
	mpisize = MPI::COMM_WORLD.Get_size();
	for (size_t i=0; i < mpisize-1 ; i++)
	{
		int j ; 
		j = mpirank + i + 1 ; 
		if (j>mpisize-1) j = j -mpisize  ; 
		num_proc.push_back(j) ;
	}
	// ahora voy a meter los punteros a los nodos de cada frontera en nodos_part
	// voy a recorrer p-1 veces boundaryConnectivity (no es muy eficiente
	// pero solo se hace una vez
	vector<node *> tmpnodosboun ; 
	bsend = new double * [mpisize-1] ; 
	brec = new double * [mpisize-1] ; 
	for (size_t i=0;i<mpisize-1;i++) {
		tmpnodosboun.clear();
		itcon = linkedModel->boundaryConnectivity.begin() ;
		while(itcon!=linkedModel->boundaryConnectivity.end()) {
			int ipart = 0 ;
			while(ipart<(*itcon).second.size()) 
			{
				if (((*itcon).second[ipart])==num_proc[i] )
					tmpnodosboun.push_back((*itcon).first) ;	
				++ipart ; 
			}
			++itcon ;	
		}
		nodos_part.push_back(tmpnodosboun) ; 
		num_nodos_part.push_back(tmpnodosboun.size());
		bsend[i] = new double  [tmpnodosboun.size()*3] ; 
		brec[i] = new double   [tmpnodosboun.size()*3] ; 
	}
	assembleMasses2();
	
}


MPIExplicitAssembler :: MPIExplicitAssembler(longvector& v, sparseMatrix& m) : 
explicitAssembler(v, m), initmesh(0)
{
	mpirank = MPI::COMM_WORLD.Get_rank();
}


MPIExplicitAssembler :: MPIExplicitAssembler() :
explicitAssembler(), initmesh(0)
{
	mpirank = MPI::COMM_WORLD.Get_rank();
} 


void MPIExplicitAssembler :: initMesh(model &m )
{
	if (initmesh ==0) 
	{
		linkedModel = &m ;  
		// de momento no puedo poner la malla en linkedModel (debe ser porque es una referencia 
		list< pair<node*, vector<int> > >::iterator itcon = m.boundaryConnectivity.begin() ;
		nconnections = 0 ;
        
		//   logger::mainlog << "\n boundaryConnect .size().. = " << m.boundaryConnectivity.size() ;
    	while(itcon!=m.boundaryConnectivity.end()) {
			nconnections += (*itcon).second.size();
			++itcon ;
		}
		
		transferTable = new mpinodebuffer [nconnections] ;
		requestTableRec = new MPI_Request [nconnections] ;
		requestTableSend = new MPI_Request [nconnections] ;
		itcon = m.boundaryConnectivity.begin() ;
		int iconnection = 0 ;
		while(itcon!=m.boundaryConnectivity.end()) {
			int ipart = 0 ;
			while(ipart<(*itcon).second.size()) 
			{
				transferTable[iconnection].theNode   = (*itcon).first;
				transferTable[iconnection].partition = (*itcon).second[ipart] ;	
				transferTable[iconnection].label  = (*itcon).first->getLabel() ;	
				
				++ipart ; 
				++iconnection ;
			}
			++itcon ;
			
		}
		//   logger::mainlog << "\n initMesh nconnections " <<  nconnections ; 
		
		initmesh = 1 ; 
	}
	
}

MPIExplicitAssembler :: ~MPIExplicitAssembler()
{
	delete [] requestTableSend ;
	delete [] requestTableRec;
	delete [] transferTable;
}



void MPIExplicitAssembler :: assembleMasses()
{	
	for (size_t i=0; i<nconnections; i++)
	{
		node& nd = linkedModel->getNode(transferTable[i].receive.label);
		nd.incrementMass( transferTable[i].receive.f[0] );
	}	
}


void MPIExplicitAssembler :: requestBoundaryMasses()
{
	int tama ;
	tama = sizeof(iobuffer);
	for (size_t i=0; i<nconnections; i++)
		requestTableRec[i] = MPI::COMM_WORLD.Irecv( (void *) &(transferTable[i].receive), tama, MPI::BYTE, transferTable[i].partition, 2) ;
}



void MPIExplicitAssembler :: sendBoundaryMasses()
{
	int tama ;
	tama = sizeof(iobuffer);
	for (size_t i=0; i<nconnections; i++)
	{
		transferTable[i].send.label = transferTable[i].label;
		transferTable[i].send.f[0] = transferTable[i].theNode->getMass();
		requestTableSend[i] = MPI::COMM_WORLD.Isend( (void *) &(transferTable[i].send), tama, MPI::BYTE, transferTable[i].partition, 2);
	}
}




void MPIExplicitAssembler :: setLinkedMesh(model& m)
{
	linkedModel = &m;
	initialize();
}


void MPIExplicitAssembler :: assembleNodalForces(model &m, const double time, const bool bodies, const bool contacts, const bool external) const
{
	assembler::assembleNodalForces(m, time, bodies, contacts, external);
	
	// aqui va a ir el control de el proceso y que envios hay que hacer
	// voy a probar el caso de solo dos procesos
	// esto es mas sofisticado hay que recorrer el numero de fronteras ...
	// es el tamanho del salto 
	// salto=ipart+ 1
	// prueba para mpisize potencia de 2  
	int salto ;   
	int paridad ; 
	
    for (size_t ipart = 0 ; ipart<mpisize/2 ; ipart++) {
		salto = mpirank/2 ;  
		if (ipart==0) paridad = mpirank % 2 ;  
		if (ipart != 0 ) paridad = salto % (ipart) ;   
		
		// prueba
		paridad = mpirank / (ipart+1) ; paridad = paridad % 2 ; 
		
		if (paridad==0) 
		{
			if(ipart+1==mpisize/2) 
			{
				sendBoundaryForces(ipart) ; 
				recvBoundaryForces(ipart);
			} 
			else 
			{
				sendBoundaryForces(ipart) ; 
				sendBoundaryForces(mpisize-ipart-2) ; 
				recvBoundaryForces(ipart);
				recvBoundaryForces(mpisize-ipart-2);
			}
		}
		else 
		{
			if(ipart+1==mpisize/2) {
				recvBoundaryForces(ipart);
				sendBoundaryForces(ipart) ; 
            } else {
				recvBoundaryForces(mpisize-ipart-2);
				recvBoundaryForces(ipart);
				sendBoundaryForces(mpisize-ipart-2) ; 
				sendBoundaryForces(ipart) ; 
			}
		}
    }
	
 	
	for (size_t ipart = 0 ; ipart<mpisize/2 ; ipart++) {
		if(ipart+1==mpisize/2) {
			assembleBoundaryForces(ipart);
		}
		else {
			assembleBoundaryForces(ipart);
			assembleBoundaryForces(mpisize-ipart-2);
		}
    }
	
	
}


void MPIExplicitAssembler :: requestBoundaryForces()
{
	int tama ;
	tama = sizeof(iobuffer);
	//        logger::mainlog << "\n nconnections = " << nconnections ;
	for (size_t i=0; i<nconnections; i++)
		requestTableRec[i] = MPI::COMM_WORLD.Irecv( (void *) &(transferTable[i].receive) , tama, MPI::BYTE, transferTable[i].partition,1) ;
}



void MPIExplicitAssembler :: sendBoundaryForces()
{
	int tama ;
	tama = sizeof(iobuffer);
 	for (size_t i=0; i<nconnections; i++)
	{
		Unode* unode = dynamic_cast<Unode*>(transferTable[i].theNode);
		Dcopy( unode->force().components(), transferTable[i].send.f, 3);
		transferTable[i].send.label  =  unode->label ; 
		requestTableSend[i] = MPI::COMM_WORLD.Isend( (void *) &(transferTable[i].send), tama, MPI::BYTE, transferTable[i].partition,1) ;
	}
	
}


void MPIExplicitAssembler :: recvBoundaryForces(int ipart) const
{
	int tama ;
	tama = num_nodos_part[ipart]* 3 ;  
	if (tama!=0) MPI::COMM_WORLD.Recv( (void *) (brec[ipart]), tama, MPI::DOUBLE, num_proc[ipart],1) ;
}


void MPIExplicitAssembler :: sendBoundaryForces(int ipart) const
{
	int tama ;
	tama = num_nodos_part[ipart]* 3 ;  
	
	
	for (size_t i=0; i< num_nodos_part[ipart]; i++)
	{
		Unode* unode = dynamic_cast<Unode*>(nodos_part[ipart][i]);
		Dcopy( unode->force().components(), bsend[ipart]+i*3, 3);
	}
	
 	if (tama!=0) MPI::COMM_WORLD.Send( (void *) (bsend[ipart]), tama, MPI::DOUBLE, num_proc[ipart],1) ;
	
}


void MPIExplicitAssembler :: assembleBoundaryForces(int ipart) const
{
	
	ivector vtmp ; 
	node *pnode ;
 	for (size_t i=0; i< num_nodos_part[ipart]; i++)
	{
		//		unode = dynamic_cast<Unode*>(nodos_part[ipart][i]);
		pnode = nodos_part[ipart][i];
		Dcopy( brec[ipart]+i*3, vtmp.components(), 3);
		pnode->assembleForce( vtmp);
	}
}

void MPIExplicitAssembler :: assembleBoundaryForces() const
{
	int tama ;
	tama = sizeof(iobuffer);
	ivector vtmp;
	
 	for (size_t i=0; i<nconnections; i++)
	{
		Dcopy( transferTable[i].receive.f, vtmp.components(), 3);
		node &pnode = linkedModel->getNode(transferTable[i].receive.label) ; 
		pnode.assembleForce( vtmp);
	}	
}


void MPIExplicitAssembler :: recvBoundaryMasses(int ipart)
{
	// utilizo las mismas estructuras que para fuerzas pero con
	// menos tamanho. Por lo tanto valen. 
	int tama ;
	tama = num_nodos_part[ipart]* 1 ;  
	if (tama!=0) MPI::COMM_WORLD.Recv( (void *) (brec[ipart]), tama, MPI::DOUBLE, num_proc[ipart],2) ;
}


void MPIExplicitAssembler :: sendBoundaryMasses(int ipart)
{
	int tama ;
	tama = num_nodos_part[ipart]* 1 ;
	
	for (size_t i=0; i< num_nodos_part[ipart]; i++)
	{
		Unode* unode = dynamic_cast<Unode*>(nodos_part[ipart][i]);
		bsend[ipart][i] = unode->getMass() ; 
	}
	
 	if (tama!=0) MPI::COMM_WORLD.Send( (void *) (bsend[ipart]), tama, MPI::DOUBLE, num_proc[ipart],2) ;
	
}


void MPIExplicitAssembler :: assembleBoundaryMasses(int ipart)
{
	
	ivector vtmp ; 
	node *pnode ;
 	for (size_t i=0; i< num_nodos_part[ipart]; i++)
	{
		pnode = nodos_part[ipart][i];
		pnode->incrementMass( brec[ipart][i]);
	}
}


void MPIExplicitAssembler :: assembleMasses2()
{
	// con bloqueo 
	int salto ;   
	int paridad ; 
	
    for (size_t ipart = 0 ; ipart<mpisize/2 ; ipart++) {
		salto = mpirank/2 ;  
		
		// prueba
		paridad = mpirank / (ipart+1) ; paridad = paridad % 2 ; 
		if (paridad==0) 
		{
			if(ipart+1==mpisize/2) {
				sendBoundaryMasses(ipart) ; 
				recvBoundaryMasses(ipart);
			} else {
				sendBoundaryMasses(ipart) ; 
				sendBoundaryMasses(mpisize-ipart-2) ; 
				recvBoundaryMasses(ipart);
				recvBoundaryMasses(mpisize-ipart-2);
			}
		}
		else 
		{
			if(ipart+1==mpisize/2) {
				recvBoundaryMasses(ipart);
				sendBoundaryMasses(ipart) ; 
            } else {
				recvBoundaryMasses(mpisize-ipart-2);
				recvBoundaryMasses(ipart);
				sendBoundaryMasses(mpisize-ipart-2) ; 
				sendBoundaryMasses(ipart) ; 
			}
		}
    }
	
 	
	for (size_t ipart = 0 ; ipart<mpisize/2 ; ipart++) {
		if(ipart+1==mpisize/2) {
			assembleBoundaryMasses(ipart);
		}
		else {
			assembleBoundaryMasses(ipart);
			assembleBoundaryMasses(mpisize-ipart-2);
		}
    }
}

#endif


radicalLumperAssembler :: radicalLumperAssembler(longvector& lv)
{
	_v = &lv;
}



void radicalLumperAssembler :: assemble	(const itensor& t, vectorDofset& uda, vectorDofset& udb)
{	
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( ! udb.isConstrained(j) ) 
					_v->data [ uda.getID(i) ] += fabs( t(i,j) );
		}
	}
}


void radicalLumperAssembler :: assembleStaticTangent(const model &m)
{
	// Define assembler tasks
	_tasks.reset();
	_tasks.static_tangent = true;
	
    // Reset global lumped tangent, which is just a vector
	_v->setZero();
	
    BOOST_FOREACH(modelpart* bd, m.theParts )
        assembleGenericContribution(*bd);
	
    DebugMessage("Static tangent assembled.");
}



singleElementAssembler :: singleElementAssembler(const element& e) :
    eres( e.getNDofs() ),
    emat( e.getNDofs(), e.getNDofs() )
{
	int local=0;
	for (size_t a=0; a< e.getNNodes(); a++)
	{
		node& nn = e.getNode(a);
		for (size_t b=0; b< nn.getNDofs(); b++)
		{
			if ( nn.getID(b) >= 0) globalToLocal[nn.getID(b)] = local;
			local++;
		}
	}
}


singleElementAssembler :: ~singleElementAssembler()
{
	
}


void singleElementAssembler :: assemble	(const ivector& f, const vectorDofset& uds)
{
	if ( ! uds.isConstrained(0) ) eres[ globalToLocal[uds.getID(0)] ] += f[0];
	if ( ! uds.isConstrained(1) ) eres[ globalToLocal[uds.getID(1)] ] += f[1];
	if ( ! uds.isConstrained(2) ) eres[ globalToLocal[uds.getID(2)] ] += f[2];
}

void singleElementAssembler :: assemble	(const ivector& f, const directorDofset& uds)
{
	if ( ! uds.isConstrained(0) ) eres[ globalToLocal[uds.getID(0)] ] += f[0];
	if ( ! uds.isConstrained(1) ) eres[ globalToLocal[uds.getID(1)] ] += f[1];
}


void singleElementAssembler :: assemble	(const itensor& t, vectorDofset& uda, vectorDofset& udb)
{
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( ! udb.isConstrained(j) ) 
					emat.addTo( t(i,j), globalToLocal[uda.getID(i)], globalToLocal[udb.getID(j)] );
		}
	}
}



void singleElementAssembler :: assemble	(const double t, vectorDofset& uda, vectorDofset& udb)
{
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) && ! udb.isConstrained(i) ) 
			emat.addTo( t, globalToLocal[uda.getID(i)], globalToLocal[udb.getID(i)] );
	}
}



void singleElementAssembler :: assemble	(const itensor& t, directorDofset& uda, vectorDofset& udb)
{
	for (size_t i=0; i<2; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<3; j++)
				if ( ! udb.isConstrained(j) ) 
					emat.addTo( t(i,j), globalToLocal[uda.getID(i)], globalToLocal[udb.getID(j)] );
		}
	}
}


void singleElementAssembler :: assemble	(const itensor& t, vectorDofset& uda, directorDofset& udb)
{
	for (size_t i=0; i<3; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<2; j++)
				if ( ! udb.isConstrained(j) ) 
					emat.addTo( t(i,j), globalToLocal[uda.getID(i)], globalToLocal[udb.getID(j)] );
		}
	}
}


void singleElementAssembler :: assemble	(const itensor& t, directorDofset& uda, directorDofset& udb)
{
	for (size_t i=0; i<2; i++)
	{
		if (! uda.isConstrained(i) ) 
		{
			for (size_t j=0; j<2; j++)
				if ( ! udb.isConstrained(j) ) 
					emat.addTo( t(i,j), globalToLocal[uda.getID(i)], globalToLocal[udb.getID(j)] );
		}
	}
}



void singleElementAssembler :: assemble	(const ivector& v, vectorDofset& dsa, scalarDofset& dsb)
{
	for (size_t i=0; i<3; i++)
		if (! dsa.isConstrained(i) ) 
			if ( ! dsb.isConstrained() ) 
				emat.addTo( v[i], globalToLocal[dsa.getID(i)], globalToLocal[dsb.getID()] );
}


void singleElementAssembler :: assemble	(const ivector& v, scalarDofset& dsa, vectorDofset& dsb)
{
	if (! dsa.isConstrained() ) 
		for (size_t j=0; j<3; j++)
			if ( ! dsb.isConstrained(j) ) 
				emat.addTo( v[j], globalToLocal[dsa.getID(0)], globalToLocal[dsb.getID(j)] );
}


void singleElementAssembler :: assemble	(const double d, scalarDofset& dsa, scalarDofset& dsb)
{
	if (! dsa.isConstrained() ) 
		if ( ! dsb.isConstrained() ) 
			emat.addTo( d , globalToLocal[dsa.getID()], globalToLocal[dsb.getID()] );
}



void singleElementAssembler :: assemble	(const double f, const scalarDofset& sds)
{
	if (! sds.isConstrained() ) 
		eres.data[ globalToLocal[sds.getID(0)] ] += f;
}



void singleElementAssembler :: assemble(const longvector& v, const vector<int>& doflabels)
{
    for (size_t i=0; i<doflabels.size(); i++)
    {
        if (doflabels[i] >= 0) eres.data[ globalToLocal[doflabels[i]] ] += v[i];
    }
}



void singleElementAssembler :: assemble(const matrix& m, const vector<int>& doflabels)
{
    for (size_t i=0; i<doflabels.size(); i++)
    {
        if (doflabels[i] >= 0)
        {
            for (size_t j=0; j<doflabels.size(); j++)
                if (doflabels[j] >= 0) emat(globalToLocal[doflabels[i]], globalToLocal[doflabels[j]]) += m(i,j);
        }
    }    
}






