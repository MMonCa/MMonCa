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
 * assembler.h
 *
 * i. romero
 * july 2000
 *
 * functions that assemble the tangent matrix and right hand side vector
 */


#ifndef _assembler_h
#define _assembler_h

#include "Elements/elmttasks.h"
#include "Math/matrix.h"
#include "Math/vector.h"
#include "Math/tensor.h"

#include <vector>
#include <map>
#include <list>
#ifdef WITHMPI
#include "mpi.h"
#endif


class body;
class contactpair;
class element;
class loading;
class matrix;
class model;
class node;
class scalarDofset;
class sparseMatrix;
class modelpart;
class vectorDofset;
class Unode;

class directorDofset;
class scalarDofset;
class rotationDofset;
class vectorDofset;


class assembler
{
    
public:	
                    assembler();
                    assembler(longvector& v, sparseMatrix& m);
	virtual         ~assembler(){}

	void			setNumberOfTimeDerivatives(const int n){ _time_derivatives=n;};
	longvector&		getVector(){return *_v;}
	
	///MMM: function to reassign the _m matrix, keeping the rest of the assembler.
    ///III: we need to talk about this funtion in september
    inline void     reassignTangent(sparseMatrix& M ){_m = &M;}

	// for explicit methods
	virtual void	assembleNodalForces(model &m, const double time, const bool bodies=true, 
                                        const bool contacts=true, const bool external=true) const;
	
	
	// functions to go through elements assembling global matrices, vectors
	void			assembleGenericContribution(modelpart& theBody);
	void			assembleDampingMatrix  (model &m, sparseMatrix &damping);
	void			assembleExternalForces (const std::list<loading*>& loads, double time, double factor);
	void			assembleL2ProjectionMatrix(const model &m);
	virtual void	assembleMassMatrix     (const model &m);
	void			assembleReactions      (model &m, const double time) const;
	virtual void	assembleResidual       (const model &m, double time, char part='A');
	virtual void	assembleResidualTangent(const model &m, double time);
	
	virtual void	assembleSpecialIntegral(const model &m);
	virtual void	assembleStaticTangent  (const model &m);
	virtual void	assembleTangent        (const model &m);
	
	// residual assembler
	virtual void	assemble	(const blue::ivector& f, const directorDofset& ds);
	virtual void	assemble	(const blue::ivector& f, const rotationDofset& ds);
	virtual void	assemble	(const double   d, const scalarDofset& ds);
	virtual void	assemble	(const blue::ivector& f, const vectorDofset& ds);

	
	
	// tangent assembler
	virtual void	assemble	(const blue::itensor& m, Unode& na, Unode& nb);
	virtual void	assemble	(const blue::itensor& m, vectorDofset& dsa, vectorDofset& dsb);
	virtual void	assemble	(const double   f, vectorDofset& dsa, vectorDofset& dsb); // assemble scaled I
	virtual void	assemble	(const blue::itensor& m, vectorDofset& dsa, rotationDofset& dsb);
	virtual void	assemble	(const blue::itensor& m, rotationDofset& dsa, vectorDofset& dsb);
	virtual void	assemble	(const blue::itensor& m, rotationDofset& dsa, rotationDofset& dsb);
	virtual void	assemble	(const blue::itensor& m, vectorDofset& dsa,	directorDofset& dsb);
	virtual void	assemble	(const blue::itensor& m, directorDofset& dsa, vectorDofset& dsb);
	virtual void	assemble	(const blue::itensor& m, directorDofset& dsa, directorDofset& dsb);
	
	
	virtual void	assemble	(const blue::ivector& v, vectorDofset& dsa, scalarDofset& dsb);
	virtual void	assemble	(const blue::ivector& v, scalarDofset& dsa, vectorDofset& dsb);
	virtual void	assemble	(const double   f, scalarDofset& dsa, scalarDofset& dsb);

	
	// old fashion assembler
	void			assemble	(const node& nd, const int dof, const double f);
	virtual void    assemble	(const longvector& v, const std::vector<int>& doflabels);
	virtual void	assemble	(const matrix& m,     const std::vector<int>& doflabels);

	// for staggered methods
	inline bool		belongsToOperatorSplit1() const		{ return _first_operator_split;}
	inline bool		belongsToOperatorSplit2() const		{ return _second_operator_split;}
	
	
	inline  const	elmtTasks&		getTasks() const	{return _tasks;}
	inline			elmtTasks&		getTasks()			{return _tasks;}
	
	elmtTasks       _tasks;

	
protected:
	longvector*		_v;
	sparseMatrix*	_m;
	bool			_first_operator_split;
	bool			_second_operator_split;
	int				_time_derivatives;
	
	
	friend class BCassembler;
	
	void assembleResidualTangentContribution(modelpart& sm);
	void assembleResidualContribution(modelpart& sm);

private:
	assembler(const assembler& a);
};







class BCassembler : public assembler
	{

	public:
		BCassembler(longvector& v, sparseMatrix& m);
		BCassembler(const assembler &a);

		// tangent assembler
		virtual void	assemble(const blue::itensor& m, vectorDofset& dsa, vectorDofset& dsb);
		virtual void	assemble(const double   f, vectorDofset& dsa, vectorDofset& dsb); // assemble scaled I
		virtual void	assemble(const blue::itensor& m, vectorDofset& dsa, rotationDofset& dsb);
		virtual void	assemble(const blue::itensor& m, rotationDofset& dsa, vectorDofset& dsb);
		virtual void	assemble(const blue::itensor& m, rotationDofset& dsa, rotationDofset& dsb);
		
		virtual void	assemble(const blue::ivector& v, vectorDofset& dsa, scalarDofset& dsb);
		virtual void	assemble(const blue::ivector& v, scalarDofset& dsa, vectorDofset& dsb);
		virtual void	assemble(const double   f, scalarDofset& dsa, scalarDofset& dsb);
		
		// residual assembler
		virtual void	assemble	(const blue::ivector& f, const vectorDofset& ds){assembler::assemble(f, ds);}
		virtual void	assemble	(const double   d, const scalarDofset& ds){assembler::assemble(d, ds);}
		
		
	private:
		
	};






#ifdef WITHPETSC
#include "petscda.h"
#include "petscsnes.h"

class MPIImplicitAssembler : public assembler
{
public:
	MPIImplicitAssembler();
	MPIImplicitAssembler(longvector& v, sparseMatrix& m);
	~MPIImplicitAssembler();
	
	
	void			assembleResidual       (const model &m, double time, char part='A');
	void			assembleResidualTangent(const model &m, double time);
	void			assembleMassMatrix     (const model &m);
	
private:
	int			rank;
	Vec			Vecb;
	Vec			VecbCopy;
	VecScatter	scatterCtx;
};

#endif


class explicitAssembler : public assembler
	{
	public:
		explicitAssembler() : assembler() {}
		explicitAssembler(longvector& v, sparseMatrix& m) : assembler(v, m) {}

		
	};



class radicalLumperAssembler : public assembler
{
public:
	radicalLumperAssembler(longvector& lv);
	void assembleStaticTangent  (const model &m);
	void assemble(const blue::itensor& t, vectorDofset& uda, vectorDofset& udb);

};


class UP0assembler : public assembler
	{
	public:
		UP0assembler(){ _first_operator_split = true; }
		UP0assembler(longvector& r, sparseMatrix& m) : assembler(r,m){ _first_operator_split = true; }
	};




class UP1assembler : public assembler
	{
	public:		
		UP1assembler(){ _second_operator_split = true;}
		UP1assembler(longvector& r, sparseMatrix& m) : assembler(r,m){ _second_operator_split = true;}
	};




#ifdef WITHMPI
//#include "mpi.h"

class iobuffer  
	{
	public:
		double f[3] ;
		int    label ;
	};

class mpinodebuffer 
	{
	public:
		
		iobuffer send , receive ;
		int   partition;
		node* theNode;
		int label ; 
	} ;


class MPIExplicitAssembler : public explicitAssembler
	{
		
	public:
		MPIExplicitAssembler(const commandLine &cl) {};
		MPIExplicitAssembler(longvector& v, sparseMatrix& m);
		MPIExplicitAssembler();
		~MPIExplicitAssembler();
		
		void assembleNodalForces(model &m, const double time, const bool bodies, const bool contacts, const bool external) const;
		void initialize();
		void setLinkedMesh(model& m);
		
		
	private:
		void assembleBoundaryForces() const;
		void requestBoundaryForces();
		void sendBoundaryForces();
		void recvBoundaryForces(int) const;
		void sendBoundaryForces(int) const;
		void assembleBoundaryForces(int) const;
		void recvBoundaryMasses(int);
		void sendBoundaryMasses(int);
		void assembleBoundaryMasses(int);
		
		void requestBoundaryMasses();
		void sendBoundaryMasses();
		void assembleMasses(); 
		void assembleMasses2();   // para con bloqueo
		
		void initMesh(model &) ; // es una chapucilla 
		
		
		
		int initmesh ; // parche para inicializar linkedModel 
		// si 0 se inicializa , es una chapuza
		model		*linkedModel;
		mpinodebuffer *transferTable ;
		MPI_Request *requestTableSend ;
		MPI_Request *requestTableRec ;
		int nconnections ;
		// nuevas estructuras para comunicacion con bloqueo
		vector<int> num_nodos_part ; // numero de nodos en cada 
		// frontera
		// nodos_part[0] es el numero de nodos que tiene la frontera con
		// el proceso siguiente
		// nodos_part[p-2] es el numero de nodos con el proceso anterior
		// es un vector de tamanho p-1 
		vector<int>  num_proc ; // numero de proceso que 
		// se corresponde con el orden de frontera
		//       num_proc[0] es el proceso siguiente
		vector<vector<node*> > nodos_part ; 
		// cuales son los nodos en cada frontera 
		// num_part dice cuantos y nodos_part cuales 
		// num_proc nos dice en que proceso esta esa frontera
		// buffers para envio y para recepcion
		double **bsend ;   // bsend[num_proc][nodos_part][3]
		double **brec ;    // igual 
		// nodos de las fronteras
		// no hace falta esta en nodos_part   node *pnode ;  
		int mpirank ; 
		int mpisize ; 
	};

#endif




// this assembler gathers the data of a single element, including its
// residual and its tangent 
class singleElementAssembler : public assembler
{
public:
	singleElementAssembler(const element& e);
	~singleElementAssembler();
	
	// residual assembler
	virtual void	assemble	(const blue::ivector& f, const vectorDofset& ds);
	virtual void	assemble	(const double   d, const scalarDofset& ds);
	virtual void	assemble	(const blue::ivector& f, const directorDofset& ds);
	
	virtual void	assemble	(const blue::itensor& m, vectorDofset& dsa, vectorDofset& dsb);
    virtual void	assemble	(const double   f, vectorDofset& dsa, vectorDofset& dsb); // assemble scaled I
	virtual void	assemble	(const blue::itensor& m, directorDofset& dsa, vectorDofset& dsb);
	virtual void	assemble	(const blue::itensor& m, vectorDofset& dsa, directorDofset& dsb);
	virtual void	assemble	(const blue::itensor& m, directorDofset& dsa, directorDofset& dsb);
	virtual void	assemble	(const blue::ivector& v, vectorDofset& dsa, scalarDofset& dsb);
	virtual void	assemble	(const blue::ivector& v, scalarDofset& dsa, vectorDofset& dsb);
	virtual void	assemble	(const double   f, scalarDofset& dsa, scalarDofset& dsb);

	// old fashion
	virtual void    assemble	(const longvector& v, const std::vector<int>& doflabels);
	virtual void	assemble	(const matrix& m,     const std::vector<int>& doflabels);

public:
	longvector          eres;
	matrix              emat;
    std::map<int,int>   globalToLocal;
	
	friend class element;	
};


#endif
