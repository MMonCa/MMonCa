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
 *  dofset.h
 *  feliks
 *
 *  Created by Ignacio Romero on 6/13/09.
 *  Copyright 2009 Universidad Politécnica de Madrid. All rights reserved.
 *
 */


#ifndef _dofset_h
#define _dofset_h


#include "Math/vector.h"
#include "Math/tensor.h"
#include "idmap.h"
#include <iostream>

#ifdef WITHTBB
#include "tbb/spin_mutex.h"
#include "tbb/atomic.h"
#endif



// a enum definition to refer to the dofs in a partitioned scheme
enum partitionDofs {
	PARTITION_1,
	PARTITION_2,
	PARTITION_ALL
};



class integrator;

// class that gathers degrees of freedom that receive idential treatment
// this is a pure virtual class, for no set of degrees of freedom can
// be of no type.
class dofset
	{
		
	public:
		enum    evaluation_time{ t0, tn, tna, tn1};

#ifdef WITHTBB
		tbb::spin_mutex		assemble_mutex;
#endif		
		
		
		virtual			~dofset(){}
		virtual double& delta(const int l)=0;
		IDmap&			getIDs() 	{return ID;}
		const IDmap&	getIDs() const	{return ID;}
				int		getID(const int l=0) const {return ID.getDof(l);}
				int		getNConstrainedDofs() const;
				bool	isConstrained(const int l=0) const	{return ID.isConstrained(l);}
				bool	isUnused() const			{return ID.isUndefined();}
		virtual void	localizeSolutionIncrement(const longvector& delta);
		virtual int		ndof() const=0;
		virtual void	printDOFS(std::ostream& of=std::cout) const {}
		virtual void	recoverFromLastSolution()=0;
		virtual void	resetDelta()=0;
		virtual void	resetForce()=0;
		bool			setConstrainedDof(const int local);

	protected:
		IDmap	ID;		
	};





// vector of degrees of freedom for problems such as elasticity
class vectorDofset : public dofset
	{
		
	public:
		vectorDofset();
		
		int				ndof() const{return 3;}
		double&			delta(const int l);
		blue::ivector&		displacement(const evaluation_time when=dofset::tna);
		blue::ivector&		increment() {return dU;}
		double&			increment(const int l){return dU[l];}
		blue::ivector&		velocity(const evaluation_time when=dofset::tna);
		const blue::ivector&	velocity(const evaluation_time when=dofset::tna) const;
		blue::ivector&		acceleration(const evaluation_time when=dofset::tna);
		blue::ivector&		force()	{return F;}

		friend std::ostream&  operator<<(std::ostream &os, const vectorDofset &ds);
		void			printDOFS(std::ostream& of=std::cout) const;
		void			recoverFromLastSolution();
		void			resetDelta()	{dU.setZero();}
		void			resetForce()	{F.setZero();}
		void			setDof(const int local, const double v);	
		
		
	private:
		blue::ivector dU;
		blue::ivector U, Un, Un1;
		blue::ivector V, Vn, Vn1;
		blue::ivector A, An, An1;
		blue::ivector F;
        
        // tbb::AtomicDouble Fx, Fy, Fz;
	};



class scalarDofset : public dofset
	{
		
	public:
		scalarDofset();
		
		int			ndof() const{return 1;}
		double&		increment();
		double&		delta(const int l);
		double&		displacement(const evaluation_time when=dofset::tna);
		double&		force()	{return F;}
		double&		velocity(const evaluation_time when=dofset::tna);
		double&		acceleration(const evaluation_time when=dofset::tna);
		
		void		printDOFS(std::ostream& of=std::cout) const;
		void		recoverFromLastSolution();
		void		resetDelta()	{dp=0.0;}
		void		resetForce()	{F=0.0;}
		void		setDof(const double v);

		friend std::ostream&  operator<<(std::ostream &os, const scalarDofset &ds);


	private:
		double		dp;
		double		p, pn, pn1;
		double      V, Vn, Vn1;
		double      A, An, An1;
		double		F;
	};




class rotationDofset : public dofset
	{
	
	public:
		rotationDofset();
		
		int				ndof() const{return 3;}
		blue::ivector			director(const evaluation_time when=dofset::tna);
		blue::ivector&		force()	{return M;}
		blue::ivector&		increment()			{return dtheta;}
		double&			increment(const int l)			{return dtheta[l];}
		blue::ivector&		rotationVector();
		double&			delta(const int l);
		blue::iquaternion&	rotation(const evaluation_time when=dofset::tna);

		const blue::ivector	bodyAcceleration(const evaluation_time when=dofset::tna);
		const blue::ivector	bodyStepRotationVector();
		const blue::ivector	bodyVelocity(const evaluation_time when=dofset::tna);

		
		blue::ivector&		spatialVelocity(const evaluation_time when=dofset::tna);
		blue::ivector&		spatialAcceleration(const evaluation_time when=dofset::tna);
		blue::ivector&		spatialStepRotationVector();

		void			printDOFS(std::ostream& of=std::cout) const;
		void			recoverFromLastSolution();
		void			resetDelta()	{dtheta.setZero();}
		void			resetForce()	{M.setZero();}

		friend std::ostream&  operator<<(std::ostream &os, const rotationDofset &ds);
		
	//private:
		blue::ivector		dtheta;				// rot. increment in last iteration (spatial)
		blue::ivector		theta;				// rot. increment in current step (spatial)
		blue::iquaternion	rot0, rot, rotn, rotn1;
		blue::ivector		omega, omegan, omegan1;
		blue::ivector		alpha, alphan, alphan1;
		//blue::ivector		Omega, Omegan, Omegan1;
		//blue::ivector		Alpha, Alphan, Alphan1;
		blue::ivector		M;
	};





class directorDofset : public dofset
	{
		
	public:
		directorDofset();
		directorDofset(const blue::ivector& d);
		
		double&			delta(const int l);
		blue::ivector			director(const evaluation_time when=dofset::tna);
		blue::vector2&		force()	{return M;}
		blue::vector2&		increment()		{return dTheta;}
		int				ndof() const{return 2;}
		void			printDOFS(std::ostream& of=std::cout) const;
		blue::iquaternion&	rotation(const evaluation_time when=dofset::tna);
		blue::vector2&		bodyStepRotationVector();

		
		blue::vector2&		bodyAcceleration(const evaluation_time when=dofset::tna);
		blue::vector2&		bodyVelocity(const evaluation_time when=dofset::tna);

		blue::ivector			spatialAcceleration(const evaluation_time when=dofset::tna);
		blue::ivector			spatialVelocity(const evaluation_time when=dofset::tna);
		
				
		void			recoverFromLastSolution();
		void			resetDelta()	{dTheta.setZero();}
		void			resetForce()	{M.setZero();}
		
		friend std::ostream&  operator<<(std::ostream &os, const directorDofset &ds);


	private:
		blue::iquaternion		rot0, rot, rotn, rotn1;
		blue::vector2			Omega, Omegan, Omegan1;
		blue::vector2			Alpha, Alphan, Alphan1;
		blue::vector2			dTheta;						// convected rotation vector in last iteration
		blue::vector2			Theta;						// convected rotatton vector of current step
		blue::vector2			M;
	};




template <int N> 
class nvectorDofset : public dofset
{

	
public:
	nvectorDofset();
	nvectorDofset(const nvectorDofset<N>& d);
	
	double*		increment()			{return dU;}
	double&		delta(const int l)	{return dU[l];}
	double*		force()				{return F;}
	int			ndof() const		{return N;}
	void		printForce(std::ostream& of=std::cout)
	{
		for (int i=0; i<N; i++) std::cout << F[i];
	}
	
	void		recoverFromLastSolution()
	{
		for (int i=0; i<N; i++)
		{
			U[i] = Un1[i] = Un[i];
			V[i] = Vn1[i] = Vn[i];
			A[i] = An1[i] = An[i];
		}
		resetDelta();
	}
	
	void		resetDelta()	
	{
		for (int i=0; i<N; i++)	delta(i) = 0.0;
	}
	
	void		resetForce()	
	{
		for (int i=0; i<N; i++)	F[i] = 0.0;
	}
	
	friend std::ostream&  operator<<(std::ostream &os, const nvectorDofset &ds)
	{
		return os;
	}
	
	
	
private:
	double		dU[N];
	double		U[N], Un[N], Un1[N];
	double		V[N], Vn[N], Vn1[N];
	double		A[N], An[N], An1[N];
	double		F[N];
	
	
};



#endif

