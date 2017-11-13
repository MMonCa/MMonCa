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
 * shapefun.cpp
 */
#include "feshapefun.h"
#include "Geometry/quadpoint.h"


#include <cmath>
#include <iostream>
#include <iomanip>
#include "Elements/element.h"
#include "Model/elmtface.h"
#include "Io/message.h"
#include "Math/matrix.h"
#include "Math/vector.h"
#include "Io/logger.h"
#include "Math/tensor.h"


using namespace blue;


// this function should only be called by a derived constructor
FEshapefunbuilder :: FEshapefunbuilder()
{}


// constructor for the shapefunctions blocks. 
// the blocks are always associated with an element, it is not just a collection of shapefun's
FEshapefunbuilder :: FEshapefunbuilder(const element& e)
    :   shp( e.getNNodes() )
{
	for (int a=0; a<shp.size(); a++) 
        shp[a].parentNode = const_cast<node*>(&(e.getNode(a)));
}



// constructor for the shapefunctions blocks. 
// the blocks are always associated with an element, it is not just a collection of shapefun's
FEshapefunbuilder :: FEshapefunbuilder(const element& e, const int nshp)
:
    shp( nshp)
{ 
	for (int a=0; a<nshp; a++) 
        shp[a].parentNode = const_cast<node*>(&(e.getNode(a)));
}




/* This is the main shape function file for 1D problems (beam, rods, ...). The function works
even if it is posed in 2 or 3d, that is for rods and beams in 2 and three dimensional problems.*/
void FEshapefunbuilder :: eval1d(const quadpoint& gp, double& j, const bool dd, const int push)
{
    // reference element coordinate
    double xi ( gp.coor[0] );
	double l;
	ivector vtmp;
    
    switch( shp.size() )
    {
        // linear element with 2 nodes
        case(2):
            /* length */
			vtmp = shp[1].parentNode->getReferenceCoordinates() - shp[0].parentNode->getReferenceCoordinates();
			l = vtmp.norm();
            
            /* shape function value */
            shp[0].value = 0.5*(1.0 - xi);
            shp[1].value = 0.5*(1.0 + xi);
            
            /* constant jacobian L_real/L_ref */
            j = 0.5*l;
            
            if (push == PUSH)
            {
                /* cartesian derivatives */
                shp[0].dxyz[0] = -1.0/l;
                shp[1].dxyz[0] =  1.0/l;
            }
            else
            {
				/* natural derivatives */
                shp[0].dxyz[0] = -0.5;
                shp[1].dxyz[0] =  0.5;
			}
                
            break;
            
            
        /* the quadratic rod. Node numbers:   1----3----2 */
        case(3):
            
            /* shape function values */
            shp[0].value = 0.5*xi*(xi - 1.0);
            shp[1].value = 0.5*xi*(xi + 1.0);
            shp[2].value = 1.0 - xi*xi;
            
            /* natural derivatives */
            shp[0].dxyz[0] = xi - 0.5;
            shp[1].dxyz[0] = xi + 0.5;
            shp[2].dxyz[0] = -2.0*xi;
            
            /* jacobian */
			vtmp.setZero();
			for (int a=0; a<shp.size(); a++)	
                vtmp += shp[a].dxyz[0]*shp[a].parentNode->getReferenceCoordinates();
				j = vtmp.norm();
            
            /* derivatives w/r to global coordinates if needed, otherwise, return the natural ones */
            if (push == PUSH)
                for (int a=0; a<shp.size(); a++) shp[a].dxyz[0] = shp[a].dxyz[0]/j;
                
				/* second derivatives should go here */
                
				break;
    }
    
    /* reset all other values to zero */
    for (int a=0; a<shp.size(); a++)
    {
        shp[a].dxyz[1] = shp[a].dxyz[2] = 0.0;
        shp[a].secondDflag = false;
    }
}




void FEshapefunbuilder :: eval2d(const quadpoint& gp, double& j, const bool dd, const int push)
{
    switch ( size() )
    {        
        case(3): /* triangle */
            evalCST(gp, j, dd, push);
            break;
            
        case(4): // four node quad
            evalQUAD(gp, j, dd, push);
            break;
            
        case(6): //quadratic triangle 
            evalTRIQUAD(gp, j, dd, push);
            break;
            
        default:
            logger::warnings << "2d Shape functions not programmed for " << shp.size() << " nodes.";
    }
}





void FEshapefunbuilder :: eval3d(const quadpoint& gp, double& jacobian, const bool secondd, const int push)
{
    switch (shp.size() )
    {
        case(4):
            evalTET(gp, jacobian, secondd, push);
            break;
            
        case(8):
            evalBRICK(gp, jacobian, secondd, push);
            break;
			
		case(10):
            evalTET10(gp, jacobian, secondd, push);
            break;
			
        default:
            logger::warnings << "3d Shape functions not programmed for " << shp.size() << "nodes.";
    }
}




/*
 * the brick has the nodes in the following order:
 *
 *              /\ y                   /\ y 
 *         2 ----|----- 1         6 ----|----- 5
 *         |     |      |         |     |      |
 *         |   (z=-1)------>x     |   (z=+1)---|--> x
 *         |            |         |            |
 *         3 ---------- 4         7 ---------- 8 
 *
 */
void FEshapefunbuilder :: evalBRICK(const quadpoint& gp, double& jacobian, const bool secondd, const int push)
{
    double xi, eta, zeta;
    static double xisgn[]={+1.0, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0, +1.0};
    static double etsgn[]={+1.0, +1.0, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0};
    static double zesgn[]={-1.0, -1.0, -1.0, -1.0, +1.0, +1.0, +1.0, +1.0};    
    
    // reference element coordinates
    xi    =  gp.coor[0];   
    eta   =  gp.coor[1];
    zeta  =  gp.coor[2];
        
    // the value of the shape functions at the quad point
    for (int a=0; a<8; a++) 
        shp[a].value = 0.125 * (1.0 + xisgn[a]*xi) * (1.0 + etsgn[a]*eta) *(1.0 + zesgn[a]*zeta);
    
    
    /* the shape function natural derivatives */
	ivector Nd[8];
    for (int a=0; a<8; a++)
    {
        Nd[a][0] = 0.125 *        xisgn[a]      * (1.0 + etsgn[a]*eta) *(1.0 + zesgn[a]*zeta);
        Nd[a][1] = 0.125 * (1.0 + xisgn[a]*xi ) *        etsgn[a]      *(1.0 + zesgn[a]*zeta);
        Nd[a][2] = 0.125 * (1.0 + xisgn[a]*xi ) * (1.0 + etsgn[a]*eta )*       zesgn[a]      ;  
    }
    
    
    /* jacobian matrix of the mapping from the biunit square to the current element 
	 [ x,xi x,eta x,zeta; y,xi y,eta y,zeta; z,xi z,eta z,zeta] */
	itensor F;
	F.setZero();
	for (int a=0; a<8; a++)
		F.addDyadic(shp[a].parentNode->getReferenceCoordinates(),	Nd[a]);
	
	// inverse jacobian matrix and determinant
	itensor Finv(F);
	jacobian = Finv.invert();
	
    if (jacobian < 0.0) 
		Message("\n Warning, negative jacobian in BRICK");
    
    
    // global derivatives of the shape function -- push forward of the
	// natural derivatives with f^(-t) */
    if (push == PUSH) 
		for (int a=0; a<8; a++) 
			shp[a].dxyz = Finv.transpose() * Nd[a]; 
	
	
	if ( secondd == true)
	{
		/* form second derivatives if needed, pushing forward tensor N,xi_i xi_j
		 as in D^2 N = F^(-T) d^2 N F^(-1) */
		double d2x[3][3][3];
		double tmp33[3][3];
		double Ndd[8][3][3];
		
		
		// Construct natural hessians of shape functions: N,xi_i xi_j
		for (int a=0; a<8 ; a++)
		{
			Ndd[a][0][0] = 0.0;
			Ndd[a][0][1] = 0.125*       xisgn[a]    *       etsgn[a]      *(1.0 + zesgn[a]*zeta);
			Ndd[a][0][2] = 0.125*       xisgn[a]    *(1.0 + etsgn[a]*eta) *       zesgn[a];
			
			Ndd[a][1][0] = 0.125*       xisgn[a]    *       etsgn[a]      *(1.0 + zesgn[a]*zeta);
			Ndd[a][1][1] = 0.0;
			Ndd[a][1][2] = 0.125*(1.0 + xisgn[a]*xi)*       etsgn[a]      *       zesgn[a];
			
			Ndd[a][2][0] = 0.125*       xisgn[a]    *(1.0 + etsgn[a]*eta) *       zesgn[a];
			Ndd[a][2][1] = 0.125*(1.0 + xisgn[a]*xi)*       etsgn[a]      *       zesgn[a]; 
			Ndd[a][2][2] = 0.0;
		}
		
		
		/* construct second derivatives of coordinates x w/r to xi: d^2 x_i / d xi_j d xi_k */
		for (int i=0; i<3 ; i++)
			for (int j=0; j<3 ; j++)
				for (int k=0; k<3 ; k++)
				{
					d2x[i][j][k] = 0.0;
					for (int a=0; a<8; a++) 
						d2x[i][j][k] += Ndd[a][j][k]*shp[a].parentNode->getReferenceCoordinates()[i];
				}
		
		// loop over the shape functions to obtain second derivatives
		for (int a=0; a<8; a++)
		{
			
			// First part: push forward N,xi xi as a covariant-covariant tensor
			for (int i=0; i<3; i++) 
				for (int j=0; j<3; j++)
				{
					tmp33[i][j] = 0.0;
					for (int k=0; k<3; k++) 
						tmp33[i][j] += Finv(k,i)*Ndd[a][k][j];
				}
			
			for (int i=0; i<3; i++) 
				for (int j=0; j<3; j++)
				{
					shp[a].dd(i,j) = 0.0;
					for (int k=0; k<3; k++) 
						shp[a].dd(i,j) += tmp33[i][k] * Finv(k,j);
				}
			
			// second part, push forward of rank-three tensor
			for (int i=0; i<3; i++) 
				for (int j=0; j<3; j++)
					for (int k=0; k<3; k++)
						for (int l=0; l<3; l++)
							for (int m=0; m<3; m++)
								shp[a].dd(i,j) -= shp[a].dxyz[m]*Finv(l,i)*Finv(k,j)*d2x[m][l][k];
			
			shp[a].secondDflag = true;
		}
	}
}




/* shape functions for quadratic tetrahedra. There is a choice to be made here because there is
* not a unique way to number the nodes in a quadratic tet. Different references use different
* conventions. 
*
* 
* The following faces are defined (outward normal corresponds to outward normal to the screen)
* (This is the numbering of FEAP and GID)
*
*  1               1               1            4
*  |\              |\              |\           |\
*  7 5             8 7             5  8         9 10
*  |  \            |  \            |   \        |  \
*  3-6-2           4-10-3          2-9-4        2-6-3
*/
void FEshapefunbuilder :: evalTET10(const quadpoint& gp, double& jacobian, const bool dd, const int push)
{
	double r, s, t, u;
	extern   double global_macheps;
	
	// reference element coordinates
	r   = gp.coor[0];
    s   = gp.coor[1];
    t   = gp.coor[2];
    u   = 1.0 - r - s - t;
	
	
	// the value of the shape functions at the gauss point (Zienkiewicz & Taylor)
    shp[0].value = u * (2.0*u-1.0);
	shp[1].value = r * (2.0*r-1.0);
    shp[2].value = s * (2.0*s-1.0);
	shp[3].value = t * (2.0*t-1.0);
	
    shp[4].value = 4.0*u*r;
    shp[5].value = 4.0*r*s;
    shp[6].value = 4.0*s*u;
	
	shp[7].value = 4.0*u*t;
	shp[8].value = 4.0*r*t;
	shp[9].value = 4.0*s*t;
	
	
	// There is a faster way to do this, but we take advantage of the routine already
	// programmed for the hexahedra
    // natural derivatives
	ivector Nd[10];
	double u4 = 4.0*u, r4 = 4.0*r, s4 = 4.0*s, t4 = 4.0*t;
	
	Nd[0][0] = Nd[0][1] = Nd[0][2] = 1.0 - u4;
	
	Nd[1][1] = Nd[1][2] = 0.0;
	Nd[1][0] = r4 - 1.0;
	
	Nd[2][0] = Nd[2][2] = 0.0;
	Nd[2][1] = s4 - 1.0;
	
	Nd[3][0] = Nd[3][1] = 0.0;
	Nd[3][2] = t4 - 1.0;
	
	Nd[4][0] = -r4 + u4;
	Nd[4][1] = -r4;
	Nd[4][2] = -r4;
	
	Nd[5][0] = s4;
	Nd[5][1] = r4;
	Nd[5][2] = 0.0;
	
	Nd[6][0] = -s4;
	Nd[6][1] = -s4 + u4;
	Nd[6][2] = -s4;
	
	Nd[7][0] = -t4;
	Nd[7][1] = -t4;
	Nd[7][2] = -t4 + u4;
	
	Nd[8][0] = t4;
	Nd[8][1] = 0.0;
	Nd[8][2] = r4;
	
	Nd[9][0] = 0.0;
	Nd[9][1] = t4;
	Nd[9][2] = s4;
	
	
	/* jacobian matrix of the mapping from the biunit square to the current element 
        [ x,xi x,eta x,zeta; y,xi y,eta y,zeta; z,xi z,eta z,zeta] */
	itensor F;
	F.setZero();
	for (int a=0; a<10; a++) 
        F.addDyadic( shp[a].parentNode->getReferenceCoordinates() , Nd[a] );
    
	
	/* inverse jacobian matrix and determinant */
	itensor Finv(F);
	jacobian = Finv.invert();
    if (jacobian < 20*global_macheps) Message("\n Warning, negative or zero jacobian in TET10");
    
    
    /*global derivatives of the shape function -- push forward of the
        natural derivatives with f^(-t) */
    if (push) for (int a=0; a<10; a++) shp[a].dxyz = Finv.transpose() * Nd[a];
	
    return;
}






/*  shape functions for constant strain triangle. 
*  implementation using area coordinates
*/
void FEshapefunbuilder :: evalCST(const quadpoint& gp, double& jacobian, const bool secondd, const int push)
{
    int    i,j,a;
    double s,t,u;
    double df[2][2], Finv[2][2];	
	
    // reference element coordinates
    s   = gp.coor[0];   
    t   = gp.coor[1];
    u   = 1.0 - s - t;
    
    /* the value of the shape functions at the quad point */
    shp[0].value = u;
    shp[1].value = s;
    shp[2].value = t;
    
    
    if (push)
    {
		const ivector& coor0( shp[0].parentNode->getReferenceCoordinates() );
		const ivector& coor1( shp[1].parentNode->getReferenceCoordinates() );
		const ivector& coor2( shp[2].parentNode->getReferenceCoordinates() );
        
		/* jacobian matrix [ x,xi   x,eta; y,xi   y,eta] and determinant */
        df[0][0] = coor1[0] - coor0[0];    /* x1-x0*/
        df[0][1] = coor2[0] - coor0[0];    /* x2-x0*/
        df[1][0] = coor1[1] - coor0[1];    /* y1-y0*/
        df[1][1] = coor2[1] - coor0[1];    /* y2-y0*/
        jacobian = df[0][0]*df[1][1]-df[0][1]*df[1][0];
        if (jacobian < 0.0) Message("\n Warning, negative jacobian in CST");
        
        /* global derivatives of the shape function -- push forward of the
            natural derivatives with f^(-t) */
        /* compute f^(-1) */
        Finv[0][0] =  df[1][1]/ (jacobian);
        Finv[1][1] =  df[0][0]/ (jacobian);
        Finv[0][1] = -df[0][1]/ (jacobian);
        Finv[1][0] = -df[1][0]/ (jacobian);
        
        shp[0].dxyz[0] = -(Finv[0][0] + Finv[1][0]);
        shp[0].dxyz[1] = -(Finv[0][1] + Finv[1][1]);
        shp[0].dxyz[2] =  0.0;
        
        shp[1].dxyz[0] =  Finv[0][0];
        shp[1].dxyz[1] =  Finv[0][1];
        shp[1].dxyz[2] =  0.0;
        
        shp[2].dxyz[0] =  Finv[1][0];
        shp[2].dxyz[1] =  Finv[1][1];
        shp[2].dxyz[2] =  0.0;
    }
    
    // derivatives with respect to isoparametric coordinates
    else
    {
        shp[0].dxyz[0] = -1.0;
        shp[0].dxyz[1] = -1.0;
        shp[0].dxyz[2] =  0.0;
        
        shp[1].dxyz[0] =  1.0;
        shp[1].dxyz[1] =  0.0;
        shp[1].dxyz[2] =  0.0;
        
        shp[2].dxyz[0] =  0.0;
        shp[2].dxyz[1] =  1.0;
        shp[2].dxyz[2] =  0.0;
    }      
    
    
    /* calculate second derivatives. Are zero for linear element */
    if (secondd == 1)
    {
        for (a=0; a<3; a++)
        {
            for (i=0; i<2; i++)
                for (j=0; j<2; j++)
                    shp[a].dd(i,j) = 0.0;
            
            shp[a].secondDflag = true;
        }
    }
}




/* note: do not change node ordering. It affects Bathe-Dvorkin interpolations
* the brick has the nodes in the following order:
*
*              /\ y       
*         2 ----|----- 1
*         |     |      |
*         |     --------->x
*         |            |
*         3 ---------- 4
*
*/
void FEshapefunbuilder :: evalQUAD(const quadpoint& gp, double& jacobian, const bool secondd, const int push)
{
    int    a, i, j , k, l, m;
    double xi, eta;
    double F[2][2] ,  Finv[2][2];
    double xisgn[]={+1.0 , -1.0 , -1.0 , +1.0};
    double etsgn[]={+1.0 , +1.0 , -1.0 , -1.0};
    double Nd[4][2];
    double tmp22[2][2];
    
    
    /* reference element coordinates */
    xi  = gp.coor[0];   
    eta = gp.coor[1];
    
    /* the value of the shape functions at the quad point */
    for (a=0; a<4; a++)
        shp[a].value = 0.25 * (1.0 + xisgn[a]*xi) * (1.0 + etsgn[a]*eta);
    
    /* the shape function natural derivatives */
    for (a=0; a<4; a++)
    {
        Nd[a][0] = 0.25 *      xisgn[a]     * (1.0+etsgn[a]*eta);
        Nd[a][1] = 0.25 * (1.0+xisgn[a]*xi) *      etsgn[a]     ;
    }
    
    
    /* global derivatives of the shape function -- push forward of the
        natural derivatives with f^(-t) */
	ivector coor[4];
	for (a=0; a<4; a++) coor[a] = shp[a].parentNode->getReferenceCoordinates();
	
	if (push == PUSH)
    {
        /* jacobian matrix [ x,xi x,eta; y,xi y,eta] and determinant */
        for (i=0; i<2; i++)
        {
            for (j=0; j<2; j++)
            {
                F[i][j] = 0.0;
                for (a=0; a<4; a++) F[i][j] += Nd[a][j]*coor[a][i];
            }
        }
        
        jacobian = F[1][1]*F[0][0]-F[0][1]*F[1][0];
        if (jacobian < 0.0) 
			Message("\n Warning, negative jacobian in QUAD");
        
        Finv[0][0] =  F[1][1]/ (jacobian);
        Finv[0][1] = -F[0][1]/ (jacobian);
        Finv[1][0] = -F[1][0]/ (jacobian);
        Finv[1][1] =  F[0][0]/ (jacobian);
        
        for (a=0; a<4; a++)
		{
			for (i=0; i<2; i++)
            {
                shp[a].dxyz[i] = 0.0;
                for (j=0; j<2; j++) shp[a].dxyz[i] += Finv[j][i]*Nd[a][j];
            }
			shp[a].dxyz[2] = 0.0;
		}
			
    }
        
        /* natural derivatives */
        else
        {
            for (a=0; a<4; a++)
			{
				for (i=0; i<2; i++)
                    shp[a].dxyz[i] = Nd[a][i];
				shp[a].dxyz[2] = 0.0;
			}
        }
        
        
        /* second derivatives w/r to isoparametric coordinates */
        if (secondd == SECONDD)
        {
			double Ndd[4][2][2];
			double d2x[2][2][2];

            // Construct natural hessian N,xi_i xi_j
            for (a=0; a<4 ; a++)
            {
                Ndd[a][0][0] = 0.0;
                Ndd[a][0][1] = 0.25*  xisgn[a] * etsgn[a];
                
                Ndd[a][1][0] = 0.25*  xisgn[a] * etsgn[a];
                Ndd[a][1][1] = 0.0;
            }
        
        
			// form second derivatives if needed, pushing forward tensor N,xi_i xi_j
            // as in D^2 N = F^(-T) d^2 N F^(-1)
			if (push == PUSH)		
			{
            /* construct second derivatives of coordinates x w/r to xi: d^2 x_i / d xi_j d xi_k */
            for (i=0; i<2 ; i++)
                for (j=0; j<2 ; j++)
                    for (k=0; k<2 ; k++)
                    {
                        d2x[i][j][k] = 0.0;
                        for (a=0; a<4; a++) d2x[i][j][k] += Ndd[a][j][k]*coor[a][i];
                    }
                        
                        
                        /* loop over the shape functions to obtain second derivatives */
                        for (a=0; a<4; a++)
                        {
                            for (k=0; k<2; k++)
                                for (l=0; l<2; l++)
                                {
                                    tmp22[k][l] = Ndd[a][k][l];
                                    for(m=0; m<2; m++)  tmp22[k][l] -= shp[a].dxyz[m]*d2x[m][k][l];
                                }
                                    
                                    for (i=0; i<2; i++)
                                        for (j=0; j<2; j++)
                                        {
                                            shp[a].dd(i,j) = 0.0;
                                            
                                            for (k=0; k<2; k++)
                                                for (l=0; l<2; l++)
                                                    shp[a].dd(i,j) += tmp22[k][l] * Finv[k][j] * Finv[l][i];
                                        }
                                            
                                            shp[a].secondDflag = true;
                        }
        }    
            else
                for (a=0; a<4; a++) shp[a].secondDflag = false; 
		}
}




/* shape functions for tetrahedra.
* Looking from node 4, the other 3 read 1, 2, 3 counterclockwise
* In gp come the 3 first volume coordinates. The 4th is 1 - (the sum of the other 3).
*/
void FEshapefunbuilder :: evalTET(const quadpoint& gp, double& jacobian, const bool dd, const int push)
{
    int    perm[]={0, 1, 2, 3, 0, 1, 2};
    double xi[4], a[4][4], det, idet, x[4], y[4], z[4];
    int i, j, k, l;
    
    /* shortnames for nodal coordinates */
    for (i=0; i<4; i++)
    {
		const ivector& c( shp[i].parentNode->getReferenceCoordinates() );
        x[i] = c[0];
        y[i] = c[1];
        z[i] = c[2];
    }
	
	
    /* computation of the minors of matrix A */
    for (i=0; i<4; i++)
    {
        j = perm[i+1];
        k = perm[i+2];
        l = perm[i+3];
        
        /* first, ignore the signs */
        
        /* we can do without these minors
            a[i][0] = x[j]*(y[k]*z[l] - y[l]*z[k])
            +     x[k]*(y[l]*z[j] - y[j]*z[l])
            +     x[l]*(y[j]*z[k] - y[k]*z[j]);*/
        
        a[i][1] =  y[k]*z[l] - y[l]*z[k]
            +      y[l]*z[j] - y[j]*z[l]
            +      y[j]*z[k] - y[k]*z[j];
        
        a[i][2] =  x[k]*z[l] - x[l]*z[k]
            +      x[l]*z[j] - x[j]*z[l]
            +      x[j]*z[k] - x[k]*z[j];
        
        a[i][3] =  x[k]*y[l] - x[l]*y[k]
            +      x[l]*y[j] - x[j]*y[l]
            +      x[j]*y[k] - x[k]*y[j];
        
        /* now correct the signs */
        a[i][j] = -a[i][j];
        a[i][l] = -a[i][l];
    }
    
    /* now the determinant, using the first column of A */
    det = x[0]*a[0][1] + x[1]*a[1][1] + x[2]*a[2][1] + x[3]*a[3][1];
    
    /* the volume coordinates */
    for (i=0; i<3; i++) xi[i] = gp.coor[i];
    xi[3] = 1.0 - xi[0] - xi[1] - xi[2];
    
    /* the shape functions */
    idet = 1.0/det;
    for (i=0; i<4; i++)
    {
        shp[i].value   = xi[i];
        shp[i].dxyz[0] = a[i][1]*idet;
        shp[i].dxyz[1] = a[i][2]*idet;
        shp[i].dxyz[2] = a[i][3]*idet;
    }
	
    /* the volume of the reference tetrahedra is 1/6 */
    jacobian = det/6.0;
    
    return;
}





/* shape functions for quadratic triangle. The node ordering is as follows:
*      3
*      | \
*      6   5
*      |    \
*      1--4--2
*
* see Bathe pg 373 for details
*/
void FEshapefunbuilder :: evalTRIQUAD(const quadpoint& gp, double& jacobian, const bool dd, const int push)
{
    double s, t, u;
    int i, j;
    
    // reference element coordinates
    s   = gp.coor[0];
    t   = gp.coor[1];
    u   = 1.0 - s - t;
    
    
    /* the value of the shape functions at the gauss point */
    shp[3].value = 4.0*s*u;
    shp[4].value = 4.0*s*t;
    shp[5].value = 4.0*t*u;
    shp[0].value = u - 0.5*shp[3].value                    - 0.5*shp[5].value;
    shp[1].value = s - 0.5*shp[3].value - 0.5*shp[4].value;
    shp[2].value = t                    - 0.5*shp[4].value - 0.5*shp[5].value;
    
    /* the natural derivatives */
	ivector Nd[6];
    Nd[0][0] = -3.0 + 4.0*s + 4.0*t;
    Nd[1][0] = -1.0 + 4.0*s;
    Nd[2][0] =  0.0;
    Nd[3][0] =  4.0 - 8.0*s - 4.0*t;
    Nd[4][0] =  4.0*t;
    Nd[5][0] = -4.0*t;
    
    Nd[0][1] = -3.0 + 4.0*s + 4.0*t;
    Nd[1][1] =  0.0;
    Nd[2][1] = -1.0         + 4.*t;
    Nd[3][1] = -4.0*s;
    Nd[4][1] =  4.0*s;
    Nd[5][1] =  4.0 - 4.0*s - 8.0*t;
    
    
    /* global derivatives of the shape function -- push forward of the
        natural derivatives with f^(-t) */
    if (push == PUSH)
    {
		/* jacobian matrix of the mapping from the biunit square to the current element 
        [ x,xi x,eta x,zeta; y,xi y,eta y,zeta; z,xi z,eta z,zeta] */
		itensor F;
		F.setZero();
		for (int a=0; a<6; a++) F.addDyadic( shp[a].parentNode->getReferenceCoordinates() , Nd[a] );
		
        jacobian = F(1,1)*F(0,0)-F(0,1)*F(1,0);
        if (jacobian < 0.0) Message("\n Warning, negative jacobian in TRIQUAD");
        
		itensor Finv;
        Finv(0,0) =  F(1,1) / (jacobian);
        Finv(0,1) = -F(0,1) / (jacobian);
        Finv(1,0) = -F(1,0) / (jacobian);
        Finv(1,1) =  F(0,0) / (jacobian);
        
        for (int a=0; a<6; a++)
            for (i=0; i<2; i++)
            {
                shp[a].dxyz[i] = 0.0;
                for (j=0; j<2; j++) shp[a].dxyz[i] += Finv(j,i)*Nd[a][j];
            }
    }
        
        /* natural derivatives */
        else
        {
            for (int a=0; a<6; a++)
                for (i=0; i<2; i++)
                    shp[a].dxyz[i] = Nd[a][i];
        }
        
        
        /* 
		//second derivatives w/r to isoparametric coordinates 
		//double tmp22[2][2], Ndd[6][2][2];
		//double d2x[2][2][2];
	
			if (dd == SECONDD)
        {
				// Construct natural hessian N,xi_i xi_j
				Ndd[0][0][0] = Ndd[0][0][1] = Ndd[0][1][0] = Ndd[0][1][1] = 4.0;
				Ndd[1][0][0] = 4.0;
				Ndd[2][1][1] = 4.0;
				Ndd[3][0][0] = -8.0; Ndd[3][1][0] = Ndd[3][0][1] = -4.0;
				Ndd[4][0][1] = Ndd[4][1][0] = 4.0;
				Ndd[5][0][1] = Ndd[5][1][0] = -4.0; Ndd[5][1][1] = -8.0;
				
				Ndd[1][0][1] = Ndd[1][1][0] = Ndd[1][1][1] = 0.0;
				Ndd[2][0][1] = Ndd[2][1][0] = Ndd[2][0][0] = 0.0;
				Ndd[3][1][1] = 0.0;
				Ndd[4][0][0] = Ndd[4][1][1] = 0.0;
				Ndd[5][0][0] = 0.0;
				
				for (a=0; a<6; a++) shp[a].secondDflag = true;
        }
        
        // form second derivatives if needed, pushing forward tensor N,xi_i xi_j
        //    as in D^2 N = F^(-T) d^2 N F^(-1) 
			if (dd == SECONDD && push == PUSH)
        {
				// construct second derivatives of coordinates x w/r to xi: d^2 x_i / d xi_j d xi_k 
				for (i=0; i<2 ; i++)
					for (j=0; j<2 ; j++)
						for (k=0; k<2 ; k++)
						{
							d2x[i][j][k] = 0.0;
							for (a=0; a<6; a++) d2x[i][j][k] += Ndd[a][j][k]*coor->data[a][i];
						}
							
							
							// loop over the shape functions to obtain second derivatives
							for (a=0; a<6; a++)
							{
								for (k=0; k<2; k++)
									for (l=0; l<2; l++)
									{
										tmp22[k][l] = Ndd[a][k][l];
										for(m=0; m<2; m++)  tmp22[k][l] -= shp[a].dxyz[m]*d2x[m][k][l];
									}
										
										for (i=0; i<2; i++)
											for (j=0; j<2; j++)
											{
												shp[a].dd(i,j) = 0.0;
												
												for (k=0; k<2; k++)
													for (l=0; l<2; l++)
														shp[a].dd(i,j) += tmp22[k][l] * Finv[k][j] * Finv[l][i];
											}
							}
        }
		*/
}






/* this is as above, for for rods or beams in a space of dimension 1, 2 or 3
void shapefun :: EvaluateCurveShapeFunctions(matrix coor, quadpoint gp, shapefun *shp, double *jacobian)
{
    int nnode;
    nnode = coor->rows;
    Eval1dShapeFunctions(nnode, coor, gp, shp, jacobian, NOSECONDD, PUSH);
}
*/




/* This function evaluates the shape functions and their natural derivatives (with
respect to the coordinates in the reference element
*/
void FEshapefunbuilder :: evaluateRef(const quadpoint& gp)
{
	double dummy;
    evaluate(gp, dummy, NOPUSH);
}




/* This function computes the shape functions, its derivatives and its second derivatives too
*/
void FEshapefunbuilder :: evaluateDD(const quadpoint& gp, double& jacobian)
{
    evaluate(gp, jacobian, SECONDD, PUSH);
}



FEshapefunbuilder& FEshapefunbuilder ::	operator=(const FEshapefunbuilder &sh)
{
	shp.clear();
	
	this->shp.reserve(sh.size());
	for (int a=0; a<sh.size(); a++) 
	{
		shp.push_back( sh[a] );
	}
	return *this;
}



// constructor for the shapefunctions blocks. 
// the blocks are always associated with an element, it is not just a collection of shapefun's
// theOwner is just a ressource manager that takes care of eliminating the shapefunctions contained in the object
surfaceShapefuns :: surfaceShapefuns(const element& e)
: FEshapefunbuilder(e)
{ 
}



surfaceShapefuns :: surfaceShapefuns(const element& e, const int nshp)
: FEshapefunbuilder(e, nshp)
{ 
}




/* this is as above, for for shells or similar */
void surfaceShapefuns :: evaluate(const quadpoint& gp, double& jacobian)
{ 
    eval2d(gp, jacobian, NOSECONDD, PUSH);
}


void surfaceShapefuns :: evaluateRef(const quadpoint& gp)
{
	double tmp;
    eval2d(gp, tmp, NOSECONDD, NOPUSH);
}




/* 
   Given a surface element, the isoparametric coordinates define their corresponding directional
   derivatives and the tangent vectors to the element along the isoparametric curves. These
	might be inappropriate for certain applications where we need the orthogonal tangent vectors,
	their directions and the shape functions derivative along these curves. 
 
	Given a point on a surface with its shape functions and isoparametric derivatives, and
	the tangent vectors (a0,a1) to the surface along the isoparametric curves, we find the
	pair of unit vectors (e0,e1) that span the same plane as (a0,a1) and are the closest to (a0,a1)
	among all the posible basis for the tangent plane.
 
	Also the shape function derivatives along those same directions are computed
 
*/
void surfaceShapefuns :: evaluateOrthogonal(const ivector& a0, const ivector& a1, surfaceShapefuns& shpx, double& sj)
{
    double x_xi[2][2], xi_x[2][2];
    double det, idet;
    
    // compute the normal to the pair a_xi, a_eta
	ivector nu = a0.cross(a1);
	sj = nu.norm();
	nu.normalize();
	
    //  This is the unique rotation that maps E3 into nu.
    //  The first two columns are e1 & e2
	irotation rot;
	rot.beRotationWithoutDrill(nu);
	ivector e0(rot.col(0)), e1(rot.col(1));
	
	
    // the matrix [ d(x_j) / (xi_b) ] = [e_j * a_b]
    x_xi[0][0] = e0.dot(a0);
    x_xi[0][1] = e0.dot(a1);
    x_xi[1][0] = e1.dot(a0);
    x_xi[1][1] = e1.dot(a1);
	
    // the matrix [ d(xi_b) / d(x_j) ], as the inverse of [ d(x_j) / (xi_b) ]
    det = x_xi[0][0]*x_xi[1][1] - x_xi[0][1]*x_xi[1][0];
    if (det <= 0.0)
	{
		logger :: warnings << "\nNon-positive jacobian in membrane";
		Dzero(xi_x[0], 4);
	}
    else
    {
        idet = 1.0/det;
        xi_x[0][0] =  x_xi[1][1] * idet;
        xi_x[1][1] =  x_xi[0][0] * idet;
        xi_x[0][1] = -x_xi[0][1] * idet;
        xi_x[1][0] = -x_xi[1][0] * idet;
    }
	
	
    /* finally, compute the cartesian shape function derivatives pushing forward
        the isoparametric ones with the inverse transposed jacobian */
    for (int a=0; a<shp.size(); a++)
    {
        shpx[a].value   = shp[a].value;
        shpx[a].dxyz[0] = shp[a].dxyz[0] * xi_x[0][0] + shp[a].dxyz[1] * xi_x[1][0];
        shpx[a].dxyz[1] = shp[a].dxyz[0] * xi_x[0][1] + shp[a].dxyz[1] * xi_x[1][1];
		shpx[a].dxyz[2] = 0.0;
    }
}







// constructor for the shapefunctions blocks. 
// the blocks are always associated with an element, it is not just a collection of shapefun's
// theOwner is just a ressource manager that takes care of eliminating the shapefunctions contained in the object
faceShapefuns :: faceShapefuns(const elmtface& aFace) 
{ 
    shp.resize( aFace.getNNodes() );
	for (int a=0; a<shp.size(); a++) 
        shp[a].parentNode =  const_cast<node*>(&(aFace.getNode(a)));
}


void faceShapefuns :: evaluate(const quadpoint& gp, double& j)
{
	eval2d(gp, j, NOSECONDD, PUSH);
}




void faceShapefuns :: evaluateRef(const quadpoint& gp)
{
	double tmp;
    eval2d(gp, tmp, NOSECONDD, NOPUSH);
}




/* 
 Given a surface element, the isoparametric coordinates define their corresponding directional
 derivatives and the tangent vectors to the element along the isoparametric curves. These
 might be inappropriate for certain applications where we need the orthogonal tangent vectors,
 their directions and the shape functions derivative along these curves. 
 
 Given a point on a surface with its shape functions and isoparametric derivatives, and
 the tangent vectors (a0,a1) to the surface along the isoparametric curves, we find the
 pair of unit vectors (e0,e1) that span the same plane as (a0,a1) and are the closest to (a0,a1)
 among all the posible basis for the tangent plane.
 
 Also the shape function derivatives along those same directions are computed
 
 */
void faceShapefuns :: evaluateOrthogonal(const ivector& a0, const ivector& a1, faceShapefuns& shpx, double& sj)
{
    double x_xi[2][2], xi_x[2][2];
    double det, idet;
    
    // compute the normal to the pair a_xi, a_eta
	ivector nu = a0.cross(a1);
	sj = nu.norm();
	nu.normalize();
	
    //  This is the unique rotation that maps E3 into nu.
    //  The first two columns are e1 & e2
	irotation rot;
	rot.beRotationWithoutDrill(nu);
	ivector e0(rot.col(0)), e1(rot.col(1));
	
	
    // the matrix [ d(x_j) / (xi_b) ] = [e_j * a_b]
    x_xi[0][0] = e0.dot(a0);
    x_xi[0][1] = e0.dot(a1);
    x_xi[1][0] = e1.dot(a0);
    x_xi[1][1] = e1.dot(a1);
	
    // the matrix [ d(xi_b) / d(x_j) ], as the inverse of [ d(x_j) / (xi_b) ]
    det = x_xi[0][0]*x_xi[1][1] - x_xi[0][1]*x_xi[1][0];
    if (det <= 0.0)
	{
		logger :: warnings << "\nNon-positive jacobian in membrane";
		Dzero(xi_x[0], 4);
	}
    else
    {
        idet = 1.0/det;
        xi_x[0][0] =  x_xi[1][1] * idet;
        xi_x[1][1] =  x_xi[0][0] * idet;
        xi_x[0][1] = -x_xi[0][1] * idet;
        xi_x[1][0] = -x_xi[1][0] * idet;
    }
	
	
    /* finally, compute the cartesian shape function derivatives pushing forward
	 the isoparametric ones with the inverse transposed jacobian */
    for (int a=0; a<shp.size(); a++)
    {
        shpx[a].value   = shp[a].value;
        shpx[a].dxyz[0] = shp[a].dxyz[0] * xi_x[0][0] + shp[a].dxyz[1] * xi_x[1][0];
        shpx[a].dxyz[1] = shp[a].dxyz[0] * xi_x[0][1] + shp[a].dxyz[1] * xi_x[1][1];
    }
}





// constructor for the shapefunctions blocks. 
// the blocks are always associated with an element, it is not just a collection of shapefun's
// theOwner is just a ressource manager that takes care of eliminating the shapefunctions contained in the object
curveShapefuns :: curveShapefuns(const element& e)
: FEshapefunbuilder(e)
{ 
}





/* this is as above, for rods or similar */
void curveShapefuns :: evaluate(const quadpoint& gp, double& jacobian)
{ 
    eval1d(gp, jacobian, NOSECONDD, PUSH);
}


void curveShapefuns :: evaluateRef(const quadpoint& gp)
{
    double tmp;
    
	eval1d(gp, tmp, NOSECONDD, NOPUSH);
}


