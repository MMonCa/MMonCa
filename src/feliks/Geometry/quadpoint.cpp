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
 * quadpoint.cpp
 * 
 * implement quadrature rules 
 *
 * i. romero, february 2002
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "quadpoint.h"
#include "Io/message.h"



static void     QuadraturePoints1d(const size_t nnode, const size_t ngp, quadpoint *q);
static void     QuadraturePoints2d(const size_t nnode, const size_t ngp, quadpoint *q);
static void     QuadraturePoints3d(const size_t nnode, const size_t ngp, quadpoint *q);



/* this fills up a list of quadrature points, that accomodates to the number
 of nodes and dimension of the reference cell. The gausslist received
 better have enough space for all the nodes!
 */
void fullQuadraturePoints(const size_t nnode, const geometryT geo, size_t *ngp, quadpoint *gausslist)
{
    
    *ngp = numberOfQuadraturePoints(nnode, geo);
    
    switch (geo)
    {
        case CURVE:     QuadraturePoints1d(nnode, *ngp, gausslist); break;
        case SURFACE:   QuadraturePoints2d(nnode, *ngp, gausslist); break;
        case VOLUME:    QuadraturePoints3d(nnode, *ngp, gausslist); break;
        default:
            ErrorMessage("Quadrature rule not implemented for geometry %d", geo);
    }
}





void  quadpoint :: print() const
{
    Message("\n\n Reference element coordinates (%f,%f,%f)", coor[0], coor[1], coor[2]);
    Message("\n Weight = %f", weight);
}






/* this fills up a list of quadrature points. The gausslist received
 better have enough space for all the nodes!
 */
void quadraturePoints(const size_t nnode, geometryT geo , const size_t ngp, quadpoint *gausslist)
{
    switch (geo)
    {
        case 1: QuadraturePoints1d(nnode, ngp, gausslist); break;
        case 2: QuadraturePoints2d(nnode, ngp, gausslist); break;
        case 3: QuadraturePoints3d(nnode, ngp, gausslist); break;
        default:
            ErrorMessage("Quadrature rule not implemented");
    }
}




void QuadraturePoints1d(const size_t nnode, const size_t ngp, quadpoint *q)
{
    static const double s3   = 0.57735026918962576450914878050195746; // 1/sqrt(3)
    static const double s35  = 0.77459666924148337703585307995647992; // sqrt(3/5)
    static const double s120 = 10.9544511501033222691393956560160427; // sqrt(120)
    int    a;
    
    switch (ngp)
    {
        case(1):
            q[0].coor[0] = 0.0;
            q[0].weight  = 2.0;
            break;
            
        case(2):
            q[0].coor[0] = -s3;
            q[1].coor[0] =  s3;
            
            q[0].weight  = 1.0;
            q[1].weight  = 1.0;
            break;
            
        case(3):
            q[0].coor[0] = -s35;
            q[1].coor[0] =  0;
            q[2].coor[0] =  s35;
            
            q[0].weight  =  5.0/9.0;
            q[1].weight  =  8.0/9.0;
            q[2].weight  =  5.0/9.0;
            break;
            
		case(4):
			q[0].coor[0] = -sqrt(3.0/7.9 + s120/35.0);
			q[1].coor[0] = -sqrt(3.0/7.9 - s120/35.0);
			q[2].coor[0] = -q[1].coor[0];
			q[3].coor[0] = -q[0].coor[0];
			
			q[0].weight = 0.5 - 5.0/(3.0*s120);
			q[1].weight = 0.5 + 5.0/(3.0*s120);
			q[2].weight = q[1].weight;
			q[3].weight = q[0].weight;
			break;
			
            
        default:
            Message("\n Quadrature rule for 1d element not implemented");
    }
    
    /* reset the second and third coordinate to zero */
    for (a=0; a<nnode; a++)
    {
        q[a].coor[1] = 0.0;
        q[a].coor[2] = 0.0;
    }
    
}




/* quadrature rule for 2d elements. See Bathe pg. 467 
 */
void  QuadraturePoints2d(const size_t nnode, const size_t ngp , quadpoint *gausslist)
{
    static const double s3   = 0.57735026918962576450914878050195746; // 1/sqrt(3)
    static const double s35  = 0.77459666924148337703585307995647992; // sqrt(3/5)

    int        a, i, j;
    double     s48;
    static double xi2sgn[]={+1.0 , -1.0 , -1.0 , +1.0};
    static double et2sgn[]={+1.0 , +1.0 , -1.0 , -1.0};
    static double wg3[]   ={0.5555555556 ,  0.8888888889 , 0.5555555556};
    static double cst[]   ={0.1666666666667, 0.6666666666667};
    static double tri2r[] ={0.1012865073235, 0.7974269853531, 0.1012865073235, 0.4701420641051,
        0.4701420641051, 0.0597158717898, 0.3333333333333};
    static double tri2s[] ={0.1012865073235, 0.1012865073235, 0.7974269853531, 0.0597158717898,
        0.4701420641051, 0.4701420641051, 0.3333333333333};
    static double tri2w[] ={0.1259391805448, 0.1259391805448, 0.1259391805448, 0.1323941527885,
        0.1323941527885, 0.1323941527885, 0.2250000000000};
    double cstw;
    double xi4[4], wg4[4];
    
    
    switch (nnode*100 + ngp)
    {
            
        // quadrature rules for the 3node triangle
        case(301):
            gausslist[0].coor[0] = 1.0/3.0;
            gausslist[0].coor[1] = 1.0/3.0;
            gausslist[0].coor[2] = 0.0;
            gausslist[0].weight  = 0.5;
            break;
            
        case(303):
        case(603):
            cstw = 1.0/6.0;
            
            gausslist[0].coor[0] = cst[0];
            gausslist[0].coor[1] = cst[0];
            gausslist[0].coor[2] = 0.0;
            gausslist[0].weight  = cstw;
            
            gausslist[1].coor[0] = cst[1];
            gausslist[1].coor[1] = cst[0];
            gausslist[1].coor[2] = 0.0;
            gausslist[1].weight  = cstw;
            
            gausslist[2].coor[0] = cst[0];
            gausslist[2].coor[1] = cst[1];
            gausslist[2].coor[2] = 0.0;
            gausslist[2].weight  = cstw;
            break;
            
            /* quadrature rule for the 4 node quad */
        case(401):
            gausslist[0].coor[0] = 0.0;
            gausslist[0].coor[1] = 0.0;
            gausslist[0].weight  = 4.0;
            break;
            
        case(404):
            for (a=0; a<4; a++)
            {
                gausslist[a].coor[0] = s3*xi2sgn[a];
                gausslist[a].coor[1] = s3*et2sgn[a];
                gausslist[a].coor[2] = 0.0;
                gausslist[a].weight  = 1.0;
            }
            break;
            
        case(409):
            a   = 0;
            for(i=0; i<3; i++)
            {
                for(j=0; j<3; j++)
                {
                    gausslist[a].coor[0] = (i-1.0)*s35;
                    gausslist[a].coor[1] = (j-1.0)*s35;
                    gausslist[a].coor[2] = 0.0;
                    gausslist[a].weight  = wg3[i]*wg3[j];
                    a++;
                }
            }
            break;
            
        case(416):
            s48    = sqrt(4.8);
            xi4[0] = sqrt((3.0+s48)/7.0);
            xi4[1] = sqrt((3.0-s48)/7.0);
            xi4[2] = -xi4[1];
            xi4[3] = -xi4[0];
            s48    = 1.0/(3.0*s48);
            wg4[0] = 0.5 - s48;
            wg4[1] = 0.5 + s48;
            wg4[2] = 0.5 + s48;
            wg4[3] = 0.5 - s48;
            a = 0;
            for(i=0; i<4; i++)
                for(j=0; j<4; j++)
                {
                    gausslist[a].coor[0] = xi4[j];
                    gausslist[a].coor[1] = xi4[i];
                    gausslist[a].coor[2] = 0.0;
                    gausslist[a].weight  = wg4[i]*wg4[j];
                    a++;
                }
            break;
            
            
            /* quadratic triangle */
        case(607):
            for(a=0; a<7; a++)
            {
                gausslist[a].coor[0] = tri2r[a];
                gausslist[a].coor[1] = tri2s[a];
                gausslist[a].coor[2] = 0.0;
                gausslist[a].weight  = 0.5*tri2w[a];
            }
            break;
            
        default:
            Message("Warning: 2d quadrature rule for %d nodes, %d gausspoints, does not exist", nnode, ngp);
    }
}




/* quadrature rule for 3d elements. This allocates the memory for
 * the list of quadrature points, so it must be deallocated after use.
 *
 */
void QuadraturePoints3d(const size_t nnode, const size_t ngp, quadpoint *gausslist)
{
    static const double s3   = 0.57735026918962576450914878050195746; // 1/sqrt(3)

    int        a;
    static double xisgn[]={+1.0, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0, +1.0};
    static double etsgn[]={+1.0, +1.0, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0};
    static double zesgn[]={-1.0, -1.0, -1.0, -1.0, +1.0, +1.0, +1.0, +1.0};
    static double tet4[] ={0.5854101966249685, 0.138196601125015, 0.1381966011250105, 0.1381966011250105};
    static double tet10[]={0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.5e0};
    static int    eps4[] ={0, 1, 2, 3, 0, 1};
    
    
    switch (100*nnode + ngp)
    {
            /* quadrature rules for the tetrahedra. 
             see "Symmetric Gaussian Quadrature formulae for tetrahedron regions",
             Y. Jinyun, CMAME 43 (1984) 349-353 */
            
            /* one point quadrature, error = O(h^2) */
        case(401):
            gausslist[0].coor[0] = 0.25;
            gausslist[0].coor[1] = 0.25;
            gausslist[0].coor[2] = 0.25;
            gausslist[0].weight  = 1.0;
            break;
            
            /* 4 point quadrature, error = O(h^3) */
        case(404):
            for(a=0; a<4; a++)
            {
                gausslist[a].coor[0] = tet4[eps4[a  ]];
                gausslist[a].coor[1] = tet4[eps4[a+1]];
                gausslist[a].coor[2] = tet4[eps4[a+2]];
                gausslist[a].weight  = 0.25;
            }
            break;
            
            /* five point quadrature, error = O(h^4) */
        case(405):
        case(1005):
            for(a=0; a<4; a++)
            {
                gausslist[a].coor[0] = tet10[eps4[a  ]];
                gausslist[a].coor[1] = tet10[eps4[a+1]];
                gausslist[a].coor[2] = tet10[eps4[a+2]];
                gausslist[a].weight  = 0.45e0;
            }
            gausslist[4].coor[0] = 0.25;
            gausslist[4].coor[1] = 0.25;
            gausslist[4].coor[2] = 0.25;
            gausslist[4].weight  = -0.8e0;
            break;
            
            /* quadrature rules for the hex */
        case(801):
            gausslist[0].coor[0] = 0.0;
            gausslist[0].coor[1] = 0.0;
            gausslist[0].coor[2] = 0.0;
            gausslist[0].weight  = 8.0;
            break;
            
        case(808):
            for (a=0; a<8; a++)
            {
                gausslist[a].coor[0] =  xisgn[a]*s3;
                gausslist[a].coor[1] =  etsgn[a]*s3;
                gausslist[a].coor[2] =  zesgn[a]*s3;
                gausslist[a].weight  =  1.0;
            }
            break;
            
        default:
            WarningMessage("Quadrature rule for 3d elements with this number of nodes not implemented.");
    }
}




/* this fills up a list of quadrature points, that accomodates to the number
 of nodes and dimension of the reference cell. The gausslist received
 better have enough space for all the nodes!
 */
void SPRQuadraturePoints(const size_t nnode, const geometryT geo, const bool corner, size_t *ngp, quadpoint *gausslist)
{
	if (corner)
		*ngp = numberOfQuadraturePoints(nnode, geo);
	else
		*ngp = numberOfSPRQuadraturePoints(nnode, geo);
	
	switch (geo)
    {
		case CURVE:		QuadraturePoints1d(nnode, *ngp, gausslist); break;
		case SURFACE:	QuadraturePoints2d(nnode, *ngp, gausslist); break;
		case VOLUME:
        case COHESIVE:
            QuadraturePoints3d(nnode, *ngp, gausslist); break;
		default:
			ErrorMessage("SPR Quadrature rule not implemented for geometry %d, nnode = %d, ngp = %d", geo, nnode, *ngp);
    }	
}






/* return the default number of (full) quadrature points for each element
 * type.
 */
size_t numberOfQuadraturePoints(const size_t nnode, const geometryT geo)
{
    size_t ngp=0;
    
    if (geo == CURVE)
    {
        switch (nnode)
        {
            case 2:  ngp = 2; break;
            case 3:  ngp = 3; break;
            case 4:  ngp = 4; break;
            case 5:  ngp = 5; break; /* any 1d element */
        }
    }
    
    else if (geo == SURFACE)
    {
        switch (nnode) 
        {
            case 3:  ngp = 1; break; /* cst */
            case 4:  ngp = 4; break; /* quad */
            case 6:  ngp = 3; break; /* quadratic triangle */
            case 8:  ngp = 9; break; /* serendipity */
            case 9:  ngp = 9; break; /* q2 quad */
        }   
    }
    
    else if (geo == VOLUME)
    {
        switch (nnode)
        {
            case 4:   ngp = 4; break; /* tet*/
            case 8:   ngp = 8; break; /* brick */
            case 10:  ngp = 5; break; /* quadratic tet */
        }
    }
    
    if (ngp == 0)
        WarningMessage("\nQuadrature rule for %d nodes in %d geometry type", nnode, geo);
    
    return ngp;
}




/* return the default number of superconvergent quadrature points for each element
 * type.
 */
size_t numberOfSPRQuadraturePoints(const size_t nnode, const geometryT geo)
{
    size_t ngp=0;
    
    if (geo == CURVE)
    {
        switch (nnode)
        {
            case 2:  ngp = 1; break;
            case 3:  ngp = 2; break;
            case 4:  ngp = 3; break;
            case 5:  ngp = 4; break; /* any 1d element */
        }
    }
    
    else if (geo == SURFACE)
    {
        switch (nnode) 
        {
            case 3:  ngp = 1; break; /* cst */
            case 4:  ngp = 1; break; /* quad */
            case 6:  ngp = 3; break; /* quadratic triangle */
            case 8:  ngp = 9; break; /* serendipity */
            case 9:  ngp = 9; break; /* q2 quad */
        }   
    }
    
    else if (geo == VOLUME)
    {
        switch (nnode)
        {
            case 4:   ngp = 4; break; /* tet*/
            case 8:   ngp = 8; break; /* brick */
            case 10:  ngp = 5; break; /* quadratic tet */
        }
    }
    
    if (ngp == 0)
        WarningMessage("SPR Quadrature rule for %d nodes in %d geometry type", nnode, geo);
    
	return ngp;
}





