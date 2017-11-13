 /*
 * Original author: jesus.hernandez.mangas@tel.uva.es
 *
 * Copyright 2014 University of Valladolid, Valladolid, Spain.
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
 *
 * Adapted to MMonCa: I. Martin-Bragado. IMDEA Materials Institute.
 *
 */ 
///////////////////////////////////////////////////////////////////////////
//
// VARHST.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_VARHST_H_
#define _INCLUDE_VARHST_H_

#define t_double float
//#define t_double double

#include <iostream>
#include <assert.h>
#include "vector.h"
#include "defs.h"
 
#ifdef MODEL2
#define N_BINS   100
#define N_SECTOR 21
#else
#define N_BINS   100
#define N_SECTOR 21
#endif
namespace BCA {
class VarHst
{
 public:
  t_double M[N_BINS];
  t_double mPrimero, mUltimo, mAncho;
 //------------------------------------------------------------------------
 void Set( t_double mP=0.0, t_double mU=10.0 )
  {
   mPrimero = mP;
   mUltimo  = mU;
   mAncho   = (mU-mP)/N_BINS;
   assert(mAncho>0);
  }
 //------------------------------------------------------------------------
 void Clear(t_double Val=0.0)
  {
   for(int i=0;i<N_BINS;i++) M[i]=Val;
  }
 //------------------------------------------------------------------------
 t_double Max()
  {
   t_double max = 0.0;
   for(int i=0;i<N_BINS;i++) if(max<M[i]) max=M[i];
   return(max);
  }
 //------------------------------------------------------------------------
 t_double Area()
  {
   t_double Area = 0.0;
   for(int i=0;i<N_BINS;i++) Area = Area + M[i];
   Area = Area * mAncho;
   return(Area);
  }
  
 //------------------------------------------------------------------------
 int GetIndex( t_double Value )
  {
    int Index = (int) (Value/mAncho);
    if(Value<0) return -1;
    if(Index>N_BINS-1)
     {
      ReScale(mPrimero,(Index+2)*mAncho);
      Index = (int) (Value/mAncho);
     }
    return Index;
  };
 //------------------------------------------------------------------------
 int ReScale( t_double nP, t_double nU )
  {
    t_double N[N_BINS];
    int j, i=0, I=0, Plus=1;

    t_double nAncho = (nU-nP)/N_BINS;
    assert(nAncho>0);
    if(nAncho<=mAncho) return(0);
    t_double mnDelta  = nP - mPrimero;
    for(j=0;j<N_BINS;j++) N[j]=0.0;
    do
    {
      if( mnDelta + nAncho*(I+1) < mAncho*i )
        while( mnDelta + nAncho*(I+1) < mAncho*i ) I ++;
      else
       {
        N[I] = M[i] * ( mAncho*(i+1) - mnDelta - nAncho*(I) );
        i ++;   if( i > N_BINS-1 ) Plus = 0;
       };
      while( (mAncho*(i+1) < mnDelta + nAncho*(I+1)) && Plus )
       {
          if( Plus ) N[I] = N[I] + mAncho * M[i];
          i ++;   if( i > N_BINS-1 ) Plus = 0;
       };
      if( Plus )
         N[I] = N[I] + M[i] * ( mnDelta + nAncho*(I+1) - mAncho*i );
      N[I] = N[I] / nAncho;
      I ++;
    } while( (I < N_BINS) && (Plus) );
    for(i=0;i<N_BINS;i++) M[i]=N[i];
    mPrimero=nP; mUltimo=nU; mAncho=nAncho;
    return(1);
  };
 //------------------------------------------------------------------------
  void PrintAll()
  {
    int i,j;
    std::cout << "("<< mPrimero << "," << mUltimo << "), " << mAncho << std::endl;
    for( i = 0; i < N_BINS ; i++ )
     {
    	std::cout.width(10);
    	std::cout << i*mAncho << " ";
    	std::cout.width(10);
    	std::cout << M[i] << " ";
      for ( j = 1; j < M[ i ]; j++ ) std::cout << "*";
      std::cout << std::endl;
     }
    std::cout << std::endl;
  };
};

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// VarHst3D ===============================================================
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
class VarHst3D
{
 public:
#ifdef INTEL_C
  t_double (*M)[N_SECTOR][N_SECTOR];  
#else
  t_double M[N_BINS][N_SECTOR][N_SECTOR];
#endif
  t_double mPrimero, mUltimo, mAncho;
  t_double mExtremoY, mExtremoZ, mAnchoY, mAnchoZ;

#ifdef INTEL_C
  VarHst3D()
   {
     M=new t_double [N_BINS][N_SECTOR][N_SECTOR]; 
   }
  ~VarHst3D()
   {
     delete [] M;
   }
#endif

 //------------------------------------------------------------------------
 void Set3D( t_double mP=0.0, t_double mU=1000.0, t_double mY=10.0, t_double mZ=10.0)
  {
   mPrimero = mP;
   mUltimo  = mU;
   mAncho   = (mU-mP)/N_BINS;
   assert(mAncho>0);
   mExtremoY = mY;
   mExtremoZ = mZ;
   mAnchoY   = (mY+mY)/N_SECTOR;
   mAnchoZ   = (mZ+mZ)/N_SECTOR;
   assert(mY>0 && mZ>0);
  }
 //------------------------------------------------------------------------
 void Clear3D(t_double Val=0.0)
  {
   int i,j,k;
   for(i=0;i<N_BINS;i++)
    for(j=0;j<N_SECTOR;j++)
     for(k=0;k<N_SECTOR;k++)
      M[i][j][k]=Val;
  }
 //------------------------------------------------------------------------
 t_double Max()
  {
   int i,j,k;
   t_double max;
   
   for(i=0;i<N_BINS;i++)
    for(j=0;j<N_SECTOR;j++)
     for(k=0;k<N_SECTOR;k++)
      if(max<M[i][j][k]) max=M[i][j][k];
   return(max);
  }
 //------------------------------------------------------------------------
 Vector GetIndex( Vector Value )
  {
    Vector Index;
    Index.X = (int) (Value.X/mAncho);
    if(Value.X<0) Index.X=-1;
    if(Index.X>N_BINS-1)
     {
      ReScaleX(mPrimero,(Index.X+1)*mAncho);
      Index.X = N_BINS-1;
     }
    Index.Y = (int) (Value.Y/mAnchoY + N_SECTOR/2);  
    if(Index.Y>N_SECTOR-1) Index.Y = N_SECTOR-1;
    if(Index.Y<0         ) Index.Y = 0;         
    Index.Z = (int) (Value.Z/mAnchoZ + N_SECTOR/2);
    if(Index.Z>N_SECTOR-1) Index.Z = N_SECTOR-1;
    if(Index.Z<0         ) Index.Z = 0;
    return Index;
  };
 //------------------------------------------------------------------------
 int ReScaleX( t_double nP, t_double nU )
  {
#ifdef INTEL_C
 #ifdef WINDOWS
    t_double (*N)[N_SECTOR][N_SECTOR]=new t_double [N_BINS][N_SECTOR][N_SECTOR]; 
 #else    
    t_double N[N_BINS][N_SECTOR][N_SECTOR];
 #endif
#else
    t_double N[N_BINS][N_SECTOR][N_SECTOR];
#endif
    int j, i, I, Plus, a, b;
    t_double nAncho = (nU-nP)/N_BINS;
    assert(nAncho>0);
    if(nAncho<mAncho) return(0);
    t_double mnDelta  = nP - mPrimero;
    for(a=0;a<N_SECTOR;a++)
     for(b=0;b<N_SECTOR;b++)
      for(j=0;j<N_BINS;j++)
       N[j][a][b]=0.0;
    for(a=0;a<N_SECTOR;a++)
     for(b=0;b<N_SECTOR;b++)
      {
       I=0; i=0; Plus=1;
       do
        {
         if( mnDelta + nAncho*(I+1) < mAncho*i )
          while( mnDelta + nAncho*(I+1) < mAncho*i ) I ++;
         else
         {
          N[I][a][b] = M[i][a][b] * ( mAncho*(i+1) - mnDelta - nAncho*(I) );
          i ++;   if( i > N_BINS-1 ) Plus = 0;
         };
         while( (mAncho*(i+1) < mnDelta + nAncho*(I+1)) && Plus )
         {
          if( Plus ) N[I][a][b] = N[I][a][b] + mAncho * M[i][a][b];
          i ++;   if( i > N_BINS-1 ) Plus = 0;
         };
         if( Plus )
           N[I][a][b] += M[i][a][b] * ( mnDelta + nAncho*(I+1) - mAncho*i );
         N[I][a][b] = N[I][a][b] / nAncho;
         I ++;
        } while( (I < N_BINS) && (Plus) );
       }
    for(a=0;a<N_SECTOR;a++)
     for(b=0;b<N_SECTOR;b++)
      for(j=0;j<N_BINS;j++)
       M[j][a][b]=N[j][a][b];
    mPrimero=nP; mUltimo=nU; mAncho=nAncho;
#ifdef INTEL_C
//    delete [] N;
#endif
    return(1);
  };
 //------------------------------------------------------------------------

  void PrintAll()
  {
    int i,j;
    std::cout << "("<< mPrimero << "," << mUltimo << "), " << mAncho << std::endl;
    std::cout << "("<< -mExtremoY << "," << mExtremoY << "), " << mAnchoY << std::endl;
    std::cout << "("<< -mExtremoZ << "," << mExtremoZ << "), " << mAnchoZ << std::endl;
    for( i = 0; i < N_BINS ; i++ )
     {
    	std::cout.width(10);
    	std::cout << i*mAncho << " ";
    	std::cout.width(10);
    	std::cout << M[i][0][0] << " ";
      for ( j = 1; j < M[ i ][0][0]; j++ ) std::cout << "*";
      std::cout << std::endl;
     }
    std::cout << std::endl;
  };

};
}
#endif //_INCLUDE_VARHST_H_
