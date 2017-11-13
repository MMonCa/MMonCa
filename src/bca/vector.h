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
// VECTOR.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_VECTOR_H_
#define _INCLUDE_VECTOR_H_

#include <assert.h>
#include <iostream>
#include <math.h>

#ifndef M_PI
#define M_PI 3.141592654
#endif

#define WIDTH 9

namespace BCA {
class Vector
 {
  public:
   double X, Y, Z;

   double GetMax();
   void Assign( double XC, double YC, double ZC );
   Vector        operator()( double XC, double YC, double ZC );
   Vector        operator+ (Vector SecondOperand );
   Vector        operator- ( Vector SecondOperand );
   Vector        operator- ();
   Vector        operator^ ( Vector &SecondOperand );
   double          operator* ( Vector SecondOperand );
   int           operator==( Vector SecondOperand );
   int           operator!=( Vector SecondOperand );
   Vector        operator* ( double K );

   friend Vector operator* ( double K, Vector V )
   {
    Vector R;

    R.X = V.X*K;
    R.Y = V.Y*K;
    R.Z = V.Z*K;
    return R;
   }

   Vector operator /( double K );
   Vector operator%( Vector LP );
   Vector Integer( double Divisor );
   Vector Integer2( Vector Divisor, Vector LP );
   operator float();
   operator double();
   operator long double();

   friend Vector Unitary( Vector &V )
   {
    Vector R;

    double M = sqrt( V.X*V.X + V.Y*V.Y + V.Z*V.Z ); // Modulus
    if( M!=0 )
     {
       M = 1 / M;           // Fastest than division
       R.X = V.X*M;
       R.Y = V.Y*M;
       R.Z = V.Z*M;
     }
    else
     {
       R = V;
     }
    return R;
   }

   friend double MixedProduct( Vector &A, Vector &B, Vector &C )
   {
      //    return (A*(B^C));
      return ( A.X*(B.Y*C.Z-B.Z*C.Y)+
           A.Y*(B.Z*C.X-B.X*C.Z)+
           A.Z*(B.X*C.Y-B.Y*C.X) ); // 60% Fastest
   }

   friend double Angle( Vector &A, Vector &B )
   {
      // return ( acos(A*B / ( (double )A*(double )B )));
      return ( acos(       ( A.X*B.X + A.Y*B.Y + A.Z*B.Z) // 30% Fastest
             /sqrt( (A.X*A.X + A.Y*A.Y + A.Z*A.Z)
               *(B.X*B.X + B.Y*B.Y + B.Z*B.Z))
           )
         );
   }
   
   double Deep( Vector Direction );

   // Change the coordinates to another reference given by NewX, NewY & NewZ
   friend Vector ChangeToReference( Vector &NewX, Vector &NewY,
                    Vector &NewZ, Vector Old  )
   {
    Vector R;
    // return NewX * Old.X + NewY * Old.Y + NewZ * Old.Z;
    // 50% Fastest
    R.X = NewX.X*Old.X + NewY.X*Old.Y +NewZ.X*Old.Z;
    R.Y = NewX.Y*Old.X + NewY.Y*Old.Y +NewZ.Y*Old.Z;
    R.Z = NewX.Z*Old.X + NewY.Z*Old.Y +NewZ.Z*Old.Z;
    return R;
   }

   // Discretize precision of Vector
   friend Vector Precision( Vector V, double Delta, double Delta2 )
   {
    if( fabs(V.X) < Delta ) V.X = Delta2;
    if( fabs(V.Y) < Delta ) V.Y = Delta2;
    if( fabs(V.Z) < Delta ) V.Z = Delta2;
    return V;
   }

   // Round to given precision
   friend Vector Round( Vector &V, double Delta )
   {
    Vector R;

    R.X = (long double) V.X;
    R.Y = (long double) V.Y;
    R.Z = (long double) V.Z;
    assert( Delta != 0 );
    R.X = Delta*( (long int)(R.X / Delta) );
    R.Y = Delta*( (long int)(R.Y / Delta) );
    R.Z = Delta*( (long int)(R.Z / Delta) );
    return R;
   }

   friend std::ostream &operator<<(std::ostream &Out, Vector V );
   friend std::istream &operator>>(std::istream &In,  Vector &V );
 };

// Get Maximum coordinate
inline double Vector::GetMax()
{
  double m;

  m = (X > Y) ? X : Y;
  return( (m>Z)?m:Z );
}

// Assign vector
inline void Vector::Assign( double XC, double YC, double ZC )
{
  X = XC;
  Y = YC;
  Z = ZC;
}

// Assign vector (natural form) VECTOR(X,Y,Z) are VECTOR=(X,Y,Z)
inline Vector Vector::operator()( double XC, double YC, double ZC )
{
  X = XC;
  Y = YC;
  Z = ZC;
  return *this;
}

// Vectorial add
inline Vector Vector::operator+(Vector SecondOperand )
{
  Vector R;

  R.X = X + SecondOperand.X;
  R.Y = Y + SecondOperand.Y;
  R.Z = Z + SecondOperand.Z;
  return R;
}


// Vectorial substract
inline Vector Vector::operator-( Vector SecondOperand )
{
  Vector R;

  R.X = X - SecondOperand.X;
  R.Y = Y - SecondOperand.Y;
  R.Z = Z - SecondOperand.Z;
  return R;
}

// Unary negation
inline Vector Vector::operator-()
{
  Vector R;

  R.X = -X;
  R.Y = -Y;
  R.Z = -Z;
  return R;
}

// Vectorial product
inline Vector Vector::operator^( Vector &SecondOperand )
{
  Vector R;

  R.X = Y*SecondOperand.Z - Z*SecondOperand.Y;
  R.Y = Z*SecondOperand.X - X*SecondOperand.Z;
  R.Z = X*SecondOperand.Y - Y*SecondOperand.X;
  return R;
}

// Scalar product
inline double Vector::operator*( Vector SecondOperand )
{
  return (X*SecondOperand.X + Y*SecondOperand.Y + Z*SecondOperand.Z);
}

// Vector comparison ==
inline int Vector::operator==( Vector SecondOperand )
{
  if ( (X == SecondOperand.X) &&
       (Y == SecondOperand.Y) &&
       (Z == SecondOperand.Z) ) return 1;
  else return 0;
}

// Vector comparison !=
inline int Vector::operator!=( Vector SecondOperand )
{
  if ( (X != SecondOperand.X) ||
       (Y != SecondOperand.Y) ||
       (Z != SecondOperand.Z) ) return 1;
  else return 0;
}

// Product by a double constant:
// Vector * Constant
inline Vector Vector::operator*( double K )
{
  Vector R;

  R.X = X*K;
  R.Y = Y*K;
  R.Z = Z*K;
  return R;
}

// Division by a double constant
inline Vector Vector::operator /( double K )
{
  Vector R;

  assert( K!=0 ); // Debugging

  double M = 1/K;  // 10% Fastest
  R.X = X*M;
  R.Y = Y*M;
  R.Z = Z*M;
  return R;
}

// Vector % Vector   . Normalize components
inline Vector Vector::operator%( Vector LP )
{
  Vector R;

  R.X = fmod( X , LP.X );
  R.Y = fmod( Y , LP.Y );
  R.Z = fmod( Z , LP.Z );
  return R;
}

// Integer part of a Vector                      // OK
inline Vector Vector::Integer( double Divisor )
{
  assert( Divisor != 0 ); // Debugging

  X = ((long int) (X / Divisor));
  Y = ((long int) (Y / Divisor));
  Z = ((long int) (Z / Divisor));
  return *this;          // Not only assign but automodified
}

// Integer part of a Vector                      // OK
inline Vector Vector::Integer2( Vector Divisor, Vector LP )
{
  X = LP.X * ((long int) (X / Divisor.X));
  Y = LP.Y * ((long int) (Y / Divisor.Y));
  Z = LP.Z * ((long int) (Z / Divisor.Z));
  return *this;
}

// Modulus of a vector
inline Vector::operator float()       { return (float)       sqrt( X*X + Y*Y + Z*Z ); }
inline Vector::operator double()      { return (double)      sqrt( X*X + Y*Y + Z*Z ); }
inline Vector::operator long double() { return (long double) sqrt( X*X + Y*Y + Z*Z ); }

// Deep function in any direction (SIGNED VALUE)
inline double Vector::Deep( Vector Direction )
{
  Vector R;

  double M = sqrt( Direction.X*Direction.X
               + Direction.Y*Direction.Y
               + Direction.Z*Direction.Z ); // Modulus
  assert( M!=0 ); // Debugging
  M = 1 / M;
  R.X = Direction.X*M*X;
  R.Y = Direction.Y*M*Y;
  R.Z = Direction.Z*M*Z;
  return R;
}

//----------------------------

void Initialize();

class Orientacion
{
 public:
   Vector EjeX, EjeY, EjeZ;

   void Orienta(Vector V, Vector Flat, double Tha_Cut, double Phi_Cut, int Show=0);
   double ProyectaX( Vector V );
   double ProyectaY( Vector V );
   double ProyectaZ( Vector V );
   Vector Proyecta( Vector V );
};
}
#endif  // _INCLUDE_VECTOR_H_
