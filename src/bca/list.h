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
// LIST.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_LIST_H_
#define _INCLUDE_LIST_H_

#ifdef UNIX
 #ifdef INTEL_C
  #define  ANSI_TEMPLATES 
  #define TEMPLATE_INSTANCES
 #else
  #ifdef GCC3
   #define ANSI_TEMPLATES
   #define TEMPLATE_INSTANCES
  #else
   #define ANSI_TEMPLATES
   #define TEMPLATE_INSTANCES
  #endif 
 #endif
#else
 #ifdef INTEL_C
  #undef  ANSI_TEMPLATES
  #define TEMPLATE_INSTANCES
 #else
  #define ANSI_TEMPLATES
  #define TEMPLATE_INSTANCES
 #endif 
#endif

#include "vector.h"

#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

// DEFINES
//#define WHOLE_LIST  // Activate all LIST capabilities (needs more memory)

#define integer int

#ifndef NDEBUG
#define bGarbage (char)0xCC
#endif

#undef  DOUBLY

// Define the template for the dynamic link to an object //////////////////
namespace BCA {
template <class GENERIC> class LinkObject
{
 public:
   LinkObject<GENERIC> *Next;       // Link to the NEXT object
   #ifdef DOUBLY
   LinkObject<GENERIC> *Previous;   // Link to the PREVIOUS object
   #endif
   GENERIC      Data;               // Content of the object

   LinkObject()                     // Constructor
   {
     Next     = NULL;   // Initialize to null pointers
     #ifdef DOUBLY
     Previous = NULL;
     #endif
   };

   #ifndef NDEBUG
   ~LinkObject();                   // Destructor
   #endif
};


template <class GENERIC> class List
{
 public:
    integer N;
 protected:
    LinkObject<GENERIC>* Begin; // First element
    LinkObject<GENERIC>* End;   // Last element
    LinkObject<GENERIC>* Last;  // Last element accessed

    void  QuickSort_(int (*MatchFunction)(GENERIC P1, GENERIC P2 ),
                   LinkObject<GENERIC> **Head,
                   LinkObject<GENERIC>  *End
                 );

 public:
    // Constructor (inlined)
    List()
     {
        Begin = NULL;
        End   = NULL;
        Last  = NULL;
        N     =    0;
     };

    #ifndef NDEBUG
    ~List();    // Destructor
    #endif

    // Get the number of elements
    integer NumberOfElements() { return N;};

    // Copy method
    integer   Set       ( List<GENERIC> Old );

    // Insertion methods
    integer   Add       ( GENERIC Element );
    integer   InsertAt  ( integer Position, GENERIC &Element );

    // Deleteing methods
    integer   DeleteAt  ( integer Position );
    integer   Delete    ( LinkObject<GENERIC> * TempDel );
    integer   DeleteLast();
    integer   DeleteAll ();

    // Exchange method
    integer   Exchange(LinkObject<GENERIC> *P1, LinkObject<GENERIC> *P2);

    // Quicksort function
    void  QuickSort ( int (*MatchFunction)(GENERIC P1, GENERIC P2 ) );

    // Selection Sort
    integer   Sort( double (*fMatch)          // Match function
                ( GENERIC &E, Vector &D ),
            Vector Distance ,   // Argument function
            int Order           //  0 Ascendent
              );                //  1 Descendent

    // Serching methods
    LinkObject<GENERIC>* PointerToFirstData()
    {
        assert ( Begin!=NULL ); // For debugging
        Last  = Begin->Next;
        return  Begin;
    };

    LinkObject<GENERIC>* PointerToNextData()
    {
        LinkObject<GENERIC>* T;

        T = Last;
        if( Last != NULL ) Last  = Last->Next;
        return  T;
    };


    LinkObject<GENERIC>* Search( integer Position );

    // Print method
    integer Save( char* FileName, long int NN,
                 double (*fMatch)( GENERIC &A, Vector V ) );

    // Search the "Nearest" given by fMatch function (typically ByDistance)
    GENERIC Nearest( double (*fMatch)(GENERIC &D, Vector &V), Vector R );

    // Output streams overloaded : prototype
#ifdef ANSI_TEMPLATES
    template <typename T>
    friend std::ostream &operator<<(std::ostream &Out, List<T> L );
#else
    friend ostream &operator<< ( ostream &Out, List<GENERIC> L );
#endif
    // Verify list (3 Sept 1998) For debugging
    int Verify();
};

template <class GENERIC>
inline std::ostream &operator<<(std::ostream &Out, List<GENERIC> L )
{
    integer i = 0;
    LinkObject<GENERIC> * Temp;

    Temp = L.Begin;

    //Out << "# Number of elements : " << L.N << endl << endl;
    Out << std::endl;

    while ( i < L.N )
    {
     Out << Temp->Data; // Print a generic element
     Temp = Temp->Next;
     i++;
    }
    Out << std::endl;

    return Out;
}

extern Vector NullVector;
extern Vector XVector;

// Define the LinkObject class ////////////////////////////////////////////

#ifndef NDEBUG
template <class GENERIC>
LinkObject<GENERIC>::~LinkObject()  // To prevent random values and random mistakes
    {
     unsigned int i;
     char* pData;

     pData = (char*) &Data;

     for( i=0; i<sizeof(Data); i++ )
      {
       *pData = bGarbage;   // Set to a known garbage value
       pData++;
      };

     pData = (char*) &Next;

     for( i=0; i<sizeof(Next); i++ )
      {
       *pData = bGarbage;   // Set to a known garbage value
       pData++;
      };

     #ifdef DOUBLY
     pData = (char*) &Previous;

     for( i=0; i<sizeof(Previous); i++ )
      {
       *pData = bGarbage;   // Set to a known garbage value
       pData++;
      };
     #endif
    }
#endif

// Define the List class //////////////////////////////////////////////////

#ifndef NDEBUG
template <class GENERIC>
List<GENERIC>::~List()      // Destructor
     {
           unsigned int i;
           char* pData;

           pData = (char*) &Begin;

           for( i=0; i<sizeof(Begin); i++ )
        {
         *pData = bGarbage;   // Set to a known garbage value
         pData++;
        };

           pData = (char*) &End;

           for( i=0; i<sizeof(End); i++ )
        {
         *pData = bGarbage;   // Set to a known garbage value
         pData++;
        };

           pData = (char*) &Last;

           for( i=0; i<sizeof(Last); i++ )
        {
         *pData = bGarbage;   // Set to a known garbage value
         pData++;
        };

           pData = (char*) &N;

           for( i=0; i<sizeof(N); i++ )
        {
         *pData = bGarbage;   // Set to a known garbage value
         pData++;
        };

     }
#endif

// Copy method
template <class GENERIC>
integer List<GENERIC>::Set( List<GENERIC> Old )
{
 integer i;
 LinkObject<GENERIC>* BestItem;

 DeleteAll();   // Be sure that's empty

 BestItem = Old.PointerToFirstData();
 for( i=0; i<Old.NumberOfElements(); i++ )
  {
   Add( BestItem->Data);
   BestItem = Old.PointerToNextData();
  };

 return i;
}

// Add an object to the last position of the list .........................
template <class GENERIC>
integer List<GENERIC>::Add( GENERIC Element )
{
    LinkObject<GENERIC> *Object;

    Object = new LinkObject<GENERIC>;

    if (!Object) { std::cerr << "# Error: List::Add, Allocate error\n"; exit(1); }
    Object->Data = Element;

    if (N==0)           // First object
    {
     Begin         = Object;
     End           = Object;
    }
    else                // Add to end of list
    {
     End->Next     = Object;
     #ifdef DOUBLY
     Object->Previous  = End;
     #endif
     End           = Object;
    };

    N++;
    return(N);
}


// Insert an object at a position of the list .............................
template <class GENERIC>
integer List<GENERIC>::InsertAt( integer Position, GENERIC &Element )
{
    LinkObject<GENERIC> * Temp, * Object;

    if (Position<=0 || Position>N+1) return (-1);
    else
    {

     Object = new LinkObject<GENERIC>; // Allocate a new LinkObject

     if (!Object) { std::cerr << "# Error: List::InsertAt, Allocate error\n"; exit(1); }
     Object->Data = Element;

     if (N==0 && Position==1)   // At the begin and with no objects
     {
       Begin         = Object;
       End           = Object;
     }
     else if (Position==1)      // First case: at the begin
     {
       #ifdef DOUBLY
       Begin->Previous   = Object;
       Object->Previous  = NULL;
       #endif
       Object->Next      = Begin;
       Begin         = Object;
     }
     else if (Position==N+1)    // Second case: at the end
     {
       End->Next         = Object;
       #ifdef DOUBLY
       Object->Previous  = End;
       #endif
       Object->Next      = NULL;
       End           = Object;
     }
     else               // General case
     {
       #ifdef DOUBLY
       Temp = Search( Position );

       Object->Next      = Temp;
       Object->Previous  = Temp->Previous;
       Temp->Previous->Next = Object;
       Temp->Previous    = Object;
       #else
       Temp = Search( Position-1 );

       Object->Next      = Temp->Next;
       Temp->Next        = Object;
       #endif
     }

    };

    N++;
    return (N);
}

// Delete an object at the position .......................................
template <class GENERIC>
integer List<GENERIC>::DeleteAt( integer Position )
{
    #ifndef DOUBLY
    integer i;
    #endif
    LinkObject<GENERIC> * Temp, * TempDel;

    if (Position<=0 || Position>N) return (-1);
    else
    {
     #ifdef DOUBLY
     if (Position==1)      // First case
     {
      if( Begin != End )
      {
       Temp          = Begin;
       Begin         = Begin->Next;
       Begin->Previous   = NULL;
       delete Temp;
      }
      else
      {
       Temp  = Begin;
       delete Temp;
       Begin = NULL;
       End   = NULL;
      }
     }
     else if (Position==N) // Second case
     {
       Temp          = End;
       End           = End->Previous;
       End->Next         = NULL;
       delete Temp;
     }
     else           // General case
     {
       TempDel = Search( Position );
       Temp            = TempDel;
       TempDel->Previous->Next = TempDel->Next;
       TempDel->Next->Previous = TempDel->Previous;
       delete Temp;
     }
     #else
     TempDel = Begin;
     i   = Position-1;
     if (i==0)
      {
       Temp  = Begin;
       Begin = Begin->Next;
       delete Temp;
      }
     else
      {
       while( i-- )
        {
         TempDel = TempDel->Next;
        };
       Temp      = TempDel->Next;
       if( TempDel->Next->Next == NULL ) End = TempDel->Next;
       TempDel->Next = TempDel->Next->Next;

       delete Temp;
      };
     #endif

    N--;
    };
    return (N);
}

// Delete an object .......................................................
template <class GENERIC>
integer List<GENERIC>::Delete( LinkObject<GENERIC> * TempDel )
{
  LinkObject<GENERIC> * Temp;

    Temp = TempDel;

    assert( Temp != NULL ); // For debugging

       #ifdef DOUBLY
    if( TempDel->Previous == NULL )        // The first element
     {
      Begin = TempDel->Next;

      if( TempDel->Next == NULL )          // Also the last element
       {
        End = NULL;
       }
      else
       {
        TempDel->Next->Previous = NULL;    // Only the first
       };
     }
    else
     {
      TempDel->Previous->Next = TempDel->Next; // Intermediate element

      if( TempDel->Next == NULL )          // The last element
       {
        TempDel->Previous->Next = NULL;
        End             = TempDel->Previous;
       }
      else                     // Intermediate
       {
        TempDel->Next->Previous = TempDel->Previous;
       };

     }
       #else
    // Search the element that points to the current element
    TempDel = Begin;

    if( TempDel == Temp )
     {
      Begin = Temp->Next;
     }
    else
     {
      while (TempDel->Next != Temp)
       {
        TempDel = TempDel->Next;
       };
     };

    // Modify the pointers
    TempDel->Next = Temp->Next;

    if( Temp->Next == NULL ) End = TempDel;
       #endif
    delete Temp;

    N--;
    return (N);
}


// Delete the last object of the list .....................................
template <class GENERIC>
integer List<GENERIC>::DeleteLast()
{
    #ifdef DOUBLY
    LinkObject<GENERIC> * Temp;
    #else
    LinkObject<GENERIC> * TempDel;
    #endif

    if ( N == 0 ) return( -1 );
    else
     #ifdef DOUBLY
     if ( End != Begin )
      {
       Temp          = End;
       End           = End->Previous;
       End->Next         = NULL;
       delete Temp;
      }
     else
      {
       delete Begin;
       Begin = NULL;
       End   = NULL;
      }
     #else

     // Search the element that points to the current element
     TempDel = Begin;

     if( TempDel->Next == NULL )
      {
       delete Begin;
       Begin = NULL;
       End   = NULL;
      }
     else
      {
       while (TempDel->Next->Next != NULL)
        {
         TempDel = TempDel->Next;
        };
       // Modify the pointers
       delete TempDel->Next;
       TempDel->Next = NULL;
       End = TempDel;

      };
     #endif

    N--;
    return(N);
}

// Delete all the objects of the list .....................................
template <class GENERIC>
integer List<GENERIC>::DeleteAll()
{
    integer i;
    LinkObject<GENERIC> * Temp;

    if( N == 0 ) return(0);
    for( i=N; i>0; i-- )
    {
     #ifdef DOUBLY
     if ( End != Begin )
        {
         Temp      = End;
         End       = End->Previous;
         #ifndef NDEBUG
         End->Next = NULL;
         #endif
         delete Temp;
        }
     else
        {
         delete Begin;
         Begin = NULL;
         End   = NULL;
        }
     #else
     if ( Begin != End )
        {
         Temp  = Begin;
         Begin = Begin->Next;
         delete Temp;
        }
     else
        {
         delete Begin;
         Begin = NULL;
         End   = NULL;
        };
     #endif
     N--;
    }

    return(N);
}

#ifdef WHOLE_LIST
// Exchange two LinkObjects ...............................................
template <class GENERIC>
integer List<GENERIC>::Exchange(LinkObject<GENERIC> *P1, LinkObject<GENERIC> *P2)
{
    // It could be faster if it work with pointers, but is too complicated

    GENERIC TempData;

    TempData = P1->Data;
    P1->Data = P2->Data;
    P2->Data = TempData;

    return(0);
};
#endif

// Search a LinkObject at a specified position ............................
// It returns the pointer to the LinkObject
template <class GENERIC>
LinkObject<GENERIC>* List<GENERIC>::Search( integer Position )
{
    LinkObject<GENERIC> *Temp;
    integer i;

    // Search the pointer in the specified Position.
    // This is a first version using a sequential search
    //  from the nearest edge

    #ifdef DOUBLY
    if( Position <= N/2 )
    {
     Temp = Begin;
     for( i=1;   i<Position ; i++) Temp = Temp->Next;
    }
    else
    {
     Temp = End;
     for( i=N-Position; i>0 ; i--) Temp = Temp->Previous;
    }
    #else
     Temp = Begin;
     for( i=1;   i<Position ; i++)
      Temp = Temp->Next;
    #endif
    Last = Temp;
    // Return the pointer
    return Temp;
}

#ifdef WHOLE_LIST
// QuickSort the List checking with a particular check function ...........

template <class GENERIC>
void List<GENERIC>::QuickSort(int (*MatchFunction)(GENERIC P1, GENERIC P2 ))
{
    integer npivot;
    LinkObject<GENERIC> *pivot, *old;

 // This routine has been modified by me to select the first
 //  time a randomly selected pivot. It's not the better solution
 //  but it enhances the performance in a 300% for the ordered
 //  lists

    // Selects a pivot, but not at either end!

    npivot = abs( rand() ) % N;
    if( npivot<2 || npivot>N-2 ) npivot = 2;

    // Run thru the list to the randomly selected point

    old = Begin;
    while( npivot-- ) old = old->Next;
    pivot = old->Next;      // Take as pivot
    old->Next = pivot->Next;    // Cut from chain

    pivot->Next = Begin;
    Begin = pivot;

    // Call the function
    QuickSort_( (*MatchFunction), &Begin, NULL );
};

// QuickSort the List checking with a particular check function ...........
//  Function called recursively
//  BE CAREFULL ! This routine handle simple linked lists!!

// Quicksort for linked lists. Eliminates tail recursion to minimize
//  stack usage.
// This routine would not do well on either forward- or reverse-ordered
//  lists. This type of lists can produce degenerate performance.
// Bibliography : Practicals Algorithms for Programmers. Andrew Binstock &
//          John Rex. Ed: Addison Wesley. ISBN 0-201-63208-X
template <class GENERIC>
void List<GENERIC>::QuickSort_(int (*MatchFunction)(GENERIC P1, GENERIC P2 ),
                   LinkObject<GENERIC> **Head,
                   LinkObject<GENERIC>  *End  )
{
    integer left_count, right_count, count;
    LinkObject<GENERIC> **left_walk, *pivot, *old;
    LinkObject<GENERIC> **right_walk, *right;

    if( *Head != End )
    do {
        pivot     = *Head;  // Take first element as pivot
        left_walk =  Head;  // Set up left & right halves
        right_walk= &right;
        left_count= right_count = 0;

        // Now walk the list

        for( old = (*Head)->Next; old != End; old = old->Next )
         {

           if( (*MatchFunction)(old->Data, pivot->Data) < 0 )
        {
         // Less than pivot, so goes on left
         left_count  += 1;
         *left_walk   = old;
         left_walk    = &(old->Next);

        }
           else
        {
         // Greater than or equal, so goes on right
         right_count += 1;
         *right_walk  = old;
         right_walk   = &(old->Next);
        }
         };

         // Now glue the halves together ...

         *right_walk  = End;    // Terminate right list
         *left_walk   = pivot;  // Put pivot after things on left
         pivot->Next  = right;  // And rigth list after that

         // Now sort the halves in more detail

         if( left_count > right_count )
          {
           QuickSort_( (*MatchFunction), &(pivot->Next), End );
           End   = pivot;
           count = left_count;
          }
         else
          {
           QuickSort_( (*MatchFunction), Head, pivot );
           Head  = &(pivot->Next);
           count = right_count;
          };
        }
    while (count > 1);

};
#endif

// Sort List by selection method //////////////////////////////////////////
template <class GENERIC>
integer List<GENERIC>::Sort( double (*fMatch)     // Match function
                ( GENERIC &E, Vector &D ),
            Vector Distance,        // Argument function
            int Order           //  1 Ascendent
              )                 // -1 Descendent
{
 List<GENERIC>      NewList;
 LinkObject<GENERIC>    *Pointer;
 LinkObject<GENERIC>    *MinPointer;
 double           Depth, MinDepth;
 integer        i;

 // For debugging
 assert( Order==1 || Order ==-1 );
 assert( fMatch != NULL );

 if( N == 0 ) return(0);            // If no elements then return

 do
 {
   MinDepth = 1.0e30;
   Pointer  = Begin;                // The first element

   // Search the deeptest element
   for( i=1; i<=N; i++ )
    {
      Depth = fMatch( Pointer->Data, Distance );

      if ( Depth < MinDepth )
      {
       MinDepth   = Depth;
       MinPointer = Pointer;
      };

      Pointer = Pointer->Next;          // Next element
    }

   // Add this element to the new list and delete it of the old list

   switch (Order)
   {
    case -1 :
         NewList.InsertAt( 1, MinPointer->Data );
         break;
    default :
         NewList.Add( MinPointer->Data );
   };

   Delete( MinPointer );    // OPTIMIZABLE

 } while (N != 0);

 *this = NewList;
 return( N );
}

// Return the NEAREST element given by the MATCH function
//  Typically the Match function will be a friend function called ByDistance
template <class GENERIC>
GENERIC List<GENERIC>::Nearest( double (*fMatch)(GENERIC &D, Vector &V),
                Vector R )
{
 LinkObject<GENERIC>* TheTarget;
 LinkObject<GENERIC>* NearestAtom;

 double D, Minimum = 1e10;

 TheTarget = PointerToFirstData();
 NearestAtom = TheTarget;
 do
 {
  D = fMatch( TheTarget->Data, R );
  if( D < Minimum )
   {
     NearestAtom = TheTarget;
     Minimum     = D;
   };

  TheTarget = PointerToNextData();
 }
 while( TheTarget!=NULL );
 return(NearestAtom->Data);
}

#ifdef WHOLE_LIST

template <class GENERIC>
integer List<GENERIC>::Save( char* FileName, long int NN,
                 double (*fMatch)( GENERIC &A, Vector V ) )
{
 integer i = 0;
 LinkObject<GENERIC> * Temp;
 ofstream Out;

 Out.open( FileName, ios::out|ios::app );
 if( !Out ) 
 {
  cerr << "# Error: List::Save, Can't open the file " << FileName << endl;
  exit(WRITE_ERROR);
 }
 else
 {
  Temp = Begin;
  while ( i < N )
  {
   Out.width(8);
   Out << (i+NN) << " " ;
   Out.width(10);
   Out << 0.1* fMatch( Temp->Data, XVector ) << endl;

   Temp = Temp->Next;
   i++;
  }
 };

 Out.close();
 return( N );
};

// Verify the integrity of the list
template <class GENERIC>
int List<GENERIC>::Verify()
{
 LinkObject<GENERIC>* TheTarget;
 int Index = 0;

 TheTarget = PointerToFirstData();
 do
 {
  Index ++;
  TheTarget = PointerToNextData();
 }
 while( TheTarget!=NULL );

 if( Index != N ) cerr << "X";
 return( Index == N  );
};
#endif
}

#endif  // _INCLUDE_LIST_H_
