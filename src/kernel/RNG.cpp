/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/
                     
#include "RNG.h"
#include <cassert>
#include "io/Diagnostic.h"

namespace Kernel
{

/* Period parameters */  
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* initializes mt[N] with a seed */
RNG::RNG(unsigned s) : _mti(RNG_N+1)
{
    _mag01[0]= 0x0UL; _mag01[1] = MATRIX_A;    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    _mt[0]= s & 0xffffffffUL;
    for (_mti=1; _mti<RNG_N; _mti++) {
        _mt[_mti] = 
	    (1812433253UL * (_mt[_mti-1] ^ (_mt[_mti-1] >> 30)) + _mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        _mt[_mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

	/* generates a random number on [0,0xffffffff]-interval */
unsigned RNG::genrand_int32(void)
{
    unsigned  y;
    if (_mti >= RNG_N) { /* generate N words at one time */
        int kk;

        if (_mti == RNG_N+1)   /* if init_genrand() has not been called, */
            ERRORMSG("RNG: Should not reach here"); /* a default initial seed is used */
        for (kk=0;kk<RNG_N-RNG_M;kk++) {
            y = (_mt[kk]&UPPER_MASK)|(_mt[kk+1]&LOWER_MASK);
	    assert(kk < RNG_N || kk+RNG_M <RNG_M || kk+1 <RNG_N || (y & 0x1UL) <2);
            _mt[kk] = _mt[kk+RNG_M] ^ (y >> 1) ^ _mag01[y & 0x1UL];
	    //new char;
        }
        for (;kk<RNG_N-1;kk++) {
            y = (_mt[kk]&UPPER_MASK)|(_mt[kk+1]&LOWER_MASK);
            _mt[kk] = _mt[kk+(RNG_M-RNG_N)] ^ (y >> 1) ^ _mag01[y & 0x1UL];
        }
        y = (_mt[RNG_N-1]&UPPER_MASK)|(_mt[0]&LOWER_MASK);
        _mt[RNG_N-1] = _mt[RNG_M-1] ^ (y >> 1) ^ _mag01[y & 0x1UL];
        _mti = 0;
    }
    y = _mt[_mti++];
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
}

/* generates a random number on [0,1)-real-interval */
double RNG::rand(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

}
