/*
 * Author: ignacio.martin@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain
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

//Wrapper for having the RNG as a class

#ifndef _2_KERNEL_RNG__
#define _2_KERNEL_RNG__

#define RNG_N 624
#define RNG_M 397

namespace Kernel
{
class RNG
{
public:
    RNG(unsigned seed = 4357);
    double rand();
    unsigned genrand_int32(void);
    
private:
    unsigned _mag01[2];
    unsigned _mt[RNG_N]; //the array for the state vector 
    int _mti;
};
}

#endif
