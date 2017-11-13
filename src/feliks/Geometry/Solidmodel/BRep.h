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
//
//  BRep.h
//  feliks
//
//  Created by Romero Ignacio on 11/9/11.
//  Copyright (c) 2011 Universidad Politécnica de Madrid. All rights reserved.
//

#pragma once
#ifndef feliks___brep_h
#define feliks___brep_h

#include "Math/Topology/topology.h"
#include <map>
#include <vector>

namespace feliks
{
    namespace math
    {
        class R3embedding;
        class R3submersion;
    }
    namespace geometry
    {
        class BRep
        {
        public:
            BRep(){}
            virtual ~BRep(){}
            
            virtual const   topology::cellcomplex&                  getCellcomplex() const=0;
            virtual         topology::cellcomplex&                  getCellcomplex()=0;
            virtual int                                             getDimension() const=0;
            
            const std::map<topology::cell *, math::R3embedding *>&  getEmbeddings() const{ return f;}
            std::map<topology::cell *, math::R3embedding* >&        getEmbeddings() {return f;}
            const std::map<topology::cell*, math::R3submersion*>&   getSubmersions() const {return g;}
            std::map<topology::cell*, math::R3submersion*>&         getSubmersions() {return g;}
            
            math::R3embedding&                                      getEmbedding(topology::cell& c) {return *(f[&c]);}
            math::R3submersion&                                     getSubmersion(topology::cell& c) {return *(g[&c]);}
            

            
        protected:
            std::map<topology::cell *, math::R3embedding *>  f;
            std::map<topology::cell *, math::R3submersion *> g;
        };
    }
}


#endif
