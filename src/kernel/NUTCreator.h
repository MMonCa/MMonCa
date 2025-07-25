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

#ifndef KMCNUTCREATOR_H
#define KMCNUTCREATOR_H

#include "Coordinates.h"
#include <vector>
#include <array>

namespace Kernel
{
class Domain;

class NUTCreator
{
public:
    NUTCreator(const Domain *, Coordinates &m, Coordinates &M,
               std::vector<float> const * const aLinesX, std::vector<float> const * const aLinesY, std::vector<float> const * const aLinesZ);
    ~NUTCreator();
    const std::vector<float> & getLines(int dim) const  { return _lines[dim];  }
    float getMinLine() const;
private:
    std::array<std::vector<float>, 3u> _lines;
};
}

#endif
