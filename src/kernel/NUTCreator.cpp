/* Author: ignacio.martin@imdea.org
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

#include "NUTCreator.h"
#include "io/Diagnostic.h"
#include "Domain.h"
#include "lkmc/Lattice.h"
#include "MeshParam.h"
#include <cmath>

using std::string;
using std::vector;
using std::pair;

namespace Kernel {

NUTCreator::NUTCreator(const Domain *p, Coordinates &m, Coordinates &M,
               std::vector<float> const * const aLinesX, std::vector<float> const * const aLinesY, std::vector<float> const * const aLinesZ)
{
        if(aLinesX != nullptr && aLinesX->size() > 1u && aLinesY != nullptr && aLinesY->size() > 1u && aLinesZ != nullptr && aLinesZ->size() > 1u) {
            _lines[0] = *aLinesX;
            _lines[1] = *aLinesY;
            _lines[2] = *aLinesZ;
            LOWMSG("     X: (" << m._x << " - " << M._x << ") nm. " << aLinesX->size() - 1 << " elements. Custom line spacing.");
            LOWMSG("     Y: (" << m._y << " - " << M._y << ") nm. " << aLinesY->size() - 1 << " elements. Custom line spacing.");
            LOWMSG("     Z: (" << m._z << " - " << M._z << ") nm. " << aLinesZ->size() - 1 << " elements. Custom line spacing.");
        }
        else {
            for(int i=0; i<3; ++i)
            {
                _delta[i] = p->_pMePar->_spacing[i];
                int nBoxes = int(std::floor((M[i] - m[i])/_delta[i] + .5));
                if(nBoxes < 2)
                    ERRORMSG("Less than 2 boxes in dimension " << i);
                double inc = (M[i] - m[i])/double(nBoxes);
                for(int j=0; j<nBoxes; ++j) {
                    _lines[i].push_back(m[i] + j*inc);
                }
                _lines[i].push_back(M[i]);
                std::string axis = "X";
                axis[0] += i;
                LOWMSG("     " << axis << ": (" << m[i] << " - " << M[i] <<
                 ") nm. " << nBoxes << " elements. Delta = " << inc << " nm.");
            }
        }
	LOWMSG("Total " << (_lines[0].size()-1)*(_lines[1].size()-1)*(_lines[2].size()-1) << " elements" );
}


NUTCreator::~NUTCreator()
{
}

}
