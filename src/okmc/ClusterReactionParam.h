/*
 * ClusterReactionParam.h
 *
 *  Created on: Feb 25, 2014
  *
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
 *
 *  Defines a structure for defect interaction. _interactions[def0][def1] points to the ClusterInter
 *  ClusterInter is passed two sizes, for def0 and def1
 *  _defectType[0] contains the result when size0 == size1
 *  _defectType[1] contains the result when size0 ~= size1
 *  _defectType[2] contains the result when size0 <  size1
 *  _defectType[3] contains the result when size0 >  size1
 *
 */

#ifndef CLUSTERREACTIONPARAM_H_
#define CLUSTERREACTIONPARAM_H_

#include <vector>
#include <string>
#include "kernel/Material.h"
namespace IO   { class FileParameters; }
namespace OKMC { class ClusterParam; }

namespace OKMC {

struct ClusterInter
{
	ClusterInter() { for(unsigned i=0; i<4;++i) _defectType[i] = -1; }
	int _defectType[4];  // -1 mins not defined. 0 ==, 1 ~=, 2 <, 3 >. The value is the result if the boolean operation is done
	float _approx;
	int operator()(unsigned size0, unsigned size1) const;
};

class ClusterReactionParam
{
public:
	ClusterReactionParam(const ClusterParam * pClPar, const IO::FileParameters * pPar);
	std::vector<std::vector<ClusterInter> > _interactions[Kernel::MAX_MATERIALS];

	static void warning(std::string param, Kernel::M_TYPE mt, unsigned def0, unsigned def1, std::string op1, std::string op2);
};

} /* namespace OKMC */
#endif /* CLUSTERREACTIONPARAM_H_ */
