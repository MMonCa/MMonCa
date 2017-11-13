/*
 * FromFile.h
 *
 *  Created on: Jul 21, 2011
 *      Author: ignacio.martin@imdea.org
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

#ifndef FROMFILE_H_
#define FROMFILE_H_

#include <string>
#include <vector>
#include "kernel/Coordinates.h"

namespace OKMC   { class MobileParticle; }

namespace Cascade {

class FromFile {
public:
	FromFile();
	~FromFile();

	void operator()(const std::string &, std::vector<std::string> &format,
			float fluence,	float flux,	std::string defects,
			bool bDisplace, bool periodic, bool bReact, float temp, bool voluminic, bool bCorrectX);
};

}

#endif /* FROMFILE_H_ */
