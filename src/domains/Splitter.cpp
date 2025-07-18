/*
 * Splitter.cpp
 *
 *  Created on: Oct 15, 2013
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

#include "Splitter.h"
#include "Global.h"
#include "io/FileParameters.h"

namespace Domains
{
	Splitter::Splitter(const Kernel::Coordinates &m, const Kernel::Coordinates &M)  //global simulation size
			 : _minCell(m), _maxCell(M)
	{
		_levels = global()->getFileParameters()->getInt("MC/General/levels");
	}
}
