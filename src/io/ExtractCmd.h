/*
 * ExtractCmd.h
 *
 *  Created on: Feb 15, 2011
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
 */

#ifndef EXTRACTCMD_H_
#define EXTRACTCMD_H_

#include "Command.h"

namespace IO
{
class ExtractCmd: public Command
{
public:
	ExtractCmd(Tcl_Interp *, int argc, const char *argv[]);
	virtual ~ExtractCmd() {}
    virtual int operator()();
};
}

#endif /* EXTRACTCMD_H_ */
