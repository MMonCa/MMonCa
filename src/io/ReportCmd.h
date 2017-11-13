/*
 * ReportCmd.h
 *
 *  Created on: Jun 14, 2012
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

#ifndef REPORTCMD_H_
#define REPORTCMD_H_

#include "Command.h"

namespace IO {

class ReportCmd: public IO::Command {
public:
	ReportCmd(Tcl_Interp *, int argc, const char *argv[]);
	~ReportCmd() { }

	virtual int operator()();
};

} /* namespace IO */
#endif /* REPORTCMD_H_ */
