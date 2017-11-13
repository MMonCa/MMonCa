/*
 * TestCmd.h
 *
 *  Created on: Jul 6, 2011
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

#ifndef TESTCMD_H_
#define TESTCMD_H_

#include "Command.h"

namespace IO {

class TestCmd: public IO::Command {
public:
	TestCmd(Tcl_Interp *, int argc, const char *argv[]);
    ~TestCmd() {}
    virtual int operator()();

    static std::string reformatText(const std::string &); //removes new lines and collapses too long texts for outputs.
};

}

#endif /* TESTCMD_H_ */
