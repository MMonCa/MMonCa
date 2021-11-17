/*
 * UnitTestCmd.h
 *
 *  Created on: Nov 2021
 *
 * Author: ignacio.martin
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

#ifndef UNITTESTCMD_H_
#define UNITTESTCMD_H_

#include "Command.h"

namespace IO {

class UnitTestCmd: public IO::Command {
public:
	UnitTestCmd(Tcl_Interp *, int argc, const char *argv[]);
    ~UnitTestCmd() {}
    virtual int operator()();
private:
    bool test_FermiDirac();
};

}

#endif /* UNITTESTCMD_H_ */
