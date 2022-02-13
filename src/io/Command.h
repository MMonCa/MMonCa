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

#ifndef GENERALCOMMAND_H
#define GENERALCOMMAND_H

#include <vector>
#include <map>
#include <tcl.h>
#include "Parameters.h"

namespace Kernel
{
	class Domain;
}

namespace IO {

class Command : public Parameters
{
public:
    typedef std::map<std::string, std::string> str2str; 
    static int Anneal     (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Cascade    (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Cluster    (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Init       (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Insert     (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Save       (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Exit       (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Extract    (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int LowMsg     (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int License    (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Param      (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Profile    (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Report     (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Restart    (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Test       (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int Charge	  (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    static int UnitTest	  (ClientData, Tcl_Interp *interp, int argc, const char *argv[]);
    
    Command(Tcl_Interp *, int argc, const char *argv[], bool bPrint=true);
    virtual ~Command();

    virtual int operator()() = 0;

protected:    
    void separateTokens(const char *str, char delimiter, std::vector<std::string> &vec);
    
    Tcl_Interp *_pTcl;
};

}

#endif
