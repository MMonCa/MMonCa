/*
 * ClusterCmd.h
 *
 *  Created on: May 22, 2012
 *      Author: ignacio.martin@imdea.org
 *
 *      This command is more a "TCL" command/helper than a real KMC command
 *      That is why it is not inherited from Command. Its use does not rely on the Command=type syntax.
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

#ifndef CLUSTERCMD_H_
#define CLUSTERCMD_H_

#include <tcl.h>

namespace IO {

class ClusterCmd
{
public:
	int operator()(Tcl_Interp *, int argc, const char *argv[]);
};

} /* namespace IO */
#endif /* CLUSTERCMD_H_ */
