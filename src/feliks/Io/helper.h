/* Finite Element Method Module
 *
 * Author: ignacio.romero@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain,
 *      and       Technical University of Madrid (UPM), Madrid, Spain
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
/* helper.h
 *
 * functions to provide help in the use of feliks
 *
 * i. romero, january 2002
 */

#ifndef _helper_h
#define _helper_h

#include <string>

class helper{

private:
	void  describeFELIKS() const;
	void  helpBasics() const;
	void  commandList() const;
	
public:
	void describe(const std::string& command) const;
	void printHelpFile(const std::string& filename) const;
};


#endif
