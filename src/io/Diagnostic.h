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

#include <sstream>
#include <iostream>
#include "domains/Global.h"

#define LOWMSG(msg)     do { std::stringstream _s3; _s3 << msg << std::endl; Domains::Global::out(0, _s3.str()); } while(false)
#define MEDMSG(msg)     do { std::stringstream _s3; _s3 << msg << std::endl; Domains::Global::out(1, _s3.str()); } while(false)
#define HIGHMSG(msg)    do { std::stringstream _s3; _s3 << msg << std::endl; Domains::Global::out(2, _s3.str()); } while(false)
#define LOWMSG2(msg)    do { std::stringstream _s3; _s3 << msg; Domains::Global::out(0, _s3.str()); } while(false)
#define ERRORMSG(msg)   do { std::stringstream _s3; _s3 << msg; Domains::Global::error(_s3.str()); } while(false)
#define WARNINGMSG(msg) do { std::stringstream _s3; _s3 << msg; Domains::Global::warning(_s3.str()); } while(false)

#define DEBUGLOWMSG(msg) do{ if(Domains::global()->debug()) { std::stringstream _s3; _s3 << msg << std::endl; Domains::Global::out(0, _s3.str());} } while(false)

//the end
