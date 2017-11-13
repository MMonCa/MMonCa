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
/*
 * message.h
 *
 * i. romero, september 2001
 *
 * messages to screen and files
 */

#ifndef _message_h
#define _message_h

#include <stdio.h>

#ifdef __cplusplus
extern "C"{
#endif

/* set options for messages */
void  ActivateDebugMessages(void);
void  DeactivateDebugMessages(void);
void  MessagesToScreen(int onoff);  /* 1-on, 0-off */
void  MessagesToFile(int onoff);

/* output message */
void  CenteredMessageToFile(FILE *fp, const char *fmt, ...); 
void  CenteredMessage(const char *fmt, ...); 
void  DebugMessage(const char *fmt, ...);
void  ErrorMessage(const char *fmt, ...);
void  Message(const char *fmt, ...);
void  MessageToScreen(const char *fmt, ...);
void  MessageToFile(const char *fmt, ...);
void  WarningMessage(const char *fmt, ...);
void  WelcomeMessage(FILE *fp);


void SetDebug_Level(int level);
int  Debug_Level(void);


#ifdef __cplusplus
}
#endif
	

#endif

