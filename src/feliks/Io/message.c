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
 * message.c
 *
 * i. romero, september 2001
 *
 * coordinates message output
 *
 */


#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "Main/feliks.h"
#include "Io/message.h"
#define  LINELENGTH 80
#define IRISMAX(a,b)   ( (a>b) ? a : b )
#define IRISMIN(a,b)   ( (a<b) ? a : b )

/* module variables  and default values*/
static int toscreen=1;
static int tofile=0;
static int debug=0;
static int debuglev = 0;

void SetDebug_Level(int level)
{
    debuglev = level;
    MessageToScreen("\nDebug level set to %d", level);
}



int  Debug_Level(void)
{
  return debuglev;
}



void ActivateDebugMessages(void)
{
	Message("Debug mode on");
	debug=1;
}



/*
 * Message that is printed centered in the line of size provided
 */
void CenteredMessageToFile(FILE *fp, const char *fmt, ...)
{
	va_list args;
	static char message[100], mlen, nspaces;
		
	// write the whole message in the variable "message", which has memory allocated
	va_start(args, fmt);
	vsprintf(message, fmt, args);
	va_end(args);
	
	// rewrite the message in a new string, now centered.
	mlen    = strlen(message);
	nspaces = (LINELENGTH - mlen)/2;	
	fprintf(fp, "%*s", IRISMAX(0,nspaces)+mlen, message);
}




void CenteredMessage(const char *fmt, ...)
{
	va_list args;
	size_t   mlen, nspaces;
	static char message[100]="";

	// print the original string in variable "message"
	va_start(args, fmt);
	vsprintf(message, fmt, args);
	va_end(args);
	
	// measure the length of the message and compute the additional space
	mlen    = strlen(message);
	nspaces = (LINELENGTH - mlen)/2;
	
	//output the original string, as well as the additional spaces
	Message("\n%*s", IRISMAX(0,nspaces)+mlen, message);
}



void DeactivateDebugMessages(void)
{
  debug=0;
}


void DebugMessage(const char *fmt, ...)
{
  va_list args , argscopy;
  char *newfmt=NULL;
  int   i = 0;

  
  if (debug > 0)
    {
      newfmt = (char *) malloc((strlen(fmt) + 20)* sizeof(char));

      /* eliminate '\n' characters at the begining of the format*/
      while(fmt[i] == '\n')
	  newfmt[i++] = '\n';

      /* put fmt in newfmt, adding at begining and end some chars */
      strcpy(newfmt+i,"\n\t*** ");
      strcat(newfmt, fmt+i);
      va_start(args , fmt);
      va_start(argscopy, fmt);

		vprintf(newfmt, args);
		fflush(stdout);

      va_end(args);
      va_end(argscopy);
      free(newfmt);
    }
}



void ErrorMessage(const char *fmt, ...)
{
  va_list args , argscopy;
  char *newfmt=NULL;
  int   i = 0;

  newfmt = (char *) malloc( (strlen(fmt) + 20)* sizeof(char));
  
  /* eliminate '\n' characters at the begining of the format*/
  while(fmt[i] == '\n')
    newfmt[i++] = '\n';
  
  /* put fmt in newfmt, adding at begining and end some chars */
  strcpy(newfmt+i,"\n****** Error :");
  strcat(newfmt, fmt+i);
  strcat(newfmt, " ******\n\n\n");
  va_start(args , fmt);
  va_start(argscopy, fmt);

  vprintf(newfmt, args);
  //vToLogFile(newfmt, argscopy);
  
  va_end(args);
  va_end(argscopy);
  free(newfmt);
  
  exit(1);
}







/*
 * this function just redirects the output of a printf command to the screen
 * and the log file  that must have been already opened. Uses the standard function vfprintf
 * that works like fprintf but for a va_list argument. See K&R pg 174
 */
void Message(const char *fmt, ...)
{
  va_list args;
  
  if (toscreen == 1) 
    {
      va_start(args , fmt);
      vprintf(fmt, args);
      va_end(args);
    }

  if (tofile   == 1) 
  {
    va_start(args , fmt);
    //vToLogFile(fmt , args);
    va_end(args);
  }
}



/* sends message to file only, if file output is activated 
 */
void MessageToFile(const char *fmt, ...)
{
  va_list args;
  
  if (tofile == 1) 
  {
    va_start(args, fmt);
    //vToLogFile(fmt, args);
    va_end(args);
  }
}



/* sends message to screen only 
 */
void MessageToScreen(const char *fmt, ...)
{
  va_list args;
  
  if (toscreen == 1) 
  {
    va_start(args , fmt);
    //vToLogFile(fmt , args);
    va_end(args);
  }
}



void MessagesToFile(int onoff)
{
  if (onoff == 1) tofile = 1;
  else            tofile = 0;
}



void MessagesToScreen(int onoff)
{
  if (onoff == 1) toscreen = 1;
  else            toscreen = 0;
}




void WarningMessage(const char *fmt, ...)
{
    va_list args , argscopy;
    char *newfmt=NULL;
    int   i = 0;

    newfmt = (char *) malloc( (strlen(fmt) + 20)* sizeof(char));

    /* eliminate '\n' characters at the begining of the format*/
    while(fmt[i] == '\n') i++;

    /* put fmt in newfmt, adding at begining and end some chars */
    strcpy(newfmt+i,"\n  ** WARNING: ");
    strcat(newfmt, fmt+i);
    strcat(newfmt, "\n");
    va_start(args , fmt);
    va_start(argscopy, fmt);
    //vToLogFile(newfmt , argscopy);

    va_end(args);
    va_end(argscopy);
    free(newfmt);
}



void WelcomeMessage(FILE *fp)
{
  fprintf(fp,"\n\n\n\n\n\n");
  fprintf(fp,"   -------------------------------------------------------------------");
  fprintf(fp,"\n\n                             F    E    LI    K   S");
  fprintf(fp,"\n\n\n                               version %4.2f", FELIKS_VERSION);
  fprintf(fp,"\n                   Ignacio Romero - iromero@mecanica.upm.es");
  fprintf(fp,"\n                      Universidad Politecnica de Madrid");
  fprintf(fp,"\n\n");
  fprintf(fp,"   -------------------------------------------------------------------");
  fprintf(fp,"\n\n\n");
}
