/***

    ./src/pmmconvert.c 

    Protein Motion TK (pmtk) - A library and several tools to estimate protein motions from varying input
    Copyright (C) 2009-2014  Benjamin Falkner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


*/



#include <stdio.h>
#include <string.h>
#include "config.h"
#include "pmtk.h"
#include "calc.h"

static char help [] = 
HELP_HEAD("Convert Model")
USAGE("pmmconvert","OPTION","PDB/XTC File") 
OPTIONS_HEAD
"  --help -h                 Display this help\n"
"  --out -o <model file>     Write average model to <model file>\n"
"  --structure -s <pdb file> Structure reference\n\n";

static char version [] =
VERSION_STRING("Convert PCA");


int main(int argc,  char *argv[])
{
        int i=1;
	char *fnout=0;
	PMProtein *data = pmCreateNew(PM_TYPE_MODEL);


	while(i<argc){
                if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")){
                        fputs(help,stdout);
                        exit(0);
                }
                else if(!strcmp(argv[i],"--version")){
                        fputs(version,stdout);
                        exit(0);
                }
                else if(!strcmp(argv[i],"-s") || !strcmp(argv[i],"--structure")){
                        if(argv[i+1][0]!='-'){
                                i++;
				SAY("Reference: %s ...",argv[i]);
                        	if(!pmReadRefM((PMProteinModel*)data,argv[i])){
					 WARN("Can not read \"%s\" <- will be ignored",argv[i]);
				}
			}
                }
                else if(!strcmp(argv[i],"-o") || !strcmp(argv[i],"--out")){
                        if(argv[i+1][0]!='-')
                                fnout=argv[++i];
                }
                else if(argv[i][0]!='-'){
			SAY("try reading: %s ...",argv[i]);
                        if(!pmReadM((PMProteinModel*)data,argv[i])){
				 WARN("Can not read \"%s\" <- will be ignored",argv[i]);
                        }
                }
                i++;
	}

	if(fnout){
		pmWriteM((PMProteinModel*)data,fnout);
	}
}
