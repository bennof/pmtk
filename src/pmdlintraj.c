/***

    ./src/pmdlintraj.c 

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
HELP_HEAD("Linear Trajectory Density")
USAGE("pmmlintraj","OPTION","PDB/XTC File") 
OPTIONS_HEAD
"  --help -h                 Display this help\n"
"  --out -o <model file>     Write average model to <model file>\n"
"  --center -c <model file>  Center of the linear trajectory (i.e. mean)\n"
"  --evec -e <model file>    Eigenvector <model file>\n"
"  --idxevec -i <int>        Index of the eigenvector\n"
"  --frames -f <int>         Number of frames that will be calculated\n"
"  --peak -p <float>         Peak value of the trajectory\n\n";

static char version [] =
VERSION_STRING("Linear Trajectory Density");


int main(int argc,  char *argv[])
{
        int i=1;
	char *fnout=0;
	PMProtein *evecs = pmCreateNew(PM_TYPE_DENSITY);
	PMProtein *center = pmCreateNew(PM_TYPE_DENSITY);
	PMProtein *out = pmCreateNew(PM_TYPE_DENSITY);
	size_t cframes=10;
	float peak=10.;


	while(i<argc){
                if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")){
                        fputs(help,stdout);
                        exit(0);
                }
                else if(!strcmp(argv[i],"--version")){
                        fputs(version,stdout);
                        exit(0);
                }
                else if(!strcmp(argv[i],"-o") || !strcmp(argv[i],"--out")){
                        if(argv[i+1][0]!='-')
                                fnout=argv[++i];
                }
                else if(!strcmp(argv[i],"-c") || !strcmp(argv[i],"--center")){
                        if(argv[i+1][0]!='-'){
                                i++;
				SAY("Center: %s ...",argv[i]);
                        	if(!pmReadD((PMProteinDensity*)center,argv[i])){
					 WARN("Can not read \"%s\" <- will be ignored",argv[i]);
				}
			}
                }
                else if(!strcmp(argv[i],"-e") || !strcmp(argv[i],"--evec")){
                        if(argv[i+1][0]!='-'){
                                i++;
				SAY("Eigenvector: %s ...",argv[i]);
                        	if(!pmReadD((PMProteinDensity*)evecs,argv[i])){
					 WARN("Can not read \"%s\" <- will be ignored",argv[i]);
				}
			}
                }
                else if(!strcmp(argv[i],"-f") || !strcmp(argv[i],"--frames")){
                        if(argv[i+1][0]!='-')
                                cframes=atol(argv[++i]);
                }
                else if(!strcmp(argv[i],"-p") || !strcmp(argv[i],"--peak")){
                        if(argv[i+1][0]!='-')
                                peak=atof(argv[++i]);
                }
		else{
			WARN(" [Unkown argumant \"%s\"] <- will be ignored",argv[i]);
                }
                i++;
	}

	
	pmCreatefReff(center,cframes);
        pmLinTraj(out,center,evecs,0,peak);

	if(fnout){
		pmWriteD((PMProteinDensity*)out,fnout);
	}
}
