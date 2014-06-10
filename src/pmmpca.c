/***

    ./src/pmmpca.c 

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
HELP_HEAD("Model PCA")
USAGE("pmmpca","OPTION","PDB/XTC File") 
OPTIONS_HEAD
"  --help -h                 Display this help\n"
"  --mean -m <model file>    Write average model to <model file>\n"
"  --variance  -v <model file>\n"
"                            Write variance model to <model file>\n"
"  --evecs -e <model file>   Write eigenvectors to <model file>\n"
"  --evals -V <file>         Write eigenvalues to <file>\n"
"  --length -l <int>         Set the number of eigenvectors to be written\n"
"  --structure -s <pdb file> Structure reference\n"
"  --evec-sacle -S <float>   Scaling eigenvectors\n"
"  --align -a                Align data\n\n";

static char version [] =
VERSION_STRING("Model PCA");


int main(int argc,  char *argv[])
{
        int i=1;
	int align=0;
        size_t length=0;
        PMProtein *data,*mean,*evecs;
        float *evals=0;
        char *fnevecs=0, *fnmean=0, *fnvar=0, *fnevals=0;
	float S=1.0;
	
	data = pmCreateNew(PM_TYPE_MODEL);
        if(!data)
		FAIL("MALLOC: no memmory assigned")
	
	//loop options	
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
                else if(!strcmp(argv[i],"-m") || !strcmp(argv[i],"--mean")){
                        if(argv[i+1][0]!='-')
                                fnmean=argv[++i];
                }
                else if(!strcmp(argv[i],"-v") || !strcmp(argv[i],"--variance")){
                        if(argv[i+1][0]!='-')
                                fnvar=argv[++i];
                }
                else if(!strcmp(argv[i],"-S") || !strcmp(argv[i],"--evec-scale")){
                        if(argv[i+1][0]!='-')
                                S=atof(argv[++i]);
                }
                else if(!strcmp(argv[i],"-e") || !strcmp(argv[i],"--evecs")){
                        if(argv[i+1][0]!='-')
                                fnevecs=argv[++i];
                }
                else if(!strcmp(argv[i],"-V") || !strcmp(argv[i],"--evals")){
                        if(argv[i+1][0]!='-')
                                fnevals=argv[++i];
                }
                else if(!strcmp(argv[i],"-l") || !strcmp(argv[i],"--length")){
                        if(argv[i+1][0]!='-')
                                length=atoi(argv[++i]);
                }
                else if(!strcmp(argv[i],"-a") || !strcmp(argv[i],"--align")){
			align=1;
                }
		else if(argv[i][0]=='-'){
			WARN(" [Unkown argumant \"%s\"] <- will be ignored",argv[i]);
		}
                else if(argv[i][0]!='-'){
			SAY("try reading: %s ...",argv[i]);
                        if(!pmReadM((PMProteinModel*)data,argv[i])){
				 WARN("Can not read \"%s\" <- will be ignored",argv[i]);
                        }
                }
                i++;
	}


	pmPrintInfo(data);


	if(align){
		SAY("Alignment ...")
		pmCenter((PMProteinModel*)data,(PMProteinModel*)data);
		pmAlignTo((PMProteinModel*)data,(PMProteinModel*)data,(PMProteinModel*)data);
	}


	if(fnvar){
		SAY("Calculating mean and variance ...");
		mean = pmCreatefReff(data,2);
		if(!mean)
			FAIL("MALLOC: no memmory assigned")
		if(!pmMeanVar(mean,data))
			FAIL("CALC: calculation returned null pointer")
		if(fnmean) SAY("Writing %s ...",fnmean);
			pmWriteSM((PMProteinModel*)mean,fnmean,0);
		SAY("Writing %s ...",fnvar);
		pmWriteSM((PMProteinModel*)mean,fnvar,1);
	}
	else {
		SAY("Calculating mean ...");
		mean = pmCreatefReff(data,1);
		if(!mean)
			FAIL("MALLOC: no memmory assigned")
		if(!pmMean(mean,data))
			FAIL("CALC: calculation returned null pointer")
		if(fnmean) SAY("Writing %s ...",fnmean);
			pmWriteSM((PMProteinModel*)mean,fnmean,0);
	}

	if(fnevecs || fnevals){		
		SAY("Calculating pca ...");
		if(length<=0) length=pmGetDF(data);
		evecs = pmCreatefReff(data,length);
		evals = pmPCA(evecs,data,mean);
		if(S!=1.0){
			pmScale(evecs,evecs,S);
		}
		if(fnevecs) 
			SAY("Writing %s ...",fnevecs);
			pmWriteM((PMProteinModel*)evecs,fnevecs);
			
		if(fnevals){
			SAY("Writing %s ...",fnevals);
			writeArray(fnevals,"Eigenvalues",evals,length);
		}
	}
	return 0;
}

