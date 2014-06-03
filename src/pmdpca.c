/***

    ./src/pmdpca.c 

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
HELP_HEAD("Density PCA")
USAGE("pmdpca","OPTION","MRC File") 
OPTIONS_HEAD
"  --help -h                 Display this help\n"
"  --mean -m <density file>  Write average density to <density file>\n"
"  --variance  -v <density file>\n"
"                            Write variance density to <density file>\n"
"  --evecs -e <density file> Write eigenvectors to <density file>\n"
"  --evals -V <file>         Write eigenvalues to <file>\n"
"  --length -l <int>         Set the number of eigenvectors to be written\n\n";

static char version [] =
VERSION_STRING("Density PCA");


int main(int argc,  char *argv[])
{
        int i=1;
        size_t length=0;
        PMProtein *data,*mean,*evecs;
        float *evals=0;
        char *fnevecs=0, *fnmean=0, *fnvar=0, *fnevals=0;
	
	data = pmCreateNew(PM_TYPE_DENSITY);
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
                else if(!strcmp(argv[i],"-m") || !strcmp(argv[i],"--mean")){
                        if(argv[i+1][0]!='-')
                                fnmean=argv[++i];
			SAY("mean: %s",fnmean);
                }
                else if(!strcmp(argv[i],"-v") || !strcmp(argv[i],"--variance")){
                        if(argv[i+1][0]!='-')
                                fnvar=argv[++i];
			SAY("variance: %s",fnvar);
                }
                else if(!strcmp(argv[i],"-e") || !strcmp(argv[i],"--evecs")){
                        if(argv[i+1][0]!='-')
                                fnevecs=argv[++i];
			SAY("eigenvectors: %s",fnevecs);
                }
                else if(!strcmp(argv[i],"-V") || !strcmp(argv[i],"--evals")){
                        if(argv[i+1][0]!='-')
                                fnevals=argv[++i];
                        SAY("eigenvalues: %s",fnevals);
                }
                else if(!strcmp(argv[i],"-l") || !strcmp(argv[i],"--length")){
                        if(argv[i+1][0]!='-')
                                length=atoi(argv[++i]);
			SAY("number of eigenvectors: %lu",length);		
                }
		else if(argv[i][0]=='-'){
			WARN("Unkown argumant will be ignored \"%s\"",argv[i]);
		}
                else if(argv[i][0]!='-'){
			SAY("try reading: %s ...",argv[i]);
                        if(!pmReadD((PMProteinDensity*)data,argv[i])){
				 WARN("Can not read \"%s\" - will be ignored",argv[i]);
                        }
                }
                i++;
	}
	if(fnvar){
		SAY("Calculating mean and variance ...");
		mean = pmCreatefReff(data,2);
		if(!mean)
			FAIL("MALLOC: no memmory assigned")
		if(!pmMeanVar(mean,data))
			FAIL("CALC: calculation returned null pointer")
		if(fnmean) SAY("Writing %s ...",fnmean);
			pmWriteSD((PMProteinDensity*)mean,fnmean,0);
		SAY("Writing %s ...",fnvar);
		pmWriteSD((PMProteinDensity*)mean,fnvar,1);
		SAY("done");
	}
	else {
		SAY("Calculating mean ...");
		mean = pmCreatefReff(data,1);
		if(!mean)
			FAIL("MALLOC: no memmory assigned")
		if(!pmMean(mean,data))
			FAIL("CALC: calculation returned null pointer")
		if(fnmean) SAY("Writing %s ...",fnmean);
			pmWriteSD((PMProteinDensity*)mean,fnmean,0);
		SAY("done");
	}

	if(fnevecs || fnevals){		
		SAY("Calculating pca ...");
		if(length<=0) length=pmGetDF(data);
		evecs = pmCreatefReff(data,length);
		evals = pmPCA(evecs,data,mean);
		if(fnevecs) 
			SAY("Writing %s ...",fnevecs);
			pmWriteD((PMProteinDensity*)evecs,fnevecs);
			
		if(fnevals){
			SAY("Writing %s ...",fnevecs);
			writeArray(fnevecs,"Eigenvalues",evals,evecs->pany.records);
		}
		SAY("done");
	}
	return 0;
}

