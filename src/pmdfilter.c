/***

    ./src/pmdfilter.c 

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




#include "config.h"
#include "pmtk.h"
#include <string.h>


static char help [] = 
HELP_HEAD("Density Filter")\
USAGE("pmdfilter","FILTER/Options","MRC File") \
"Filter:\n"
#ifdef USE_FFTW3
"  cutoff <lb> <ub>          Apply a cutoff at lower boudary <lb> and upper boundary <ub> in Fourier space (the values are in Angstrom)\n"
"  gauss <sigma>             Apply a Gaussian filter with a sigma in Angstrom\n"
"  laplace                   Apply a Laplace filter\n"
"  framp                     Apply a ramp filter in Fourrier space\n"
#endif
"  rramp                     Apply a ramp filter in real space\n"
"  norm                      Normalize the density\n"\
"  noise <sigma>             Add a flat band noise\n"\
"  scale <factor>            Scale the density by <factor>\n"\
"  sqrt                      Calculate the square root\n"\
"  isqrt                     Calculate the inverse square root\n"\
"  threshold <value>         Apply a threshold to cut all density values below <value>\n"\
"  mask <density file> <value>\n"\
"                            Use <density file> as a mask with a threshold at <value>\n"\
"  add <density file>        Add another density to the current one\n"\
"  sub <density file>        Substract another density from the current one\n"\
"  mult <density file>       Multiply another density with the current one\n"\
"  div <density file>        Divide current density by another one\n"\
"\n"\
"Options:\n"\
"  --help -h                 Display this help\n"\
"  --out  -o <density file>  Write output to <density file>\n\n";



static char version [] =
VERSION_STRING("Density Filter");


typedef struct _filter filter; 

struct _filter{
        filter *next;
        int id; 
        char **arg;
};

filter *start, *cur;

static void newFilter(int id, char **arg)
{
        filter *f = (filter*)malloc(sizeof(filter));
        f->id=id;
        f->next=0;
        if(start==0){
                start=f;
                cur=f;
        }
        else{
                cur->next=f;
        }
}

static inline void frewind()
{
	cur=start;
}

static inline filter* getNext()
{
	filter *c=cur;
        cur=cur->next;
        return c;
}

#define CUTOFF  1
#define GAUSS   2
#define LAPLACE 3
#define FRAMP   4
#define RRAMP   5
#define NORM    6
#define NOISE   7
#define SCALE   8
#define SQRT    9
#define ISQRT  10
#define THRES  11
#define MASK   12
#define ADD    13
#define SUB    14
#define MULT   15
#define DIV    16

	
static PMProtein *opF(const char* name)
{
	PMProtein* p;
	p = pmCreateNew(PM_TYPE_DENSITY);
        if(!p)
		FAIL("MALLOC: no memmory assigned")
	SAY("try reading: %s ...",name);
        if(!pmReadD((PMProteinDensity*)p,name)){
		WARN("Can not read \"%s\" - will be ignored",name);
       	}
	return p;
}


int main(int argc, char *argv[])
{
        int i=1;
        const char *fnout;
        PMProtein *data,*tmp;

        //parse args
	data = pmCreateNew(PM_TYPE_DENSITY);
        if(!data)
		FAIL("MALLOC: no memmory assigned")
        while(i<argc){
                //OPTIONS
                if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")){
                        fputs(help,stdout);
                        exit(0);
                }
                else if(!strcmp(argv[i],"-v") || !strcmp(argv[i],"--version")){
                        fputs(version,stdout);
                        exit(0);
                }
                else if(!strcmp(argv[i],"-o") || !strcmp(argv[i],"--out")){
                        if(argv[i+1][0]!='-')
                                fnout=argv[i++];
                }
                //FILTER
                else if (!strcmp(argv[i],"cutoff")){
                        newFilter(CUTOFF,&argv[i]);
			i+=2;
                }
                else if (!strcmp(argv[i],"gauss")){
                        newFilter(GAUSS,&argv[i]);
			i+=1;
                }
                else if (!strcmp(argv[i],"laplace")){
                        newFilter(LAPLACE,&argv[i]);
                }
                else if (!strcmp(argv[i],"framp")){
                        newFilter(FRAMP,&argv[i]);
                }
                else if (!strcmp(argv[i],"rramp")){
                        newFilter(RRAMP,&argv[i]);
                }
                else if (!strcmp(argv[i],"norm")){
                        newFilter(NORM,&argv[i]);
                }
                else if (!strcmp(argv[i],"noise")){
                        newFilter(NOISE,&argv[i]);
			i+=1;
                }
                else if (!strcmp(argv[i],"scale")){
                        newFilter(SCALE,&argv[i]);
			i+=1;
                }
                else if (!strcmp(argv[i],"sqrt")){
                        newFilter(SQRT,&argv[i]);
                }
                else if (!strcmp(argv[i],"threshold")){
                        newFilter(THRES,&argv[i]);
                }
                else if (!strcmp(argv[i],"mask")){
                        newFilter(MASK,&argv[i]);
			i+=2;
                }
                else if (!strcmp(argv[i],"add")){
                        newFilter(ADD,&argv[i]);
			i+=1;
                }
                else if (!strcmp(argv[i],"sub")){
                        newFilter(SUB,&argv[i]);
			i+=1;
                }
                else if (!strcmp(argv[i],"mult")){
                        newFilter(MULT,&argv[i]);
			i+=1;
                }
                else if (!strcmp(argv[i],"div")){
                        newFilter(DIV,&argv[i]);
			i+=1;
                }
		else {
			SAY("try reading: %s ...",argv[i]);
        		if(!pmReadD((PMProteinDensity*)data,argv[i])){
				WARN("Can not read \"%s\" - will be ignored",argv[i]);
       			 }
			
		}
                i++;
        }
	frewind();
	
	 //loop through filters
        while(getNext()){
                switch(cur->id){
#ifdef USE_FFTW3
                case CUTOFF:
			if(data->type^PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_FORWARD);
			SAY("cutoff lb=%s ub=%s [phase space]",cur->arg[1],cur->arg[2])
                        pmFFTCutOffFilter((PMProteinDensity*)data,atof(cur->arg[1]),atof(cur->arg[2]));
                        break;
                case GAUSS:
			if(data->type^PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_FORWARD);
			SAY("gauss sigma=%s [phase spacee]",cur->arg[1])
                        pmFFTGaussFilter((PMProteinDensity*)data,atof(cur->arg[1]));
                        break;
                case LAPLACE:
			if(data->type^PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_FORWARD);
			SAY("laplace edge detection [phase space]")
                        pmFFTLaplaceFilter((PMProteinDensity*)data);
                        break;
                case FRAMP:
			if(data->type^PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_FORWARD);
			SAY("RampFilter 1/r [phase space]")
                        pmFFTRampFilter((PMProteinDensity*)data);
                        break;
#endif
		case RRAMP:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("RampFilter 1/r")
                        pmRampFilter((PMProteinDensity*)data,(PMProteinDensity*)data);
                        break;
                case NORM:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("normalize mu=0.0, sigma=1.0")
                        pmNormalize((PMProteinDensity*)data,(PMProteinDensity*)data);
                        break;
                case NOISE:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("Noise %s",cur->arg[1])
                        pmAddNoise(data,data,atof(cur->arg[1]));
                        break;
                case SCALE:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("scale by %s",cur->arg[1])
                        pmScale(data,data,atof(cur->arg[1]));
                        break;
                case SQRT:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("sqrt")
                        pmSqrt(data,data);
                        break;
                case ISQRT:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("1/sqrt")
                        pmIsqrt(data,data);
                        break;
                case THRES:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("threshold at %s",cur->arg[1])
                        pmThreshold((PMProteinDensity*)data,(PMProteinDensity*)data,atof(cur->arg[1]));
                        break;
                case MASK:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("mask with %s at %s",cur->arg[1],cur->arg[2])
			tmp = opF((const char*)cur->arg[1]);
			pmMask((PMProteinDensity*)data,(PMProteinDensity*)data,(PMProteinDensity*)tmp,1,atof(cur->arg[2]));	
			pmDelete(tmp);
			break;
		case ADD:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("add %s",cur->arg[1])
			tmp = opF((const char*)cur->arg[1]);
                        pmAdd(data,data,tmp);
			pmDelete(tmp);
                        break;
                case SUB:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("substract %s",cur->arg[1])
			tmp = opF((const char*)cur->arg[1]);
                        pmSub(data,data,tmp);
			pmDelete(tmp);
                        break;
                case MULT:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("multiply by %s",cur->arg[1])
			tmp = opF((const char*)cur->arg[1]);
                        pmMult(data,data,tmp);
			pmDelete(tmp);
                        break;
                case DIV:
#ifdef USE_FFTW3
			if(data->type&PM_DENSITY_COMPLEX)
				pmFFT((PMProteinDensity*)data,PM_FFT_BACKWARD);
#endif
			SAY("divide by %s",cur->arg[1])
			tmp = opF((const char*)cur->arg[1]);
                        pmDiv(data,data,tmp);
			pmDelete(tmp);
                        break;
		}	
	}
	//output
	if(!fnout)
		fnout="out.mrc";
	SAY("Writing %s ...",fnout);
        pmWriteD((PMProteinDensity*)data,fnout);	
	SAY("done");
	
	return 0; 
}
