/***

    ./src/lib/pmtk.h 

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





#ifndef PMTK_H_
#define PMTK_H_
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "vec.h"
#include "config.h"

#define PM_TYPE_MODEL         0x1
#define PM_TYPE_DENSITY       0x2
#define PM_DENSITY_COMPLEX    0x4

#define PM_FFT_FORWARD 1
#define PM_FFT_BACKWARD -1


typedef unsigned int  uint;
typedef unsigned long ulong;

typedef struct _PMProteinAtomDesc PMProteinAtomDesc;

typedef struct _PMProteinDensity{
	int                 type;
	size_t 		    frames;
	size_t 		    records;
	float 		    **data;
	size_t 		    cframeatt;
	float 		    **frameatt;
	size_t		    dim[3];
	float		    apix[3];
	float		    origin[3];
}PMProteinDensity;

typedef struct _PMProteinModel{
	int                 type;
	size_t 		    frames;
	size_t 		    records;
	float 		    **data;
	size_t 		    cframeatt;
	float 		    **frameatt;
	PMProteinAtomDesc   *desc;
	size_t 		    catomatt;
	float 		    **atomatt;
}PMProteinModel;

typedef struct _PMProteinAny{
	int                 type;
	size_t 		    frames;
	size_t 		    records;
	float 		    **data;
	size_t 		    cframeatt;
	float 		    **frameatt;
}PMProteinAny;

typedef union _PMProtein {
	int                 type;
	PMProteinAny        pany;
	PMProteinDensity    pdensity;
	PMProteinModel      pmodel;
} PMProtein;


PMProtein *pmCreateNew  (int type);
PMProtein *pmCreatefRef (PMProtein *ref);
PMProtein *pmCreatefReff(PMProtein *ref,size_t frames);
PMProtein *pmDelete     (PMProtein *protein);
PMProtein *pmAllocFrames(PMProtein *protein,size_t frames);
PMProtein *pmOpen       (PMProtein *protein,const char *fname);
PMProtein *pmOpenRef    (PMProtein *protein,const char *fname);
PMProtein *pmSave       (PMProtein *protein,const char *fname);
PMProtein *pmRmLast     (PMProtein *protein);

//char *pmInfo            (PMProtein *protein,char *buffer,size_t n);
size_t    pmGetDF       (PMProtein *Protein);

void pmPrintInfo        (PMProtein *ref);
void pmPrintInfoM       (PMProteinModel* ref);
void pmPrintInfoD       (PMProteinDensity* ref);

PMProtein *pmAddFrame   (PMProtein *protein,float *frame,size_t size);
PMProtein *pmAddFrames  (PMProtein *protein,float **frames,size_t cframes,size_t size);

PMProtein *pmReadP      (PMProtein        *protein,const char* fname);
PMProtein *pmReadD      (PMProteinDensity *protein,const char* fname);
PMProtein *pmReadM      (PMProteinModel   *protein,const char* fname);
PMProtein *pmReadRefM   (PMProteinModel   *protein,const char* fname);

PMProtein *pmWriteP  (PMProtein        *protein,const char* fname);
PMProtein *pmWriteD  (PMProteinDensity *protein,const char* fname);
PMProtein *pmWriteSD (PMProteinDensity *protein,const char* fname,size_t idx);
PMProtein *pmWriteM  (PMProteinModel   *protein,const char* fname);
PMProtein *pmWriteSM (PMProteinModel   *protein,const char* fname,size_t idx);

PMProtein *pmAdd       (PMProtein *out,PMProtein *a,PMProtein *b);
PMProtein *pmSub       (PMProtein *out,PMProtein *a,PMProtein *b);
PMProtein *pmMult      (PMProtein *out,PMProtein *a,PMProtein *b);
PMProtein *pmDiv       (PMProtein *out,PMProtein *a,PMProtein *b);
PMProtein *pmScale     (PMProtein *out,PMProtein *a,float s);
PMProtein *pmSqrt      (PMProtein *out,PMProtein *a);
PMProtein *pmIsqrt     (PMProtein *out,PMProtein *a);
PMProtein *pmAddNoise  (PMProtein *out,PMProtein *a,float sigma);

PMProtein *pmMean      (PMProtein *out,PMProtein *a);
PMProtein *pmMeanVar   (PMProtein *out,PMProtein *a);
float     *pmPCA       (PMProtein *out,PMProtein *a,PMProtein *mean);
PMProtein *pmLinTraj   (PMProtein *out,PMProtein *center,PMProtein *vector,size_t frame,float scale);

float     *pmCorrelation(PMProtein *out,PMProtein *a);
PMProtein *pmNormalize  (PMProteinDensity *out,PMProteinDensity *a);
PMProtein *pmRampFilter (PMProteinDensity *out,PMProteinDensity *in);
PMProtein *pmMask       (PMProteinDensity *out,PMProteinDensity *in, PMProteinDensity *mask, size_t frame, float thr);
PMProtein *pmThreshold  (PMProteinDensity *out,PMProteinDensity *in, float thr);

#ifdef USE_FFTW3 
PMProteinDensity *pmFFT               (PMProteinDensity *out,int flags);
PMProteinDensity *pmFFTShift          (PMProteinDensity *out,float dx,float dy,float dz);
PMProteinDensity *pmFFTCutOffFilter   (PMProteinDensity *out,float lb,float ub);
PMProteinDensity *pmFFTGaussFilter    (PMProteinDensity *out,float sigma);
PMProteinDensity *pmFFTLaplaceFilter  (PMProteinDensity *out);
PMProteinDensity *pmFFTRampFilter     (PMProteinDensity *out);
#endif


float *pmRMSD          (PMProteinModel *out,PMProteinModel *refi,size_t idx);
vec3  pmCOG            (PMProteinModel *out,PMProteinModel *ref);

PMProtein *pmTranslate (PMProteinModel *out,PMProteinModel *in,vec3 v);
PMProtein *pmRotate    (PMProteinModel *out,PMProteinModel *in,mat3 m);

PMProtein *pmCenter    (PMProteinModel *out,PMProteinModel *in);
PMProtein *pmAlignTo   (PMProteinModel *out,PMProteinModel *in,PMProteinModel *ref);
#ifdef __cplusplus
}
#endif
#endif
