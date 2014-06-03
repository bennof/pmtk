/***

    ./src/lib/density.c 

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



#ifdef __cplusplus
extern "C" 
#endif

#include "pmtk.h"
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "calc.h"
#include "mrc.h"

PMProteinDensity *pmInitD_ (PMProteinDensity *protein,PMProteinDensity *ref)
{
	return protein;
}
PMProteinDensity *pmDelD_ (PMProteinDensity *protein)
{
	return protein;
}

void pmPrintInfoD(PMProteinDensity* ref)
{
	SAY("Protein Info:")
	SAY("Type: Density") 
	if (ref->type&PM_DENSITY_COMPLEX)
		SAY("Space: phase space")
	else
		SAY("Space: real")
	SAY("Dimensions: %lux%lux%lu",ref->dim[0],ref->dim[1],ref->dim[2])
	SAY("Voxelsize in A/pixel: (%f,%f,%f)",ref->apix[0],ref->apix[1],ref->apix[2])
	SAY("Origin: (%f,%f,%f)",ref->origin[0],ref->origin[1],ref->origin[2])
}


PMProtein *pmReadD   (PMProteinDensity *protein,const char* fname)
{
	PMMRCHeader head;
	float *h;
	fprintf(stdout,"MRC read: %s",fname);
	if(protein->records==0){
		if(!(h=pmFReadMRC(fname,&head,0,0))){
			fputs(" [failed]\n",stdout);
			return 0;
		}
		protein->dim[0]=head.nx;
		protein->dim[1]=head.ny;
		protein->dim[2]=head.nz;
		protein->apix[0]=head.cella[0]/head.mx;
		protein->apix[1]=head.cella[1]/head.my;
		protein->apix[2]=head.cella[2]/head.mz;
		protein->origin[0]=head.origin[0];
		protein->origin[1]=head.origin[1];
		protein->origin[2]=head.origin[2];
		protein->records = protein->dim[0]*protein->dim[1]*protein->dim[2];
	}
	else {
		if(!(h=pmFReadMRC(fname,&head,protein->records*sizeof(float),0))){
			fputs(" [failed]\n",stdout);
			return 0;
		}
	}
	pmAddFrame((PMProtein*)protein,h,protein->records);
	fputs(" [ok]\n",stdout);
	return (PMProtein*)protein;
}

#ifndef STRING_BUFFER_SIZE
#define STRING_BUFFER_SIZE 2048
#endif
PMProtein *pmWriteD  (PMProteinDensity *protein,const char* fname)
{
	PMMRCHeader head;
	char buffer[STRING_BUFFER_SIZE];
	char *e,*p;
	size_t l,le,d,c;
	if(protein->frames==1){
		fprintf(stdout,"MRC write: %s\r",fname);
		setMRCHeader(&head,protein->dim,protein->apix,protein->origin);
		c=pmFWriteMRC(fname,&head, protein->data[0],protein->records*sizeof(float),0);
		if(c){
			fputs("\n",stdout);
			return 0;
		}
	}
	else {
		e=(char*)fname;
		while(e){
			if(*e=='.')p=e;
			e++;
		}
		if(p==0)
		 	return 0;
		l = (size_t) p-(size_t)fname;
		le = (size_t) e-(size_t)p;
		if(l+le+10 > STRING_BUFFER_SIZE)
			return 0;
	
		memcpy(buffer,fname,l);
		e=&buffer[l];
		d=0;
		size_t i;
		setMRCHeader(&head,protein->dim,protein->apix,protein->origin);
		for(i=0;i<protein->frames;i++){
			sprintf(&buffer[l],"%lu%s",i,p);
			fprintf(stdout,"MRC write: %s",fname);
			c=pmFWriteMRC(fname,&head,protein->data[i],protein->records*sizeof(float),0);
			if(c){
				fputs("\n",stdout);
				return 0;
			}
		}
	}
	return (PMProtein*)protein;
}
PMProtein *pmWriteSD  (PMProteinDensity *protein,const char* fname,size_t idx)
{
	PMMRCHeader head;
	size_t c;
	fprintf(stdout,"MRC write: %s\r",fname);
	setMRCHeader(&head,protein->dim,protein->apix,protein->origin);
	c=pmFWriteMRC(fname,&head, protein->data[idx],protein->records*sizeof(float),0);
	if(c)
		fputs("\n",stdout);
		return 0;
	return (PMProtein*)protein;
}


PMProtein *pmNormalize  (PMProteinDensity *out,PMProteinDensity *a)
{
	normalizeArray2D(out->data,a->data,out->frames,out->records); 
	return (PMProtein*)out;
}

PMProtein *pmRampFilter (PMProteinDensity *out,PMProteinDensity *in)
{
	rampArray2D(out->data,in->data,*(uvec3*)&in->dim,*(vec3*)&in->apix,out->frames,out->records);
	return (PMProtein*)out;
}
PMProtein *pmMask       (PMProteinDensity *out,PMProteinDensity *in, PMProteinDensity *mask, size_t frame, float thr)
{
	maskArray2D(out->data,in->data,*mask->data,thr,out->frames,out->records);
	return (PMProtein*)out;
}
PMProtein *pmThreshold  (PMProteinDensity *out,PMProteinDensity *in, float thr)
{
	thresholdArray2D(out->data,in->data,thr,out->frames,out->records);
	return (PMProtein*)out;
}

float *pmCorrelation(PMProtein *a,PMProtein *b)
{
	size_t i;
	float *h,*r = malloc(a->pany.frames*sizeof(float));
	h=r;
	for(i=0;i<a->pany.frames;i++){
		correlationArray(h,a->pany.data[i],b->pany.data[i],a->pany.records);
		h++;
	}
	return r;
}

#ifdef USE_FFTW3 
PMProteinDensity *pmFFT               (PMProteinDensity *out,int flags)
{
	if(PM_FFT_FORWARD&flags){
		dftr2cArray2D(out->data,*(uvec3*)&out->dim,&out->records,out->frames);
		out->type&=PM_DENSITY_COMPLEX;
	}
	else{
		dftc2rArray2D(out->data,*(uvec3*)&out->dim,&out->records,out->frames);
		out->type^=PM_DENSITY_COMPLEX;
	}
	return out;
}

PMProteinDensity *pmFFTShift          (PMProteinDensity *out,float dx,float dy,float dz)
{
	dftshiftArray2D(out->data,*(uvec3*)out->dim, dx, dy, dx,out->frames);
	return out;
}
PMProteinDensity *pmFFTCutOffFilter   (PMProteinDensity *out,float lb,float ub)
{
	dftFilterCutOffArray2D(out->data,*(uvec3*)out->dim,lb,ub,out->frames);
	return out;
}
PMProteinDensity *pmFFTGaussFilter    (PMProteinDensity *out,float sigma)
{
	dftFilterGaussArray2D(out->data,*(uvec3*)out->dim,*(vec3*)out->apix,sigma,out->frames);
	return out;
}
PMProteinDensity *pmFFTLaplaceFilter  (PMProteinDensity *out)
{
	dftFilterLaplaceArray2D(out->data,*(uvec3*)out->dim,out->frames);
	return out;
}
PMProteinDensity *pmFFTRampFilter     (PMProteinDensity *out)
{
	dftFilterRampArray2D(out->data,*(uvec3*)out->dim,*(vec3*)out->apix,out->frames);
	return out;
}
#endif
