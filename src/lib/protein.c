/***

    ./src/lib/protein.c 

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
#include "calc.h"
#include <stdlib.h>
#include <string.h>


extern PMProteinDensity *pmInitD_ (PMProteinDensity *protein,PMProteinDensity *ref);
extern PMProteinDensity *pmDelD_ (PMProteinDensity *protein);
extern PMProteinModel   *pmInitM_ (PMProteinModel *protein,PMProteinModel *ref);
extern PMProteinModel   *pmDelM_ (PMProteinModel *protein);

void pmPrintInfo(PMProtein *ref)
{
	switch(ref->type){
	case PM_TYPE_MODEL:
		pmPrintInfoM((PMProteinModel*)ref);
		break;
	case PM_TYPE_DENSITY:
	case PM_TYPE_DENSITY&PM_DENSITY_COMPLEX:
		pmPrintInfoD((PMProteinDensity*)ref);
		break;
	default:
		SAY("Protein Info:")
		SAY("Type: ???")
		SAY("Frames: %lu",ref->pany.frames)
		SAY("Records: %lu",ref->pany.records)
	}

}

PMProtein *pmCreateNew  (int type)
{
	PMProtein *p;
	p = (PMProtein *) calloc(1,sizeof(PMProtein));
	if(p)p->type=type;
	return p;
}

PMProtein *pmCreatefRef (PMProtein *ref)
{
	PMProtein *p; 
	p = (PMProtein *) calloc(1,sizeof(PMProtein));
	memcpy(p,ref,sizeof(PMProtein));

	p->pany.frames = 0;
	p->pany.data   = 0; 

	//type special
	switch(p->type){
		case PM_TYPE_MODEL:
			pmInitM_((PMProteinModel*)p,(PMProteinModel*)ref);
		case PM_TYPE_DENSITY:
			pmInitD_((PMProteinDensity*)p,(PMProteinDensity*)ref);
		break;
	}
	return p;
}

PMProtein *pmCreatefReff (PMProtein *ref,size_t frames)
{
	size_t i;
	PMProtein *p; 
	p = (PMProtein *) calloc(1,sizeof(PMProtein));
	memcpy(p,ref,sizeof(PMProtein));
        
	p->pany.frames = frames;

	p->pany.data = calloc(p->pany.frames,sizeof(float*));
	for(i=0;i<frames;i++)
		p->pany.data[i] = (float*)malloc(p->pany.records*sizeof(float));

	//type special
	switch(p->type){
		case PM_TYPE_MODEL:
			pmInitM_((PMProteinModel*)p,(PMProteinModel*)ref);
		case PM_TYPE_DENSITY:
			pmInitD_((PMProteinDensity*)p,(PMProteinDensity*)ref);
		break;
	}
	return p;
}

PMProtein *pmAllocFrames(PMProtein *protein,size_t frames)
{	
	size_t i;
	float **h;
	h = (float**) realloc(protein->pany.data,frames*sizeof(float*));
	if(!h)
		return 0;
	protein->pany.data=h;
	for(i=protein->pany.frames;i<frames;i++)
		protein->pany.data[i] = (float*)malloc(protein->pany.records*sizeof(float));
	return protein;
}

PMProtein *pmDelete     (PMProtein *p)
{
	size_t i;
	float **pp;
	//type special
	switch(p->type){
		case PM_TYPE_MODEL:
			pmDelM_((PMProteinModel*)p);
		case PM_TYPE_DENSITY:
			pmDelD_((PMProteinDensity*)p);
		break;
	}
	pp=p->pany.data;
	if(pp){
		for(i=p->pany.frames;i>0;i--)
			if(*pp) free(*pp++);
		free(p->pany.data);
	}
	free(p);
	return 0;
}

PMProtein *pmAddFrame  (PMProtein *p,float *frame,size_t size)
{
	float **h;
	if(p->pany.records==size){
		h = (float**) realloc(p->pany.data,(p->pany.frames+1)*sizeof(float*));
		if(!h){
			WARN(" [malloc failed]");
			return 0;
		}
		h[p->pany.frames]=frame;
		p->pany.frames++;
		p->pany.data=h;
		return p;
	}
	else if (p->pany.records==0){
		p->pany.records=size;
		h = (float**) realloc(p->pany.data,(p->pany.frames+1)*sizeof(float*));
		if(!h){
			WARN(" [malloc failed]");
			return 0;
		}
		h[p->pany.frames]=frame;
		p->pany.frames++;
		p->pany.data=h;
		return p;
	}
	WARN(" [do not match](%lu/%lu)",p->pany.records,size);
	return 0;
}

PMProtein *pmAddFrames (PMProtein *p,float **frames,size_t cframes,size_t size)
{
	float **h;
	size_t i;
	if(p->pany.records==size){
		h = (float**) realloc(p->pany.data,(p->pany.frames+cframes)*sizeof(float*));
		if(!h)
			WARN(" [malloc failed]");
			return 0;
		p->pany.data=h;
		h=&h[p->pany.frames];
		for(i=cframes;i>0;i++){
			(*h++)=(*frames++);	
		}
		p->pany.frames+=cframes;
		return p;
	}
	else if (p->pany.records==0){
		p->pany.records=size;
		h = (float**) realloc(p->pany.data,(p->pany.frames+cframes)*sizeof(float*));
		if(!h)
			WARN(" [malloc failed]");
			return 0;
		p->pany.data=h;
		h=&h[p->pany.frames];
		for(i=cframes;i>0;i++){
			(*h++)=(*frames++);	
		}
		p->pany.frames+=cframes;
		return p;
	}
	WARN(" [do not match](%lu/%lu)",p->pany.records,size);
	return 0;
}



PMProtein *pmAdd       (PMProtein *out,PMProtein *a,PMProtein *b)
{
	addArray2D(out->pany.data,a->pany.data,b->pany.data,out->pany.frames,out->pany.records);
	return out;
}

PMProtein *pmSub       (PMProtein *out,PMProtein *a,PMProtein *b)
{
	subArray2D(out->pany.data,a->pany.data,b->pany.data,out->pany.frames,out->pany.records);
	return out;
}

PMProtein *pmMult      (PMProtein *out,PMProtein *a,PMProtein *b)
{
	multArray2D(out->pany.data,a->pany.data,b->pany.data,out->pany.frames,out->pany.records);
	return out;
}

PMProtein *pmDiv       (PMProtein *out,PMProtein *a,PMProtein *b)
{
	divArray2D(out->pany.data,a->pany.data,b->pany.data,out->pany.frames,out->pany.records);
	return out;
}

PMProtein *pmScale     (PMProtein *out,PMProtein *a,float s)
{
	scaleArray2D(out->pany.data,a->pany.data,s,out->pany.frames,out->pany.records);
	return out;
}

PMProtein *pmSqrt      (PMProtein *out,PMProtein *a)
{
	sqrtArray2D(out->pany.data,a->pany.data,out->pany.frames,out->pany.records);
	return out;
}

PMProtein *pmIsqrt     (PMProtein *out,PMProtein *a)
{
	invsqrtArray2D(out->pany.data,a->pany.data,out->pany.frames,out->pany.records);
	return out;
}

PMProtein *pmAddNoise  (PMProtein *out,PMProtein *a,float sigma)
{
	addNoiseArray2D(out->pany.data,a->pany.data,sigma,out->pany.frames,out->pany.records);
	return out;
}

PMProtein *pmMean      (PMProtein *out,PMProtein *a)
{
	meanArray2D(a->pany.data, out->pany.data[0], a->pany.frames,out->pany.records);
	return out;
}

PMProtein *pmMeanVar   (PMProtein *out,PMProtein *a)
{
	meanVarArray2D(a->pany.data,out->pany.data[0],out->pany.data[1],a->pany.frames,out->pany.records);
	return out;
}

size_t pmGetDF(PMProtein *out)
{
	if(out->pany.frames<out->pany.records){
		SAY("DOF: %zu",out->pany.frames)
		return out->pany.frames; 
	}
	else{ 
		SAY("DOF: %zu",out->pany.records)
		return out->pany.records;
	}
}

float *pmPCA       (PMProtein *out,PMProtein *a, PMProtein *mean)
{
	float *m,*e,*cov,h;
	size_t i;
	if(!mean){
		m=(float*) malloc(out->pany.records*sizeof(float));
		meanArray2D(a->pany.data, m,a->pany.frames,out->pany.records);
	}	
	else 
		m=mean->pany.data[0];

	if(a->pany.frames<out->pany.records){//backward
		cov = inverseCovarianceArray2D(a->pany.data,m,a->pany.frames,out->pany.records);
		e = (float*)malloc(a->pany.frames*sizeof(float));
		calcEigenVectors(cov,e,a->pany.frames);
		//scale eigenvalues 
		h=1./a->pany.frames; 
		for(i=0;i<out->pany.frames;i++){
			e[i]*=h;
		}
		//backtransform and cpy 
		for(i=0;i<out->pany.frames;i++){
			covarianceBTArray(out->pany.data[i],cov,m,a->pany.data, i, out->pany.frames,out->pany.records);
	   	}
	}
	else {
		cov = covarianceArray2D(a->pany.data,m,a->pany.frames,out->pany.records);
		e = (float*)malloc(out->pany.records*sizeof(float));
		calcEigenVectors(cov,e,out->pany.records);
		h=1./a->pany.frames;
		for(i=0;i<out->pany.records;i++)
			e[i]*=h;  
		for(i=0;i<out->pany.frames;i++) 
			memcpy(out->pany.data[i],cov,out->pany.records*sizeof(float));
	}
	if(!mean) free(m); 
	free(cov);
	return e;
}

PMProtein *pmRmLast     (PMProtein *protein)
{
	protein->pany.frames--;
	return protein;
}

PMProtein *pmLinTraj   (PMProtein *out,PMProtein *center,PMProtein *vector,size_t frame,float scale)
{
	trajectoryArray2D(out->pany.data, center->pany.data[0],vector->pany.data[frame], scale, out->pany.frames,out->pany.records);
	return out;
}
