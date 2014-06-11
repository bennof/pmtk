#ifdef __cplusplus
extern "C"
#endif

#include "pmtk.h"
#include "config.h"
#include <stdlib.h>


/*
 * function from xdrfile_xtc.h,v 1.5
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Coordinate compression (c) by Frans van Hoesel.
 * LGPL3
 *
 */
 typedef struct _XDRFILE XDRFILE;
 XDRFILE* xdrfile_open     (const char* filename,const char* mode);
 void xdrfile_close        (XDRFILE *file);
 int read_xtc_natoms       (const char *fn,int *natoms);
 int read_xtc              (XDRFILE *xd,int natoms,int *step,float *time,float *box,float *x,float *prec);
 int write_xtc             (XDRFILE *xd,int natoms,int step,float time,float *box,float *x,float prec);


PMProteinModel *xtcopen(PMProteinModel *protein,const char * filename, int mode)
{
	int N,S,e,i;
	float t,H;
	float *buffer,*b;
	float h[9]={1,0,0,0,1,0,0,0,1};
	size_t size;
	XDRFILE *f;
	PMProteinModel *p=protein;

	SAY("Open: %s",filename);
	read_xtc_natoms((char*)filename,&N);
	size=N;

	if(!p)
		p=(PMProteinModel*)pmCreateNew(PM_TYPE_MODEL);
	if(p->type!=PM_TYPE_MODEL){
		WARN("Mixing different representations (abort): %s",filename);
		if(p!=protein) pmDelete((PMProtein*)p);
		return 0;
	}
	if(p->records==0){
		p->records=3*size;
	}
	else {
		if(protein->records!=3*size){
			WARN("Model differs: %s",filename);
			if(p!=protein) pmDelete((PMProtein*)p);
			return 0; 
		}
	}

	if(mode==1){
		return p;
	}

	if(!(f=xdrfile_open(filename,(const char*)"r"))){
		WARN("Failed open: %s",filename);
		if(p!=protein) pmDelete((PMProtein*)p);
		return 0;
	}

	//read loop
	do {
		buffer=(float*)malloc(p->records*sizeof(float));
	        e=read_xtc(f,N,&S,&t,h,buffer,&H);
		if(e) {
			free(buffer);
			break;
		}
		b=buffer;
		for(i=p->records;i>0;i--)
			(*b++)*=10.0f;

		if(!pmAddFrame((PMProtein*)protein,buffer,protein->records)){
			if(p!=protein) pmDelete((PMProtein*)p);
			xdrfile_close(f);
			return 0;
		}

	} while(1);
	if(e!=11) {
		WARN("Failed reading [code %i]:%s",e,filename);	
		if(p!=protein) pmDelete((PMProtein*)p);
		xdrfile_close(f);
		return 0;
	}
	xdrfile_close(f);
	return p;
}



PMProteinModel *xtcsave(PMProteinModel *protein,const char * filename, int mode)
{
	XDRFILE *f;
	float **H,*b;
	float h[9]={1,0,0,0,1,0,0,0,1};
	size_t i,j;
	int e;

	SAY("Open: %s",filename);

	if(protein->type!=PM_TYPE_MODEL){
		WARN("Mixing different representations (aborted): %s",filename);
		return 0;
	}
	H=(protein->data);
	if((f=xdrfile_open(filename,(const char*)"w"))){
		for(i=0;i<protein->frames;i++){
	        	printf("[%zu]\r",i);
			b=*H;
			for(j=protein->records;j>0;j--)
			(*b++)*=0.1f;
			e=write_xtc(f,protein->records/3,i,0.0,h,*H,1000.0);
			b=*H;
			for(j=protein->records;j>0;j--)
				(*b++)*=10.0f;
			H++;
		}
	}
	else {
		WARN("Failed open: %s",filename);
		return 0;
	}
	xdrfile_close(f);
	return protein;
}
