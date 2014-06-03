/***

    ./src/lib/model.c 

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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pdb.h"
#include "calc.h"



struct _PMProteinAtomDesc{
	size_t refcount;
	size_t natoms;
	PMProteinAtomDesc *next;
	PMAtomDesc *atoms;
};

PMProteinAtomDesc *PM_INFO_STACK=0;

static PMProteinAtomDesc *new_info(size_t natoms, PMAtomDesc *atoms)
{
	PMProteinAtomDesc* i;
	i = (PMProteinAtomDesc*)malloc(sizeof(PMProteinAtomDesc));
	i->refcount=1;
	i->natoms=natoms;
	i->atoms=atoms;
	i->next=PM_INFO_STACK;
	PM_INFO_STACK=i->next;
	return i;
}
static PMProteinAtomDesc *del_info(PMProteinModel *protein)
{
	PMProteinAtomDesc *i,*p;
	p = protein->desc;

	if(!p)
		return 0;

	p->refcount--;
	if(!p->refcount){
		//remove info
		i=PM_INFO_STACK;
		if(i==p){
			PM_INFO_STACK = i->next;
		}
		else{
			while(i && (i->next != p))
				i=i->next;
			if(i)
				i->next=p->next;
		}
		free(p->atoms);
		free(p);
	}
	return 0;
}

static PMProteinAtomDesc *dub_info(PMProteinModel *protein)
{
	if(protein->desc)
		protein->desc->refcount++;
	return protein->desc;
}


void pmPrintInfoM(PMProteinModel* ref)
{
	SAY("Protein Info:")
	SAY("Type: Model")
  	SAY("Atoms: %lu",ref->records/3);
	SAY("PDB Info: %s",(ref->desc)?"yes":"no")
	SAY("Frames: %lu",ref->frames)
	SAY("Records: %lu",ref->records)
}

PMProteinModel *pmInitM_ (PMProteinModel *protein,PMProteinModel *ref)
{	
	protein->desc = dub_info(ref);
	return protein;
}

PMProteinModel *pmDelM_ (PMProteinModel *protein)
{
	protein->desc = del_info(protein);
	return protein;
}

PMProtein *pmReadM   (PMProteinModel   *protein,const char* fname)
{
	XDRFILE *f;
	int N,S,e,i,j=0;
	float t,p;
	float *buffer,*b;
	float h[9]={1,0,0,0,1,0,0,0,1};
	FILE *F;
	PMAtomDesc *info_r;
	size_t iatoms = protein->records/3;

	switch(getType(fname)){
	case PM_ATOM_FILE_XDR:
		if((f=xdrfile_open(fname,(const char*)"r"))){
		    read_xtc_natoms((char*)fname,&N);
		    if(!protein->records){
		    	protein->records=3*N;
		    }
		    if(protein->records==3*N){
		    	buffer=(float*)malloc(protein->records*sizeof(float));
		    	while(!(e=read_xtc(f,N,&S,&t,h,buffer,&p))){
	        		printf("[%i]\r",j++);
		    		b=buffer;
		    		for(i=protein->records;i>0;i--)
		    			(*b++)*=10.0f;
		    		if(!pmAddFrame((PMProtein*)protein,buffer,protein->records)){
					WARN(" [malloc failed]");
		    			return 0;
				}
		        	buffer=(float*)malloc(protein->records*sizeof(float));
		    	}
		    	free(buffer);
		    	if(e && e!=11) {
				WARN(" [XDR returned %i!=11]",e);	
				return 0;
			}
		    }
		    else {
			WARN(" [failed number of atoms]");
		    	return 0;
		    }
		    xdrfile_close(f);
		  }
		else{
			WARN(" [failed open file %s]",fname);
			return 0;
		}
		return (PMProtein*) protein;
		break;
	case PM_ATOM_FILE_PDB:
		F = fopen(fname,"r");
		if(!protein->desc || protein->desc->atoms){
			buffer = pmFReadPDBFrameFull(F,&info_r,&iatoms);
			if(buffer){
				pmAddFrame((PMProtein*)protein,buffer,protein->records);
				protein->desc = new_info(iatoms,info_r);
				protein->records=iatoms*3;
			}
		}
		while(!feof(F)){
			buffer = pmFReadPDBFramePlain(F,&protein->desc->natoms);
			if(buffer)
				pmAddFrame((PMProtein*)protein,buffer,protein->records);
		}

		fclose(F);
		break;
	default:
	        WARN("[Unkown file format: %s]",fname);
		return 0;
		break;
	}
	return (PMProtein*)protein;
}

PMProtein *pmReadRefM   (PMProteinModel   *protein,const char* fname)
{
	FILE *F;
	float *buffer;
	PMAtomDesc *info_r;
	size_t iatoms = protein->records/3;

	switch(getType(fname)){
	case PM_ATOM_FILE_PDB:
		F = fopen(fname,"r");
		if(!protein->desc || protein->desc->atoms){
			buffer = pmFReadPDBFrameFull(F,&info_r,&iatoms);
			if(buffer){
				free(buffer);
				protein->desc = new_info(iatoms,info_r);
				protein->records=iatoms*3;
			}
		}
		fclose(F);
		break;
	default:
	        WARN("[Unkown file format: %s]",fname);
		return 0;
		break;
	}
	return (PMProtein*)protein;
}

PMProtein *pmWriteM  (PMProteinModel   *protein,const char* fname)
{
	XDRFILE *f;
	float **H,*b;
	float h[9]={1,0,0,0,1,0,0,0,1};
	size_t i,j;
	int e;
	FILE *F;

	switch(getType(fname)){
	case PM_ATOM_FILE_XDR:
		H=(protein->data);
		if((f=xdrfile_open(fname,(const char*)"w"))){
			for(i=0;i<protein->frames;i++){
	        		printf("[%zu]\r",i);
				b=*H;
				for(j=protein->records;j>0;j--)
					(*b++)*=0.1f;
				e=write_xtc(f,protein->records/3,i,0.0,h,*H,1000.0);
				b=*H;
				for(j=protein->records;j>0;j--)
					(*b++)*=1.0f;
				H++;
			}
			xdrfile_close(f);
		}
		else {
			fputs(" [failed open]\n",stdout);
			return 0;
		}
		return (PMProtein*)protein;
		break;
	case PM_ATOM_FILE_PDB:
            	if((F = fopen(fname,"w"))){
	    		pmFWriteFramesPDB(F,(protein->desc)?protein->desc->atoms:0,protein->data,protein->frames,protein->records/3);
	    		fclose(F);
		}
		else{
			fputs(" [failed open]\n",stdout);
			return 0;
		}
	    	return (PMProtein*)protein;
		break;
	default:
		return 0;
		break;
	}
	return (PMProtein*)protein;

}

PMProtein *pmWriteSM  (PMProteinModel   *protein,const char* fname,size_t idx)
{
	XDRFILE *f;
	float *H,*b;
	float h[9]={1,0,0,0,1,0,0,0,1};
	unsigned int i;
	int e;
	FILE *F;

	switch(getType(fname)){
	case PM_ATOM_FILE_XDR:
	        fprintf(stdout,"XTC write: %s\r",fname);
		H=(protein->data[idx]);
		if((f=xdrfile_open(fname,(const char*)"w"))){
			b=H;
			for(i=protein->records;i>0;i--)
				(*b++)*=0.1f;
			e=write_xtc(f,protein->records/3,i,0.0,h,H,1000.0);
			b=H;
			for(i=protein->records;i>0;i--)
				(*b++)*=10.0f;
			xdrfile_close(f);
		}
		else {
			fputs("\n",stdout);
			return 0;
		}
		return (PMProtein*)protein;
		break;
	case PM_ATOM_FILE_PDB:
	    	fprintf(stdout,"PDB write: %s\r",fname);
            	if((F = fopen(fname,"w"))){
	    		pmFWriteFramesPDB(F,protein->desc->atoms,&protein->data[idx],1,protein->records/3);
	    		fclose(F);
		}
		else{
			fputs(" [failed open]\n",stdout);
			return 0;
		}
		fputs(" [ok]\n",stdout);
	    	return (PMProtein*)protein;
		break;
	default:
		return 0;
		break;
	}
	return (PMProtein*)protein;

}

float *pmRMSD        (PMProteinModel *out,PMProteinModel *ref,size_t idx)
{
	size_t i;
	float *h,*r = malloc(out->frames*sizeof(float));
	h=r;
	for(i=0;i<out->frames;i++){
		*h = rmsdArray(out->data[i], ref->data[idx], out->records, 0);
		h++;
	}
	return r;	
}

vec3  pmCOG         (PMProteinModel *out,PMProteinModel *ref)
{
	return meanVec3Array2D(out->data,out->frames,out->records);
}

PMProtein *pmTranslate   (PMProteinModel *out,PMProteinModel *in,vec3 v)
{
	translateVec3Array2D(out->data,in->data,v,out->frames,out->records);
	return (PMProtein*) out;
}

PMProtein *pmRotate      (PMProteinModel *out,PMProteinModel *in,mat3 m)
{
	rotateVec3Array2D(out->data,in->data,m,out->frames,out->records);
	return (PMProtein*) out;
}

PMProtein *pmCenter      (PMProteinModel *out,PMProteinModel *in)
{
	centerAtVec3Array2D(out->data,in->data,VEC3(0.,0.,0.),out->frames,out->records);
	return (PMProtein*) out;
}

PMProtein *pmAlignTo     (PMProteinModel *out,PMProteinModel *in,PMProteinModel *ref)
{
	alignAtVec3Array2D(out->data,in->data,ref->data[0],out->frames,out->records);
	return (PMProtein*)out;
}
