#ifdef __cplusplus
extern "C"
#endif

#include "pmtk.h"
#include "config.h"
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


/**
 * Modes of a mrc file data section
 * @{
 */
#define CEM_MRC_MODE_CHAR				0 /**< image : signed 8-bit bytes range -128 to 127 */
#define CEM_MRC_MODE_SHORT				1 /**< image : 16-bit halfwords */
#define CEM_MRC_MODE_FLOAT				2 /**< image : 32-bit reals */
#define CEM_MRC_MODE_COMPLEX_INT16	   		3 /**< transform : complex 16-bit integers */
#define CEM_MRC_MODE_COMPLEX_FLOAT			4 /**< transform : complex 32-bit reals */
#define CEM_MRC_MODE_UNKOWN				5 /**< not used */
#define CEM_MRC_MODE_USHORT				6 /**< image : unsigned 16-bit range 0 to 65535 */
/** @} */

/**
 * MRC header structure
 */
typedef struct _PMMRCheader {
	int nx;					/**< number of columns (fastest changing in map) */
	int ny; 				/**< number of rows */
	int nz;					/**< number of sections (slowest changing in map) */
	int mode;				/**< data type @see MRC Modes */
	int nxstart;			/**< number of first cols in map (Default = 0) */
	int nystart;			/**< number of first rows in map */
	int nzstart;			/**< number of first section in map */
	int mx;					/**< number of intervals along X */
	int my;					/**< number of intervals along Y */
	int mz;					/**< number of intervals along Z */
	float cella[3];			/**< cell dimensions in angstroms */
	float cellb[3];			/**< cell angles in degrees */
	int mapc;				/**< axis corresp to cols (1,2,3 for X,Y,Z) */
	int mapr;				/**< axis corresp to rows (1,2,3 for X,Y,Z) */
	int maps;				/**< axis corresp to sections (1,2,3 for X,Y,Z) */
	float min;				/**< minimum density value */
	float max;				/**< maximum density value */
	float mean;				/**< mean density value */
	int ispg;				/**< pace group number 0 or 1 (default=0) */
	int nsymbt;				/**< number of bytes used for symmetry data (0 or 80) */
	int extra[25];			/**< extra space used for anything   - 0 by default */
	float origin[3];		/**< origin in X,Y,Z used for transforms */
	char map[4];			/**< character string 'MAP ' to identify file type */
	int machst;				/**< machine stamp */
	int rms;				/**< rms deviation of map from mean density */
	int nlabl;				/**< number of labels being used */
	char label[10][80];		/**< 10 80-character text labels */
} PMMRCHeader;


static void setMRCHeader(PMMRCHeader *header,size_t dim[], float apix[], float origin[])
{
	header->nx=dim[0];
	header->ny=dim[1];
	header->nz=dim[2];
	header->mode=CEM_MRC_MODE_FLOAT;
	header->nxstart=0;
	header->nystart=0;
	header->nzstart=0;
	header->mx=dim[0];
	header->my=dim[1];
	header->mz=dim[2];
	header->cella[0]=apix[0]*dim[0];
	header->cella[1]=apix[1]*dim[1];
	header->cella[2]=apix[2]*dim[2];
	header->cellb[0]=90.0f;
	header->cellb[1]=90.0f;
	header->cellb[2]=90.0f;
	header->mapc=1;
	header->mapr=2;
	header->maps=3;
	header->min=0.0f;
	header->max=0.0f;
	header->mean=0.0f;
	header->ispg=0;
	header->nsymbt=0;
	memset(header->extra,0,4*25);
	header->origin[0]=origin[0];
	header->origin[1]=origin[1];
	header->origin[2]=origin[2];
	header->map[0]='M';
	header->map[1]='A';
	header->map[2]='P';
	header->map[3]=0;
	header->machst=0x12345678;
	header->rms=0.0f;
	header->nlabl=0;
}


PMProteinDensity *mrcopen(PMProteinDensity* protein,const char *filename,int mode)
{
	int fd;
	size_t c,i,l;
	float *r,*R;
	char *a;
	short *s;
	unsigned short *us;
	PMProteinDensity* p=protein;
	PMMRCHeader head;

	SAY("Open: %s",filename);
	fd = open(filename,O_RDONLY); //open file
	if(fd<0){ 
		WARN("Failed open: %s",filename);
		return 0;
	}

	c = read(fd,&head,1024); //read header
	if(c!=1024){
		WARN("Missing Header: %s",filename);
		close(fd);
		return 0;
	}
	if(!p)
		p=(PMProteinDensity*)pmCreateNew(PM_TYPE_DENSITY);

	if(!p->type!=PM_TYPE_DENSITY){
		WARN("Mixing different representations (aborted): %s",filename);
		if(p!=protein)pmDelete((PMProtein*)p); 
		return 0;
	}

	if(p->records==0){ 
		p->dim[0]=head.nx;
		p->dim[1]=head.ny;
		p->dim[2]=head.nz;
		p->apix[0]=head.cella[0]/head.mx;
		p->apix[1]=head.cella[1]/head.my;
		p->apix[2]=head.cella[2]/head.mz;
		p->origin[0]=head.origin[0];
		p->origin[1]=head.origin[1];
		p->origin[2]=head.origin[2];
		p->records = p->dim[0]*p->dim[1]*p->dim[2];
	}
	else {
		if(p->records!=head.nx*head.ny*head.nz){
			WARN("Density differs: %s",filename);
			if(p!=protein)pmDelete((PMProtein*)p); 
			return 0; 
		}
	}
	if(mode == 1)
		return p;

	r = malloc(p->records*sizeof(float));
	if(r==0){
		close(fd);
		WARN("Allocation failed: %s",filename);
		if(p!=protein)pmDelete((PMProtein*)p); 
		return 0;
	}

	switch(head.mode){
	case CEM_MRC_MODE_CHAR:
		l=p->records;
		c = read(fd,r,l);
		if(c!=l){
			free(r);
			r=0;
		}
		else {
		    a=(char*)r;
		    a=&a[l];
		    R=&r[l];
		    for(i=0;i<l;i++){
		    	*R = ((float)*a)/128.;
		    	R--;
		    	a--;
		    }
		}
		break;
	case CEM_MRC_MODE_SHORT:
		l=p->records*sizeof(short);
		c = read(fd,r,l);
		if(c!=l){
			free(r);
			r=0;
		}
		else {
		    s=(short*)r;
		    s=&s[l];
		    R=&r[l];
		    for(i=0;i<l;i++){
		    	*R = ((float)*s)/65536.;
		    	R--;
		    	s--;
		    }
		}
		break;
	case CEM_MRC_MODE_FLOAT:
		c = read(fd,r,p->records*sizeof(float));
		if(c!=p->records*sizeof(float)){
			free(r);
			r=0;
		}
		break;
	case CEM_MRC_MODE_COMPLEX_INT16:
	case CEM_MRC_MODE_COMPLEX_FLOAT:
		r=0;
		break;
	case CEM_MRC_MODE_USHORT:
		l=p->records*sizeof(short);
		c = read(fd,r,l);
		if(c!=l){
			free(r);
			r=0;
		}
		else {
		    us=(unsigned short*)r;
		    us=&us[l];
		    R=&r[l];
		    for(i=0;i<l;i++){
		    	*R = ((float)*us)/65536.;
		   		R--;
		   		us--;
		    }
		}
		break;
	default:
		r = 0;
		break;
	}
    	close(fd);
    	if(r)
    		return (PMProteinDensity*)pmAddFrame((PMProtein*)p,r,p->records);
    	else {
        	WARN("No density read: %s",filename);
		if(p!=protein)pmDelete((PMProtein*)p); 
    		return 0;
    	}
}


#ifndef STRING_BUFFER_SIZE
#define STRING_BUFFER_SIZE 2048
#endif

static float* write_(const char* filename,PMMRCHeader *head, float * data, size_t size)
{
	int fd;
	SAY("Save: %s",filename);
	fd = open(filename,O_WRONLY|O_CREAT|O_TRUNC,0644); //open file
	if(fd<0) {
		WARN("Failed open: %s",filename);
		return 0;
	}
	if(1024!=write(fd,(const void*)head,1024)){  //write header
		WARN("Header not written: %s",filename);
		close(fd);
		return 0;
	}
	if(size!=write(fd,data,size)){ //write data
		WARN("Data not written: %s",filename);
		close(fd);
		return 0;
	}
	close(fd);
	return data;
}


PMProteinDensity* mrcsave(PMProteinDensity* protein, const char* filename,int mode)
{
	PMMRCHeader head;
	char buffer[STRING_BUFFER_SIZE];
	char *e,*p;
	size_t l,le,d;
	size_t i;
	
	if(!protein->type!=PM_TYPE_DENSITY){
		WARN("Mixing different representations (aborted): %s",filename);
		return 0;
	}

	setMRCHeader(&head,protein->dim,protein->apix,protein->origin);

	if(protein->frames==1){
		if(write_(filename,&head,protein->data[0],protein->records*sizeof(float)))
			return protein;
		else{
			return 0;
		}
	}
	else {
		e=(char*)filename;
		while(e){
			if(*e=='.')p=e;
			e++;
		}
		if(p==0)
		 	return 0;
		l = (size_t) p-(size_t)filename;
		le = (size_t) e-(size_t)p;
		if(l+le+10 > STRING_BUFFER_SIZE)
			return 0;
	
		memcpy(buffer,filename,l);
		e=&buffer[l];
		d=0;

		for(i=0;i<protein->frames;i++){
			sprintf(&buffer[l],"%lu%s",i,p);
			if(!write_(filename,&head,protein->data[i],protein->records*sizeof(float))){
				return 0;	
			}
		}
		return protein;
	
	}
}


