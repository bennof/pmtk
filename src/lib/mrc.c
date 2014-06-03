/***

    ./src/lib/mrc.c 

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

#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "mrc.h"

/**
 * Modes of a mrc file data section
 * @{
 */
#define CEM_MRC_MODE_CHAR				0 /**< image : signed 8-bit bytes range -128 to 127 */
#define CEM_MRC_MODE_SHORT				1 /**< image : 16-bit halfwords */
#define CEM_MRC_MODE_FLOAT				2 /**< image : 32-bit reals */
#define CEM_MRC_MODE_COMPLEX_INT16	    3 /**< transform : complex 16-bit integers */
#define CEM_MRC_MODE_COMPLEX_FLOAT		4 /**< transform : complex 32-bit reals */
#define CEM_MRC_MODE_UNKOWN				5 /**< not used */
#define CEM_MRC_MODE_USHORT				6 /**< image : unsigned 16-bit range 0 to 65535 */
/** @} */

void setMRCHeader(PMMRCHeader *header,size_t dim[], float apix[], float origin[])
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


float *pmFReadMRC(const char *filename, PMMRCHeader* header,size_t size,int mode)
{
	int fd;
	size_t c,i,l;
	float *r,*p;
	char *a;
	short *s;
	unsigned short *us;

	fd = open(filename,O_RDONLY); //open file
	if(fd<0) return 0;

	c = read(fd,header,1024); //read header
	if(c!=1024){
		close(fd);
		return 0;
	}
	if(size){
		if(size!=header->nx*header->ny*header->nz*sizeof(float))
			return 0;
	}
	else{
		size = header->nx*header->ny*header->nz*sizeof(float);
	}

	r = malloc(size);
	if(r==0){
		close(fd);
		return 0;
	}

	switch(header->mode){
	case CEM_MRC_MODE_CHAR:
		l=size/sizeof(float);
		c = read(fd,r,l);
		if(c!=l){
			free(r);
			r=0;
		}
		else {
		    a=(char*)r;
		    a=&a[l];
		    p=&r[l];
		    for(i=0;i<l;i++){
		    	*p = ((float)*a)/128.;
		    	p--;
		    	a--;
		    }
		}
		break;
	case CEM_MRC_MODE_SHORT:
		l=size/sizeof(short);
		c = read(fd,r,l);
		if(c!=l){
			free(r);
			r=0;
		}
		else {
		    s=(short*)r;
		    s=&s[l];
		    p=&r[l];
		    for(i=0;i<l;i++){
		    	*p = ((float)*s)/65536.;
		    	p--;
		    	s--;
		    }
		}
		break;
	case CEM_MRC_MODE_FLOAT:
		c = read(fd,r,size);
		if(c!=size){
			free(r);
			r=0;
		}
		break;
	case CEM_MRC_MODE_COMPLEX_INT16:
	case CEM_MRC_MODE_COMPLEX_FLOAT:
		r=0;
		break;
	case CEM_MRC_MODE_USHORT:
		l=size/sizeof(short);
		c = read(fd,r,l);
		if(c!=l){
			free(r);
			r=0;
		}
		else {
		    us=(unsigned short*)r;
		    us=&us[l];
		    p=&r[l];
		    for(i=0;i<l;i++){
		    	*p = ((float)*us)/65536.;
		   		p--;
		   		us--;
		    }
		}
		break;
	default:
		r = 0;
		break;
	}
    close(fd);

	return r;
}



int pmFWriteMRC(const char *filename, PMMRCHeader* header, float *data,size_t size,int mode)
{
	int fd;
	fd = open(filename,O_WRONLY|O_CREAT|O_TRUNC,0644); //open file
	if(fd<0) return 1;
	if(1024!=write(fd,(const void*)header,1024)){  //write header
		close(fd);
		return 1;
	}
	if(size!=write(fd,data,size)){ //write data
		close(fd);
		return 1;
	}
	close(fd);
	return 0;
}


