/***

    ./src/lib/calc.c 

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




#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "config.h"
#include "vec.h"

#ifdef __cplusplus
extern "C"
#endif

float rmsdArray(float *X, float *Y, size_t size, mat3 *rot);


int addArray2D(float **r,float **a,float**b,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vb,*vr;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		vb=b[j];
		for(i=rows;i>0;i--)
			(*vr++)=(*va++)+(*vb++);
	}
	return 0;
}

int subArray2D(float **r,float **a,float**b,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vb,*vr;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		vb=b[j];
		for(i=rows;i>0;i--)
			(*vr++)=(*va++)-(*vb++);
	}
	return 0;
}

int multArray2D(float **r,float **a,float**b,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vb,*vr;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		vb=b[j];
		for(i=rows;i>0;i--)
			(*vr++)=(*va++)*(*vb++);
	}
	return 0;
}

int divArray2D(float **r,float **a,float**b,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vb,*vr;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		vb=b[j];
		for(i=rows;i>0;i--){
			if((*vb)!=0)
				(*vr++)=(*va++)/(*vb++);
			else {
				(*vr++)=0.;
				va++;
				vb++;
			}
		}
	}
	return 0;
}



int scaleArray2D(float **r,float **a,float s,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vr;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		for(i=rows;i>0;i--)
			(*vr++)=(*va++)*s;
	}
	return 0;
}

int translateArray2D(float **r,float **a,float s,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vr;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		for(i=rows;i>0;i--)
			(*vr++)=(*va++)+s;
	}
	return 0;
}

int scaleVec3Array2D(float **r,float **a,vec3 s,size_t cols,size_t rows)
{
	size_t i,j,ro;
	vec3 *va,*vr;
	ro=rows/3;
	for(j=0;j<cols;j++){
		vr=(vec3*)r[j];
		va=(vec3*)a[j];
		for(i=ro;i>0;i--)
			(*vr++) = vec3mult((*va++),s);
	}
	return 0;
}

int translateVec3Array2D(float **r,float **a,vec3 s,size_t cols,size_t rows)
{
	size_t i,j,ro;
	vec3 *va,*vr;
	ro=rows/3;
	for(j=0;j<cols;j++){
		vr=(vec3*)r[j];
		va=(vec3*)a[j];
		for(i=ro;i>0;i--)
			(*vr++) = vec3add((*va++),s);
	}
	return 0;
}

int centerAtVec3Array2D(float **r,float **a,vec3 s,size_t cols,size_t rows)
{
	register size_t i,j,ro;
	register vec3 m,*p,*h;

	ro=rows/3;
	m=VEC3(0.,0.,0.);

	for(i=0;i<cols;i++){
		p=(vec3*)a[i];
		for(j=0;j<ro;j++)
			m=vec3add(m,(*p++));
		m=vec3divs(m,ro);
		m=vec3mults(m,-1.);
		m=vec3add(m,s);
		p=(vec3*)a[i];
		h=(vec3*)r[i];
		for(j=0;j<ro;j++)
			(*h++)=vec3add((*p++),m);

	}
	m=vec3divs(m,ro*cols);
	return 0;
}

float alignAtVec3Array2D(float **r,float **a,float *ref,size_t cols,size_t rows)
{
	register size_t i,j,ro;
	register vec3 *va,*vr;
	mat3 m;
	float h=0.;

	ro=rows/3;

	for(i=0;i<cols;i++){
		h += rmsdArray(a[i],ref,rows,&m);
		va=(vec3*)a[i];
		vr=(vec3*)r[i];
		for(j=0;j<ro;j++)
			(*vr++) = vec3prodmat3r(m,(*va++));
	}
	return h/cols;
}



int sqrtArray2D(float **r,float **a,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vr;
	register float h;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		for(i=rows;i>0;i--){
			h=(*va++);
			(*vr++)=(h>0.)?h*invsqrt(h):0.0;
		}
	}
	return 0;
}

int invsqrtArray2D(float **r,float **a,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vr;
	register float h;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		for(i=rows;i>0;i--){
			h=(*va++);
			(*vr++)=(h>0.)?invsqrt(h):0.0;
		}
	}
	return 0;
}

int addNoiseArray2D(float **r,float **a,float sigma,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vr;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		for(i=rows;i>0;i--){
			(*vr++)+=(*va++)+boxmullerrand()*sigma;
		}
	}
	return 0;
}

vec3 meanVec3Array2D(float **data, size_t cols, size_t rows)
{
	register size_t i,j,ro;
	register vec3 m,*p;

	ro=rows/3;
	m=VEC3(0.,0.,0.);

	for(i=0;i<cols;i++){
		p=(vec3*)data[i];
		for(j=0;j<ro;j++)
			m=vec3add(m,(*p++));
	}
	m=vec3divs(m,ro*cols);
	return m;
}

int rotateVec3Array2D(float **r,float **a,mat3 m,size_t cols,size_t rows)
{
	size_t i,j,ro;
	ro=rows/3;
	vec3 *va,*vr;

	for(j=0;j<cols;j++){
		vr=(vec3*)r[j];
		va=(vec3*)a[j];
		for(i=ro;i>0;i--)
			(*vr++) = vec3prodmat3r(m,(*va++));
	}
	return 0;
}



int normalizeArray2D(float **r,float **a,size_t cols,size_t rows)
{
	size_t i,j;
	float *vr,*va;
	register float m,s,h;
	for(j=0;j<cols;j++){
		va=a[j];
		m=0.;
		s=0.;
		for(i=rows;i>0;i--){
			h=(*va++);
			m+=h;
			s+=h*h;
		}
		m/=rows;
		s/=rows;
		s-=m*m;
		s=invsqrt(s);
		vr=r[j];
		va=a[j];
		for(i=rows;i>0;i--){
			h=(*va);
			(*vr++)=(h-m)*s;
		}
	}
	return 0;
}
int rampArray2D(float **r,float **a,uvec3 dim, vec3 voxel,size_t cols,size_t rows){
	size_t i,j,k,f;
	float *va,*vr;
	float h;
	vec3  s,x;
	ivec3 c;

	s.x=voxel.x*voxel.x;
	s.y=voxel.y*voxel.y;
	s.z=voxel.z*voxel.z;

	c.x=dim.x/2;
	c.y=dim.y/2;
	c.z=dim.z/2;

	//iterate through frames
	for(f=0;f<cols;f++){
		va=a[f];
		vr=r[f];
		for(k=0;k<dim.z;k++){
			x.z=((float)((k-c.z)*(k-c.z)))*s.z;
			for(j=0;j<dim.y;j++){
				x.y=((float)((k-c.y)*(k-c.y))) *s.y;
				for(i=0;i<dim.x;i++){
					h= ((float)((i-c.x)*(i-c.x)))  *s.x+x.y+x.z;
					(*vr++)=(*va++)*invsqrt(h);
				}
			}
		}
	}
	return 0;
}
int maskArray2D(float **r,float **a,float *mask,float thr,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vr,*vm;
	register float h;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		vm = mask;
		for(i=rows;i>0;i--){
			h=(*va++);
			(*vr++)=((*vm++)>thr)?h:0.;
		}
	}
	return 0;
}
int thresholdArray2D(float **r,float **a,float thr,size_t cols,size_t rows)
{
	size_t i,j;
	float *va,*vr;
	register float h;
	for(j=0;j<cols;j++){
		vr=r[j];
		va=a[j];
		for(i=rows;i>0;i--){
			h=(*va++);
			(*vr++)=(h>thr)?h:0.;
		}
	}
	return 0;
}



int correlationArray(float *r,float *a,float *b,size_t l)
{
	size_t i;
	register float ma,mb,sxa,sxb,sxy;
	register float x,y;
	float *A,*B;
	ma=mb=sxa=sxb=sxy=0.;

	A=a;
	B=b;
	for(i=l;i>0;i--){
		x=(*A++);
		y=(*B++);
		ma+=x;
		mb+=y;
		sxa+=x*x;
		sxb+=y*y;
		sxy+=x*y;
	}
	ma/=l;
	mb/=l;
	sxa/=l;
	sxb/=l;
	sxy/=l;
	sxa-=ma*ma;
	sxb-=mb*mb;
	sxy-=ma*mb;
	*r=sxy*invsqrt(sxa*sxb);
	return 0;
}

int meanArray2D(float **data,float *mean, size_t cols, size_t rows)
{
	register size_t i,j;
	register float *m,*p;


	memset(mean,0,rows*sizeof(float));
	for(i=0;i<cols;i++){
		p=data[i];
		m=mean;
		for(j=0;j<rows;j++)
			(*m++)+=(*p++);
	}
	for(j=0;j<rows;j++)
		(*m++)/=cols;
	return 0;
}

int meanVarArray2D(float **r,float *mean,float *variance,size_t cols,size_t rows)
{
	size_t i,j;
	register float h,hh,*m,*v,*p;

	memset(mean,0,rows*sizeof(float));
	memset(variance,0,rows*sizeof(float));

	for(j=0;j<cols;j++){
		p=r[j];
		m=mean;
		v=variance;
		for(i=rows;i>0;i--){
			h=(*p++);
			(*m)+=h;
			(*v)+=h*h;
			m++;
			v++;
		}
	}

	m=mean;
	v=variance;
	for(i=rows;i>0;i--){
		h=(*m);
		hh=(*v);
		h/=cols;
		hh/=cols;
		hh-=h*h;
		(*m++)=h;
		(*v++)=hh;
	}
	return 0;
}

float* covarianceArray2D(float **data,float *mean,size_t cols, size_t rows)
{
	register size_t i,j,k;
	register float h1,h2,s, m1, m2;
	float *cov = malloc(rows*rows*sizeof(float));
	if(cov==0) return cov;
	memset(cov,0,rows*rows*sizeof(float));
	for(i=0;i<rows;i++){
		m1=mean[i];
		for(j=i;j<rows;j++){
			m2=mean[j];
			s=0.;
			for(k=0;k<cols;k++){
				h1=data[k][i];
				h2=data[k][j];
				s += (h1-m1) * (h2-m2);
			}
			cov[i*rows+j]=s;
			cov[i+j*rows]=s;
		}
	}
	return cov;
}

float *inverseCovarianceArray2D(float **data, float *mean, size_t cols, size_t rows)
{
	register size_t i,j,k;
	register float h1,h2,s, m;
	float *cov = (float*)malloc(cols*cols*sizeof(float));
	if(cov==0) return cov;
	memset(cov,0,cols*cols*sizeof(float));
	for(i=0;i<cols;i++){
		for(j=i;j<cols;j++){
			s=0.;
			for(k=0;k<rows;k++){
				h1=data[i][k];
				h2=data[j][k];
				m=mean[k];
				s += (h1-m) * (h2-m);
			}
			cov[i*cols+j]=s;
			cov[i+j*cols]=s;
		}
	}
	return cov;
}








#ifdef USE_FFTW3
#include <fftw3.h>
static inline ssize_t dft_c(ssize_t x,ssize_t NX){return (int)x-((int)NX*((int)x/((int)NX/2+1)));}


fftwf_plan dftfree(fftwf_plan plan)
{
  if(plan) fftwf_destroy_plan(plan);
  return 0;
}


int dftr2c(float **inout,fftwf_plan plan,uvec3 dim,size_t *size)
{
	size_t s;
	float *h=0;
	if(plan==0) {
		plan = fftwf_plan_dft_r2c_3d(dim.z,dim.y,dim.x,*inout,(fftwf_complex*)h,FFTW_ESTIMATE);
		s = (dim.x+2)*dim.y*dim.z;
	}
	h = fftwf_malloc(s*sizeof(float));
	if(h==0) return -1;
	fftwf_execute_dft_r2c(plan,*inout,(fftwf_complex*)h);
	*size=s;
	*inout=h;
	return 0;
}

int dftc2r(float **inout,fftwf_plan plan,uvec3 dim,size_t *size)
{
	float *h=0,*hh;
	hh=*inout;
	if(plan==0) plan = fftwf_plan_dft_c2r_3d(dim.z,dim.y,dim.x,(fftwf_complex*)hh,h,FFTW_ESTIMATE);
	size_t s = dim.x*dim.y*dim.z;
	h = fftwf_malloc(s*sizeof(float));
	if(h==0) return -1;
	fftwf_execute_dft_c2r(plan,(fftwf_complex*)hh,h);
	*size=s;
	*inout=h;
	return 0;
}


int dftr2cArray2D(float **inout,uvec3 dim,size_t *size,size_t frames)
{
	fftwf_plan plan=0;
	size_t i;
	for(i=0;i<frames;i++){
		if(-1==dftr2c(inout,plan,dim,size))
			return -1;
		inout++;
	}
	return 0;
}

int dftc2rArray2D(float **inout,uvec3 dim,size_t *size,size_t frames)
{
	fftwf_plan plan;
	size_t i;
	for(i=0;i<frames;i++){
		if(-1==dftc2r(inout,plan,dim,size))
			return -1;
		inout++;
	}
	return 0;
}

int dftshiftArray2D(float **data,uvec3 dim, float dx,float dy, float dz,size_t cols)
{
	float *d;
	size_t i,j,k,xx,yy,zz,p,f;
	double re[6],im[6],h;

	re[0] = cos(2.0f*PI/dim.x*dx);
	im[0] = -1.0f*sin(2.0f*PI/dim.x*dx);

	re[1] = cos(2.0f*PI/dim.y*dy);
	im[1] = -1.0f*sin(2.0f*PI/dim.y*dy);

	re[2] = cos(-1.0f*PI*dy);
	im[2] = -1.0f*sin(-1.0f*PI*dy);

	re[5] = cos(2.0f*PI/dim.z*dz);
	im[5] = -1.0f*sin(2.0f*PI/dim.z*dz);

	re[4] = cos(-1.0f*PI*dz);
	im[4] = -1.0f*sin(-1.0f*PI*dz);

	xx=dim.x/2+1;

	//iterate through frames
	for(f=0;f<cols;f++){
		d = data[f];
		re[4] = cos(-1.0f*PI*dz);
		im[4] = -1.0f*sin(-1.0f*PI*dz);
		for(k=0;k<dim.z;k++){
		    zz=(k+dim.z/2)%dim.z;
		    re[2]=re[4];
		    im[2]=im[4];

		    for(j=0;j<dim.y;j++){
		    	yy=(j+dim.y/2)%dim.y;
		        re[3]=re[2];
		        im[3]=im[2];
		        for(i=0;i<xx;i++){
		        	p=2*(i+yy*xx+zz*dim.x*xx);

		        	h=d[p];
		        	d[p  ] = h*re[3]-d[p+1]*im[3];
		        	d[p+1] = d[p+1]*re[3]+h*im[3];

		        	h=re[3];
		        	re[3] = h*re[0]-im[3]*im[0];
		        	im[3] = h*im[0]+im[3]*re[0];
		        }
		        h=re[2];
		    	re[2] = h*re[1]-im[2]*im[1];
		    	im[2] = h*im[1]+im[2]*re[1];
		    }
		    h=re[4];
			re[4] = h*re[5]-im[4]*im[5];
			im[4] = h*im[5]+im[4]*re[5];
		}
	}
	return 0;
}



int dftFilterCutOffArray2D(float **data,uvec3 dim,float lb,float ub,size_t frames)
{
	size_t i,j,k,f;
	ssize_t u,v;
	const float k0 = 1./(dim.x*dim.x);
	const float k1 = 1./(dim.y*dim.y);
	const float k2 = 1./(dim.z*dim.z);
	float h, *m;
	//lb/=protein->apix[0];
	//ub/=protein->apix[0];

	lb=1./(lb*lb);
	ub=1./(ub*ub);

	for(f=0;f<frames;f++){
		m = data[j];
		for(k=0;k<dim.z;k++){
			v=dft_c(k,dim.z);
			for(j=0;j<dim.y;j++){
				u=dft_c(j,dim.y);
				for(i=0;i<(dim.x/2+1);i++){
					h=(float)(i*i)*k0+(float)(u*u)*k1+(float)(v*v)*k2;
					if(h>=ub && h<=lb){
						(m++) ;
						(m++) ;
					}
					else{
						(*m++) = 0.0;
						(*m++) = 0.0;
					}
				}
			}
	    }
	}
	return 0;
}

int dftFilterGaussArray2D(float **data,uvec3 dim,vec3 apix,float sigma,size_t frames)
{
	size_t i,j,k,f;
	ssize_t u,v;
	float h, *m;
	const float k0 = 1.0/(apix.x*apix.x*dim.x*dim.x);
	const float k1 = 1.0/(apix.y*apix.y*dim.y*dim.y);
	const float k2 = 1.0/(apix.z*apix.z*dim.z*dim.z);
	const float C = sqrt(2*PI)*sigma;
	const float c = -2.0*PI*PI*sigma*sigma;

	//iterate through frames
	for(f=0;f<frames;f++){
		m=data[f];
		for(k=0;k<dim.z;k++){
			v=dft_c(k,dim.z);
			for(j=0;j<dim.y;j++){
				u=dft_c(j,dim.y);
				for(i=0;i<(dim.x+1);i++){
					h=(float)(i*i)*k0+(float)(u*u)*k1+(float)(v*v)*k2;
					(*m++) *= C*exp(c*h);
					(*m++) *= C*exp(c*h);
				}
			}
		}
	}
	return 0;
}

int dftFilterLaplaceArray2D	(float **data,uvec3 dim,size_t frames)
{
	size_t i,j,k,f;
	ssize_t u,v;
	float  *m;
	const float k0 = -4.0*PI*PI/(dim.x*dim.x);
	const float k1 = -4.0*PI*PI/(dim.y*dim.y);
	const float k2 = -4.0*PI*PI/(dim.z*dim.z);
	float h;

	for(f=0;f<frames;f++){
		m=data[f];
	  	for(k=0;k<dim.z;k++){
	  		v=dft_c(k,dim.z);
	  		for(j=0;j<dim.y;j++){
	  			u=dft_c(j,dim.y);
	  			for(i=0;i<(dim.x/2+1);i++){
	  				h=((float)(i*i)*k0+(float)(u*u)*k1+(float)(v*v)*k2);
	  				(*m++) *= h;
	  				(*m++) *= h;
	  			}
	  		}
	  	}
	}
	return 0;
}

int dftFilterRampArray2D(float **data,uvec3 dim,vec3 voxel,size_t frames)
{
	size_t i,j,k,f;
	ssize_t u,v;
	float h, *m;
	const float k0 = 1.0/(voxel.x*voxel.x*dim.x*dim.x);
	const float k1 = 1.0/(voxel.y*voxel.y*dim.y*dim.y);
	const float k2 = 1.0/(voxel.z*voxel.z*dim.z*dim.z);

	for(f=0;f<frames;f++){
		m=data[f];
		for(k=0;k<dim.z;k++){
			v=dft_c(k,dim.z);
			for(j=0;j<dim.y;j++){
				u=dft_c(j,dim.y);
				for(i=0;i<(dim.x/2+1);i++){
					h=(float)(i*i)*k0+(float)(u*u)*k1+(float)(v*v)*k2;
					h=(h==0.0)?0.0:h*invsqrt(h);
					(*m++) *= h;
					(*m++) *= h;
				}
			}
	    }
    }
	return 0;
}

#endif


int covarianceBTArray(float *out, float *cov, float *mean, float **data,size_t id, size_t cols, size_t rows)
{
	size_t i,j;
	register float h,s,m,*c;

	c=&cov[id*cols];
	h=0.;
	for(i=0;i<rows;i++){
		m=mean[i];
		s=0.;
		for(j=0;j<cols;j++){
			s+=(data[j][i]-m)*c[j];
		}
		out[i]=s;
		h+=s*s;
	}
	h=invsqrt(h);
	for(i=0;i<rows;i++)
		out[i]*=h;
	return 0;
}

int trajectoryArray2D(float **out, float *center,float* vector, float scale, size_t cols,size_t rows)
{
	float *h,*c,*v,s;
	size_t i,j;

	scale/=cols;
	s=(cols/2-cols)*scale;

	for(i=0;i<cols;i++){
		h=out[i];
		c=center;
		v=vector;
		s += scale;
		for(j=0;j<rows;j++){
			(*h++)=(*c++)+s*(*v++);
		}
	}
	return 0;
}


extern void ssyev_(char* JOBZp,char*  UPLOp,int* Np, float* Ap, int* LDAp, float* Wp, float*  WORKp, int* LWORKp,int *IWORKp, int * LIWORKp, int* INFOp);
extern void dsyev_(char* JOBZp,char*  UPLOp,int* Np, double* Ap, int* LDAp, double* Wp, double*  WORKp, int* LWORKp,int *IWORKp, int * LIWORKp, int* INFOp);

void calcEigenVectors(float *matrix, float *vals, int dim)
{
  char u='U';
  char v='V';
  int lwork=-1, liwork=-1;
  int *iwork;
  float *work;

  iwork = (int*) malloc((dim)*sizeof(int));
  work = (float*) malloc((dim)*sizeof(float));
  int info;

  //fprintf(stderr,"dim=%i\n",dim);
  ssyev_(&v, &u, &dim, matrix, &dim, vals, work, &lwork, iwork, &liwork,&info);

  lwork  = (1>(int) work[0])?1:(int) work[0];
  liwork = (1>(int)iwork[0])?1:(int)iwork[0];
  if(iwork) free((void*)iwork);
  if( work) free((void*)work);
  iwork = (int*) malloc((liwork)*sizeof(int));
  work = (float*) malloc((lwork)*sizeof(float));

  ssyev_(&v, &u, &dim, matrix, &dim, vals, work, &lwork, iwork, &liwork,&info);
  if(iwork)free((void*)iwork);
  if( work)free((void*)work);
}

void calcEigenVectors_d(double *matrix, double *vals, int dim)
{
  char u='U';
  char v='V';
  int lwork=-1, liwork=-1;
  int *iwork;
  double *work;

  (iwork = (int*) malloc((dim)*sizeof(int)));
  (work = (double*) malloc((dim)*sizeof(double)));
  int info;

  dsyev_(&v, &u, &dim, matrix, &dim, vals, work, &lwork, iwork, &liwork,&info);

  lwork  = (1>(int) work[0])?1:(int) work[0];
  liwork = (1>(int)iwork[0])?1:(int)iwork[0];
  if(iwork) free((void*)iwork);
  if( work) free((void*)work);
  (iwork = (int*) malloc((liwork)*sizeof(int)));
  (work = (double*) malloc((lwork)*sizeof(double)));

  dsyev_(&v, &u, &dim, matrix, &dim, vals, work, &lwork, iwork, &liwork,&info);
  if(iwork)free((void*)iwork);
  if( work)free((void*)work);
}



float rmsdArray(float *X, float *Y, size_t size, mat3 *rot){
	  size_t i, n;
	  double GA,GB;
	  double A[16];
	  double xx,xy,xz,yx,yy,yz,zx,zy,zz;
	  double val[4];
	  mat3 B;

	  n=size/3;

	  GA=GB=xx=xy=xz=yx=yy=yz=zx=zy=zz=0.0f;
	  for(i=0;i<n;i++){
	    xx += (double) X[i*3  ]*Y[i*3  ];
	    xy += (double) X[i*3  ]*Y[i*3+1];
	    xz += (double) X[i*3  ]*Y[i*3+2];
	    yx += (double) X[i*3+1]*Y[i*3  ];
	    yy += (double) X[i*3+1]*Y[i*3+1];
	    yz += (double) X[i*3+1]*Y[i*3+2];
	    zx += (double) X[i*3+2]*Y[i*3  ];
	    zy += (double) X[i*3+2]*Y[i*3+1];
	    zz += (double) X[i*3+2]*Y[i*3+2];
	    GA += (double) X[i*3  ]*X[i*3  ]+X[i*3+1]*X[i*3+1]+X[i*3+2]*X[i*3+2];
	    GB += (double) Y[i*3  ]*Y[i*3  ]+Y[i*3+1]*Y[i*3+1]+Y[i*3+2]*Y[i*3+2];
	  }

	  A[0]  = xx+yy+zz;
	  A[5]  = xx-yy-zz;
	  A[10] = -xx+yy-zz;
	  A[15] = -xx-yy+zz;

	  A[1]  = yz-zy;
	  A[2]  = zx-xz;
	  A[3]  = xy-yx;
	  A[6]  = xy+yx;
	  A[7]  = xz+zx;
	  A[11] = yz+zy;
	  A[4]=A[1];
	  A[8]=A[2];
	  A[12]=A[3];
	  A[9]=A[6];
	  A[13]=A[7];
	  A[14]=A[11];

	  calcEigenVectors_d(A,val,4);

	  if(rot){
	  	static int w=12;
	  	static int x=13;
	  	static int y=14;
	  	static int z=15;
	  	B.xx = (float)(A[w]*A[w]+A[x]*A[x]-A[y]*A[y]-A[z]*A[z]);
	  	B.xy = (float)( 2*(A[w]*A[z]+A[x]*A[y] ));
	  	B.xz = (float)( 2*(A[x]*A[z]-A[w]*A[y] ));

	  	B.yx = (float)( 2*(A[x]*A[y]-A[w]*A[z] ));
	  	B.yy = (float)(A[w]*A[w]-A[x]*A[x]+A[y]*A[y]-A[z]*A[z]);
	  	B.yz = (float)( 2*(A[w]*A[x]+A[y]*A[z] ));

	  	B.zx = (float)( 2*(A[w]*A[y]+A[x]*A[z] ));
	  	B.zy = (float)( 2*(A[y]*A[z]-A[w]*A[x] ));
	  	B.zz = (float)(A[w]*A[w]-A[x]*A[x]-A[y]*A[y]+A[z]*A[z]);

	  	if(rot) *rot=B;
	  }

	  xx=(((GA+GB-2.0f*val[3])/n)<=0.0f)?0.0f:sqrt(((GA+GB-2.0f*val[3])/n));
	  return (float) xx;
}









