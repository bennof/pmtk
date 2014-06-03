/***

    ./src/lib/vec.h 

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





#ifndef VEC_H_
#define VEC_H_
#ifdef __cplusplus
extern "C" {
#endif

#include <sys/types.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.1415926f

static inline float invsqrt(float x){
	long i;
	float y  = x*0.5f;
	i  = *(long*)&x;
	i  = 0x5f3759df-(i>>1);
	x  = *(float*)&i;
	return x*(1.5f-(y*x*x));
}

/** Fast inverse square root (Quake) */
static inline double dinvsqrt(double x){
  long long i;
  double y;
  y = x*0.5f;
  i = *(long long*)&x;
  i = 0x5fe6ec85e7de30daLL-(i>>1);
  x = *(double*)&i;
  return x*(1.5f-(y*x*x));
}

#define randf() (((float)rand())/RAND_MAX)
static inline float boxmullerrand() {
	float x1, x2, w;
	w=1.1;
	while(w>=1.0){
		x1 = 2.0*randf()-1.0;
		x2 = 2.0 *randf()-1.0;
		w = x1 * x1 + x2 * x2;
	}
	w = sqrt( (-2.0 * log( w ) ) / w );
	return x1 *w;
}


/** Fast gauss  */
static inline float gauss(float x){
  int j = 0x40000000 - 0x0004f4*((int) (x*(x*5000+50)));
  return (*(float*)&j)*0.5f;
}



typedef struct _vec2{
	float x;
	float y;
} vec2;

typedef struct _vec3{
	float x;
	float y;
	float z;
} vec3;

typedef struct _vec4{
	float x;
	float y;
	float z;
	float w;
} vec4;

typedef struct _ivec2{
	ssize_t x;
	ssize_t y;
} ivec2;

typedef struct _ivec3{
	ssize_t x;
	ssize_t y;
	ssize_t z;
} ivec3;

typedef struct _ivec4{
	ssize_t x;
	ssize_t y;
	ssize_t z;
	ssize_t w;
} ivec4;

typedef struct _uvec2{
	size_t x;
	size_t y;
} uvec2;

typedef struct _uvec3{
	size_t x;
	size_t y;
	size_t z;
} uvec3;

typedef struct _uvec4{
	size_t x;
	size_t y;
	size_t z;
	size_t w;
} uvec4;


typedef struct _quat{
	float r;
	float i;
	float j;
	float k;
} quat;

typedef struct{ //invert the matrix (column major order - Fortran like)
	float xx; float xy;
	float yx;float yy;
}mat2;

typedef struct _mat3 { //invert the matrix (column major order - Fortran like)
	float xx; float xy; float xz;
	float yx; float yy; float yz;
	float zx; float zy; float zz;
} mat3;

typedef struct{	//invert the matrix (column major order - Fortran like)
	float xx; float xy;float xz;float xw;
	float yx; float yy;float yz;float yw;
	float zx; float zy;float zz;float zw;
	float wx; float wy;float wz;float ww;
}mat4;


static inline vec2 VEC2(float x, float y)
{
	register vec2 v;
	v.x=x;
	v.y=y;
	return v;
}

static inline vec3 VEC3(float x, float y,float z)
{
	register vec3 v;
	v.x=x;
	v.y=y;
	v.z=z;
	return v;
}

static inline vec4 VEC4(float x, float y,float z,float w)
{
	register vec4 v;
	v.x=x;
	v.y=y;
	v.z=z;
	v.w=w;
	return v;
}

static inline vec3 vec3add(vec3 a, vec3 b)
{
	register vec3 r;
	r.x=a.x+b.x;
	r.y=a.y+b.y;
	r.z=a.z+b.z;
	return r;
}

static inline vec3 vec3sub(vec3 a, vec3 b)
{
	register vec3 r;
	r.x=a.x-b.x;
	r.y=a.y-b.y;
	r.z=a.z-b.z;
	return r;
}

static inline vec3 vec3mult(vec3 a, vec3 b)
{
	register vec3 r;
	r.x=a.x*b.x;
	r.y=a.y*b.y;
	r.z=a.z*b.z;
	return r;
}

static inline vec3 vec3div(vec3 a, vec3 b)
{
	register vec3 r;
	r.x=a.x/b.x;
	r.y=a.y/b.y;
	r.z=a.z/b.z;
	return r;
}

static inline vec3 vec3adds(vec3 a, float b)
{
	register vec3 r;
	r.x=a.x+b;
	r.y=a.y+b;
	r.z=a.z+b;
	return r;
}

static inline vec3 vec3subs(vec3 a, float b)
{
	register vec3 r;
	r.x=a.x-b;
	r.y=a.y-b;
	r.z=a.z-b;
	return r;
}

static inline vec3 vec3mults(vec3 a, float b)
{
	register vec3 r;
	r.x=a.x*b;
	r.y=a.y*b;
	r.z=a.z*b;
	return r;
}

static inline vec3 vec3divs(vec3 a, float b)
{
	register vec3 r;
	r.x=a.x/b;
	r.y=a.y/b;
	r.z=a.z/b;
	return r;
}




static inline vec2 vec2prodmat2r(mat2 m, vec2 v){
	register vec2 h;
	h.x=m.xx*v.x+m.yx*v.y;
	h.y=m.xy*v.x+m.yy*v.y;
  	return h;
}

static inline vec3 vec3prodmat3r(mat3 m, vec3 v){
	register vec3 h;
	h.x=m.xx*v.x+m.yx*v.y+m.zx*v.z;
	h.y=m.xy*v.x+m.yy*v.y+m.zy*v.z;
	h.z=m.xz*v.x+m.yz*v.y+m.zz*v.z;
	return h;
}

static inline vec4 vec4prodmat4r(mat4 m, vec4 v){
	register vec4 h;
	h.x=m.xx*v.x+m.yx*v.y+m.zx*v.z+m.wx*v.w;
	h.y=m.xy*v.x+m.yy*v.y+m.zy*v.z+m.wy*v.w;
	h.z=m.xz*v.x+m.yz*v.y+m.zz*v.z+m.wz*v.w;
	h.w=m.xw*v.x+m.yw*v.y+m.zw*v.z+m.ww*v.w;
	return h;
}

static inline vec2 vec2prodmat2l(vec2 v, mat2 m){
	register vec2 h;
	h.x=m.xx*v.x+m.xy*v.y;
	h.y=m.yx*v.x+m.yy*v.y;
	return h;
}

static inline vec3 vec3prodmat3l(vec3 v,mat3 m){
	register vec3 h;
	h.x=m.xx*v.x+m.xy*v.y+m.xz*v.z;
	h.y=m.yx*v.x+m.yy*v.y+m.yz*v.z;
  	h.z=m.zx*v.x+m.zy*v.y+m.zz*v.z;
  	return h;
}

static inline vec4 vec4prodmat4l(vec4 v, mat4 m){
	register vec4 h;
	h.x=m.xx*v.x+m.xy*v.y+m.xz*v.z+m.xw*v.w;
	h.y=m.yx*v.x+m.yy*v.y+m.yz*v.z+m.yw*v.w;
	h.z=m.zx*v.x+m.zy*v.y+m.zz*v.z+m.zw*v.w;
	h.w=m.wx*v.x+m.wy*v.y+m.wz*v.z+m.ww*v.w;
	return h;
}

#ifdef __cplusplus
}
#endif
#endif /* VEC_H_ */
