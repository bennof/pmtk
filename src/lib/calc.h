/***

    ./src/lib/calc.h 

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





#ifndef CALC_H_
#define CALC_H_
#ifdef __cplusplus
extern "C" {
#endif


#include "vec.h"


int writeArray(const char *fname, const char *comment, float *array, size_t l);


int addArray2D  (float **r,float **a,float**b,size_t cols,size_t rows);
int subArray2D  (float **r,float **a,float**b,size_t cols,size_t rows);
int multArray2D (float **r,float **a,float**b,size_t cols,size_t rows);
int divArray2D  (float **r,float **a,float**b,size_t cols,size_t rows);

int scaleArray2D(float **r,float **a,float s,size_t cols,size_t rows);
int translateArray2D(float **r,float **a,float s,size_t cols,size_t rows);

int scaleVec3Array2D(float **r,float **a,vec3 s,size_t cols,size_t rows);
int translateVec3Array2D(float **r,float **a,vec3 s,size_t cols,size_t rows);

int sqrtArray2D(float **r,float **a,size_t cols,size_t rows);
int invsqrtArray2D(float **r,float **a,size_t cols,size_t rows);

int addNoiseArray2D(float **r,float **a,float sigma,size_t cols,size_t rows);

vec3 meanVec3Array2D(float **data, size_t cols, size_t rows);
int centerAtVec3Array2D(float **r,float **a,vec3 s,size_t cols,size_t rows);
int rotateVec3Array2D(float **r,float **a,mat3 m,size_t cols,size_t rows);
float alignAtVec3Array2D(float **r,float **a,float *ref,size_t cols,size_t rows);

int rampArray2D(float **r,float **a,uvec3 dim, vec3 voxel,size_t cols,size_t rows);
int maskArray2D(float **r,float **a,float *mask,float thr,size_t cols,size_t rows);
int thresholdArray2D(float **r,float **a,float thr,size_t cols,size_t rows);

int meanArray2D(float **r,float *mean, size_t cols, size_t rows);
int meanVarArray2D(float **r,float *mean,float *variance,size_t cols,size_t rows);

int correlationArray(float *r,float *a,float *b,size_t l);
float rmsdArray(float *X, float *Y, size_t size, mat3 *rot);

int normalizeArray2D(float **r,float **a,size_t cols,size_t rows);

float* covarianceArray2D(float **data,float *mean,size_t cols, size_t rows);
float* inverseCovarianceArray2D(float **data,float *mean,size_t cols, size_t rows);


int covarianceBTArray(float *out, float *cov, float *mean, float **data,size_t id, size_t cols, size_t rows);
int trajectoryArray2D(float **out, float *center,float* vector, float scale, size_t cols,size_t rows);

//LAPACK wrapper
void calcEigenVectors(float *matrix, float *vals, int dim);
void calcEigenVectors_d(double *matrix, double *vals, int dim);

#ifdef USE_FFTW3
int dftr2cArray2D(float **inout,uvec3 dim,size_t *size,size_t frames);
int dftc2rArray2D(float **inout,uvec3 dim,size_t *size,size_t frames);
int dftFilterCutOffArray2D(float **data,uvec3 dim,float lb,float ub,size_t frames);
int dftshiftArray2D(float **data,uvec3 dim, float dx,float dy, float dz,size_t cols);
int dftFilterRampArray2D(float **data,uvec3 dim,vec3 voxel,size_t frames);
int dftFilterLaplaceArray2D	(float **data,uvec3 dim,size_t frames);
int dftFilterGaussArray2D(float **data,uvec3 dim,vec3 apix,float sigma,size_t frames);
#endif

#ifdef __cplusplus
}
#endif
#endif /* CALC_H_ */
