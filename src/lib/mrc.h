/***

    ./src/lib/mrc.h 

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





#ifndef MRC_H_
#define MRC_H_
#ifdef __cplusplus
extern "C" {
#endif





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


void setMRCHeader(PMMRCHeader *header,size_t dim[], float apix[], float origin[]);

/**
 * read a mrc or ccp4 file
 * @param filename filename of the file
 * @param header pointer to mrc header info
 * @param size size of the return data to check in bytes, (size=0 will be ignored)
 * @param mode not used
 * @return data as float array
 */
float *pmFReadMRC(const char *filename, PMMRCHeader* header,size_t size,int mode);

/**
 * writes to a mrc or ccp4 file
 * @param filename filename of the file
 * @param header mrc header info
 * @param data data as float array
 * @param size of data in bytes
 * @param mode not used
 * @return
 */
int pmFWriteMRC(const char *filename, PMMRCHeader* header, float *data,size_t size,int mode);





typedef struct _MRC_FILE {
	int fd;
	int mode;
	size_t size;
} MRC_FILE;


MRC_FILE *pmOpenMrc(const char *filename, PMMRCHeader* header, int mode);
int pmReadMRC(MRC_FILE *fp,float *data,size_t size);
int pmWriteMRC(MRC_FILE *fp,float *data, size_t size);
int pmCloseMRC(MRC_FILE *fp);


#ifdef __cplusplus
}
#endif
#endif /* MRC_H_ */
