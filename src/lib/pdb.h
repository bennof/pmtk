/***

    ./src/lib/pdb.h 

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



#ifndef PDB_H_
#define PDB_H_
#ifdef __cplusplus
extern "C" {
#endif


#include "protein_atomic_info.h"

#define PMATOM_TYPE_ATOM       'A'
#define PMATOM_TYPE_HETATM     'H'


/**
 * structure stores pdb like atom information
 */
//typedef struct _PMAtomDesc{
//	char type; 		/**< Record name id ATOM or HETATOM @see coord types */
//	char serial[6];		/**< Atom  serial number */
//	char name[4];		/**< Atom name */
//	char altLoc;		/**< Alternate location indicator */
//	char resName[4];	/**< Residue name */
//	char chainID;		/**< Chain identifier */
//	char resSeq[4];		/**< Residue sequence number */
//	char iCode;		/**< Code for insertion of residues */
//} PMAtomDesc;


int getType(const char* name);

#define PM_ATOM_FILE_READ  0x02
#define PM_ATOM_FILE_WRITE 0x04
#define PM_ATOM_FILE_XDR   0x20
#define PM_ATOM_FILE_PDB   0x40


int    pmFWriteFramesPDB   (FILE* stream,PMAtomDesc* info,float **data,size_t models,size_t atoms);
float *pmFReadPDBFramePlain(FILE* stream,size_t *atoms);
float *pmFReadPDBFrameFull (FILE* stream,PMAtomDesc** atominfo,size_t *atoms);


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

#ifdef __cplusplus
}
#endif
#endif /* PDB_H_ */
