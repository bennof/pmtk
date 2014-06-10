#ifndef _ROTEIN_ATOMIC_INFO_H
#define _ROTEIN_ATOMIC_INFO_H

#ifdef __cplusplus
extern "C" {
#endif

#include "pmtk.h"

#define PMATOM_TYPE_ATOM       'A'
#define PMATOM_TYPE_HETATM     'H'


/**
 * structure stores pdb like atom information
 */
typedef struct _PMAtomDesc{
	char type; 		/**< Record name id ATOM or HETATOM @see coord types */
	char serial[6];		/**< Atom  serial number */
	char name[4];		/**< Atom name */
	char altLoc;		/**< Alternate location indicator */
	char resName[4];	/**< Residue name */
	char chainID;		/**< Chain identifier */
	char resSeq[4];		/**< Residue sequence number */
	char iCode;		/**< Code for insertion of residues */
} PMAtomDesc;

//typedef struct _PMProteinAtomDesc PMProteinAtomDesc;
struct _PMProteinAtomDesc{
	size_t refcount;
	size_t natoms;
	PMProteinAtomDesc *next;
	PMAtomDesc *atoms;
};

PMProteinAtomDesc* pmCreateAtomDesc(size_t natoms, PMAtomDesc *atoms);
PMProteinAtomDesc* pmDeleteAtomDesc(PMProteinModel *protein);
PMProteinAtomDesc* pmDublicateAtomDesc(PMProteinModel *protein);


#ifdef __cplusplus
}
#endif
#endif
