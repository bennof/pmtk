#include "pmtk.h"


#ifndef _PROTEIN_IO_INTERFACE_H
#define _PROTEIN_IO_INTERFACE_H
#ifdef __cplusplus
extern "C" {
#endif

typedef struct _PMProtein_IO_struct{
	const char *type;
	PMProtein * ( * open ) ( PMProtein *protein, const char * filename, int mode );
	PMProtein * ( * save ) ( PMProtein *protein, const char * filename, int mode );
}PMProtein_IO_t;

PMProtein * mrcopen ( PMProteinDensity* protein, const char * filename, int mode );
PMProtein * mrcsave ( PMProteinDensity* protein, const char * filename, int mode );
PMProtein * pdbopen ( PMProteinDensity* protein, const char * filename, int mode );
PMProtein * pdbsave ( PMProteinDensity* protein, const char * filename, int mode );
PMProtein * xtcopen ( PMProteinDensity* protein, const char * filename, int mode );
PMProtein * xtcsave ( PMProteinDensity* protein, const char * filename, int mode );

#ifdef __cplusplus
}
#endif
#endif
