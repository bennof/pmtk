#ifdef __cplusplus
extern "C" 
#endif

#include "pmtk.h"
#include "protein_io_interface.h"
#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>




//append new io interfaces here
static PMProtein_IO_t root [] = {
	{".pdb",pdbopen,pdbsave},
	{".mrc",mrcopen,mrcsave},
	{".xtc",xtcopen,xtcsave},
	{NULL}
};

static PMProtein_IO_t*  gettype(const char* name)
{
	char *h = name;
	char *p = name;
	char buffer[16];
	PMProtein_IO_t *i;
	size_t l;
	while(*h){
		if(*h=='.') 
			p=h;
		h++;
	}
	if(p==name)  	// no dot in name
		return 0;
	l=strlen(p);
	if(!l || l>15)	// dot at the end or too long type
		return 0;
	
	strcpy(buffer,p);
	p = &buffer[1];
	for ( ; *p; ++p) *p = tolower(*p); // to lower

	for(i=root;i->type;i++){
		if(!strcmp(buffer,i->type))
			return i;
	}
	return 0;
}


PMProtein *pmOpen       (PMProtein *protein,const char *fname)
{
	PMProtein_IO_t *type = gettype(fname);
	if (type)
		return type->open(protein,fname,0);
	else 
		WARN("Unkown file format: %s",fname);
	return 0;
}

PMProtein *pmOpenRef    (PMProtein *protein,const char *fname)
{
	PMProtein_IO_t *type = gettype(fname);
	if (type)
		return type->open(protein,fname,1);
	else 
		WARN("Unkown file format: %s",fname);
	return 0;
}

PMProtein *pmSave       (PMProtein *protein,const char *fname)
{
	PMProtein_IO_t *type = gettype(fname);
	if(type)
		return type->save(protein,fname,0);
	else 
		WARN("Unkown file format: %s",fname);
	return 0;
}
