#ifdef __cplusplus
extern "C" 
#endif

#include "pmtk.h"
#include "protein_atomic_info.h"
#include "config.h"




static PMProteinAtomDesc *PM_INFO_STACK=0;

PMProteinAtomDesc* pmCreateAtomDesc(size_t natoms, PMAtomDesc *atoms)
{
	PMProteinAtomDesc* i;
	i = (PMProteinAtomDesc*)malloc(sizeof(PMProteinAtomDesc));
	i->refcount=1;
	i->natoms=natoms;
	i->atoms=atoms;
	i->next=PM_INFO_STACK;
	PM_INFO_STACK=i->next;
	return i;
}

PMProteinAtomDesc* pmDeleteAtomDesc(PMProteinModel *protein)
{
	PMProteinAtomDesc *i,*p;
	p = protein->desc;

	if(!p)
		return 0;

	p->refcount--;
	if(!p->refcount){
		//remove info
		i=PM_INFO_STACK;
		if(i==p){
			PM_INFO_STACK = i->next;
		}
		else{
			while(i && (i->next != p))
				i=i->next;
			if(i)
				i->next=p->next;
		}
		free(p->atoms);
		free(p);
	}
	return 0;
}

PMProteinAtomDesc* pmDublicateAtomDesc(PMProteinModel *protein)
{
	if(protein->desc)
		protein->desc->refcount++;
	return protein->desc;
}
