#include "structure.h"
#include "output.h"

void MN_initMoc(Shell_t* s) {

	s->size = 0;
	s->atoms = NULL;
	s->cycle = NULL;
	s->bond = NULL;
}

void MN_addAlloc(Main_t* m, unsigned size) {

	int i;

	m->mocs = realloc(m->mocs, (mocSize(m)+size)*sizeof(Shell_t*));

	for (i=0; i<size; i++) {
		moc(m,mocSize(m)+i) = NULL;
	}

	mocSize(m) += size;
	//printf("size MN %d\n", mocSize(m));
}

unsigned MN_getIndiceFree(Main_t* m) {

	int i;

	for (i=0; i<mocSize(m); i++) {
		if (moc(m,i) == NULL)
			return i;
	}
	//printf("%p\n", m->mocs);
	MN_addAlloc(m, REALLOCSIZE);
	//printf("%p\n", m->mocs);
	return i;
}

unsigned MN_getIndiceFree2(Main_t* m) {

	int i;

	for (i=0; i<mocSize(m); i++) {
		if (moc(m,i) == NULL)
			return i;
	}
	//printf("%p\n", m->mocs);
	MN_addAlloc(m, 1);
	//printf("%p\n", m->mocs);
	return i;
}

unsigned MN_copyMoc(Main_t* m, Shell_t* s) {

	unsigned indice = MN_getIndiceFree(m);
	//printf("test copyMoc indice %d\n", indice);
	moc(m, indice) = SHL_copy(s);

	return indice;
}

Main_t* MN_create() {

	Main_t* m = malloc(sizeof(Main_t));
	m->substrat = NULL;
	m->envelope = NULL;
	m->envarom = NULL;
	m->mocs = NULL;
	m->mocSize = 0;

	return m;
}

void MN_delete(Main_t* m) {

	int i;

	if (substrat(m) != NULL)
		MOL_delete(substrat(m));

	if (envelope(m) != NULL)
		SHL_delete(envelope(m));

	if (envarom(m) != NULL)
		SHL_delete(envarom(m));

	if (m->mocs != NULL) {
		for (i=0; i<mocSize(m); i++) {
			if (moc(m,i) != NULL)
				SHL_delete(moc(m,i));
		}
		free(m->mocs);
	}

	free(m);
}
