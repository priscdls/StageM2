#include "structure.h"
#include "utile.h"
#include "output.h"

/**************************************/
/* SHELL ******************************/
/**************************************/

//#define SHL_addVertex(s, id) GPH_addVertex(bond(s), id)

void SHL_initAtom(AtomShl_t* a) {

	flag(a) = -1;

	atomX(a) = 0;
	atomY(a) = 0;
	atomZ(a) = 0;

	parentAtom(a) = -1;

	neighborhood(a) = LST_create();
}

int SHL_nbNeighborhood(AtomShl_t* a) {

	return LST_nbElements(neighborhood(a));
}

int SHL_getIndiceFreeNeighbor(AtomShl_t* a) {

	return LST_getIndiceFree(a->neighborhood);
}

int SHL_getIndice(AtomShl_t* a, unsigned id) {

	return LST_getIndice(a->neighborhood, id);
}

void SHL_addNeighbor(AtomShl_t* a, unsigned id) {

	LST_addElement(a->neighborhood, id);
}

void SHL_removeNeighbor(AtomShl_t* a, unsigned id) {
	
	LST_removeElement(a->neighborhood, id);
}

void SHL_addAllocAtom(Shell_t* s) {
	int i;

	s->atoms = realloc(s->atoms,
		(size(s)+REALLOCSIZE)*sizeof(AtomShl_t));

	for (i=0; i<REALLOCSIZE; i++) {
		SHL_initAtom(atom(s,size(s)+i));
	}

	size(s) += REALLOCSIZE;
}

int SHL_nbAtom(Shell_t* s) {
	int i, cpt = 0;

	for (i=0; i<size(s); i++)
		if (flag(atom(s,i)) != -1)
			cpt++;

	return cpt;
}

int SHL_nbEdges(Shell_t* s) {
	int i, cpt = 0;

	for (i=0; i<size(s); i++)
		cpt += SHL_nbNeighborhood(atom(s, i));

	return cpt/2;
}

int SHL_getIndiceFreeAtom(Shell_t* s) {
	int i;

	for (i=0; i<size(s); i++)
		if (flag(atom(s,i)) == -1)
			return i;

	SHL_addAllocAtom(s);
	return i;
}

unsigned SHL_addVertex(Shell_t* s, unsigned id) {

	return GPH_addVertex(bond(s), id);
}

void SHL_removeVertex(Shell_t* s, unsigned id) {

	GPH_removeVertex(bond(s), id);
}

void SHL_addBond(Shell_t* s, unsigned id1, unsigned id2) {

	GPH_addEdge(bond(s), id1, id2);
}

void SHL_removeBond(Shell_t* s, unsigned id1, unsigned id2) {

	GPH_removeEdge(bond(s), id1, id2);
}

void SHL_addEdge(Shell_t* s, unsigned id1, unsigned id2) {

	if (id1 < size(s) && id2 <size(s) && id1 != id2) {

		SHL_addNeighbor(atom(s, id1), id2);
		SHL_addNeighbor(atom(s, id2), id1);
	}
}

void SHL_removeEdge(Shell_t* s, unsigned id1, unsigned id2) {

	if (id1 < size(s) && id2 <size(s)) {
		
		SHL_removeNeighbor(atom(s, id1), id2);
		SHL_removeNeighbor(atom(s, id2), id1);
	}
}

unsigned SHL_addAtom(Shell_t* s, Point_t coords, unsigned parent) {

	unsigned indice = SHL_getIndiceFreeAtom(s);

	flag(atom(s,indice)) = 0;
	coords(atom(s,indice)) = coords;
	parentAtom(atom(s,indice)) = parent;

	return indice;
}

void SHL_removeAtom(Shell_t* s, unsigned id) {

	int i;

	if (id < size(s)) {
		AtomShl_t* a = atom(s,id);

		if (cycle(s,id))
			LST_removeElement(s->cycle, id);

		for (i=0; i<neighborhoodSize(a); i++)
			if (neighbor(a, i) != -1)
				SHL_removeNeighbor(atom(s, neighbor(a,i)), id);

		LST_delete(neighborhood(a));

		if (checkVertex(s,id))
			SHL_removeVertex(s, id);
		SHL_initAtom(a);
	}
}

//Relie les éléments de la liste l2 au plus proche de la liste l1
void SHL_avoir2(Shell_t* s, List_t* l1, List_t* l2) {

	int i, j, indiceMin;
	float distMin, dis;

	for (i=0; i<size(l2) && elts(l2,i)!=-1; i++) {

		indiceMin = elts(l1,0);
		distMin = dist(coords(atom(s, elts(l2,i))), coords(atom(s, elts(l1,0))));
		for (j=1; j<size(l1) && elts(l1,j)!=-1; j++) {

			dis = dist(coords(atom(s, elts(l2,i))), coords(atom(s, elts(l1,j))));
			if (dis < distMin) {
				indiceMin = elts(l1,j);
				distMin = dis;
			}
		}


		if (flag(atom(s, indiceMin)) != 1)
			SHL_addEdge(s, elts(l2,i), indiceMin);
	}
}

List_t* SHL_seekBorder(Shell_t* s, List_t* in, unsigned id) {

	int i;
	List_t* out = LST_create();
	AtomShl_t* a = atom(s, id);

	if (flag(a) == 0 || flag(a) == 1) {
		LST_addElement(out, id);
		return out;
	}

	LST_addElement(in, id);
	for (i=0; i<neighborhoodSize(a) && neighbor(a,i)!=-1; i++) {
		if (!LST_check(in, neighbor(a,i)))
			out = LST_addList(out, SHL_seekBorder(s, in, neighbor(a,i)));
	}

	return out;
}

//remplacement de avoir2
void SHL_linkBorder(Shell_t* s, unsigned id, List_t* l) {

	int i, j, indiceMin;
	float distMin, dis;
	List_t* tmp = LST_create();
	List_t* border = SHL_seekBorder(s, tmp, id);

	LST_delete(tmp);

	for (i=0; i<size(l) && elts(l,i)!=-1; i++) {

		indiceMin = elts(border,0);
		distMin = dist(coords(atom(s, elts(l,i))), coords(atom(s, elts(border,0))));
		for (j=1; j<size(border) && elts(border,j)!=-1; j++) {

			dis = dist(coords(atom(s, elts(l,i))), coords(atom(s, elts(border,j))));
			if (dis < distMin) {
				indiceMin = elts(border,j);
				distMin = dis;
			}
		}

		SHL_addEdge(s, elts(l,i), indiceMin);
	}

	LST_delete(border);
}


void SHL_addCycle(Shell_t* s, unsigned id) {

	LST_addElement(s->cycle, id);
}

void SHL_mergeAtom(Shell_t* s, unsigned eater, unsigned eaten) {

	int i;
	AtomShl_t* a;

	if (eater < size(s) && eaten < size(s) && eater != eaten) {

		a = atom(s, eaten);
		
		for (i=0; i<neighborhoodSize(a); i++) {
			if (neighbor(a, i) != -1 ) {
				SHL_addEdge(s, eater, neighbor(a,i));
			}
		}
		parentAtom(atom(s, eater)) = parentAtom(a);

		if (cycle(s, eaten))
			SHL_addCycle(s, eater);

		if (flag(atom(s, eater)) < flag(atom(s, eaten)))
			flag(atom(s, eater)) = flag(atom(s, eaten));

		SHL_removeAtom(s, eaten);
	}
}

void SHL_mergeAtom2(Shell_t* s, unsigned id1, unsigned id2) {

	AtomShl_t* a1, *a2;

	if (id1 != id2 && id1 < size(s) && id2 < size(s)) {

		a1 = atom(s, id1);

		a2 = atom(s, id2);

		while (neighbor(a2,0) != -1){
			SHL_addEdge(s, id1, neighbor(a2,0));
			SHL_removeEdge(s, id2, neighbor(a2,0));
			if (flag(a1) < flag(a2))
				flag(a1) = flag(a2);
		}

		coords(a1) = merPoint(coords(a1), coords(a2));

		SHL_removeAtom(s, id2);
	}
}

Shell_t* SHL_create() {

	Shell_t *a = malloc(sizeof(Shell_t));

	a->size = 0;
	a->atoms = NULL;
	a->cycle = LST_create();
	a->bond = GPH_create();

	return a;
}

Shell_t* SHL_copy(Shell_t* s) {

	int i;
	Shell_t* copy = malloc(sizeof(Shell_t));

	size(copy) = size(s);
	copy->atoms = malloc(size(copy)*sizeof(AtomShl_t));
	copy->cycle = LST_copy(s->cycle);
	copy->bond = GPH_copy(s->bond);


	for (i=0; i<size(s); i++) {
		
		flag(atom(copy,i)) = flag(atom(s,i));
		coords(atom(copy,i)) = coords(atom(s,i));
		parentAtom(atom(copy,i)) = parentAtom(atom(s,i));
		neighborhood(atom(copy,i)) = LST_copy(neighborhood(atom(s,i)));
	}

	return copy;
}

//Approxiation des distances
//Création d'un shell de sommets.
Shell_t* SHL_avoir(Shell_t* s) {
	int i, j, k, nb;
	float AB, dis;
	AtomShl_t* A, *B;
	Shell_t* out = SHL_create();

	//size(out) = size(s);
	out->atoms = malloc(size(s)*sizeof(AtomShl_t));

	for (i=0; i<size(s); i++) {
		SHL_addAtom(out, coords(atom(s,i)), -1);
	}

	for (i=0; i<size(s); i++) {
		A = atom(s,i);
		for (j=0; j<neighborhoodSize(A); j++) {
			B = atom(s, neighbor(A,j));
			if (i < neighbor(A,j)) {
				AB = dist(coords(A), coords(B));
				nb = roundf(AB/SIMPLE);
				dis = AB/nb;

				for (k=1; k<nb; k++)
					SHL_addAtom(
						out,
						addPoint(
							coords(A),
							normalization( vector(coords(A),coords(B)), k*dis)),
						-1);
			}
		}
	}

	return out;
}

void SHL_testDis(Shell_t* s) {

int i, j;
	for (i=0; i<size(s); i++) {

		if (flag(atom(s,i)) != -1 && flag(atom(s,i)) != 2)
			for (j=i+1; j<size(s); j++) {
				if (flag(atom(s,j)) != -1 && flag(atom(s,j)) != 2) {
					if (dist(coords(atom(s,i)), coords(atom(s,j))) < MINDIS)
						SHL_mergeAtom2(s, i, j);
				}
			}
	}
}

void SHL_deleteAtom(AtomShl_t* a) {

	LST_delete(neighborhood(a));
}

void SHL_delete(Shell_t* s) {
	
	int i;

	if (s->atoms != NULL) {
		for (i=0; i<size(s); i++)
			SHL_deleteAtom(atom(s,i));
		free(s->atoms);
	}

	if (s->cycle != NULL)
		LST_delete(s->cycle);

	if (s->bond != NULL)
		GPH_delete(s->bond);

	free(s);
}
