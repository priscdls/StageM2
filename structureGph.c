#include "structure.h"

void GPH_initVertex(Vertex_t* v) {

	v->id = -1;

	neighborhood(v) = LST_create();
	nbNeighbors(v) = 0;
}

void GPH_addNeighbor(Vertex_t* v, unsigned id) {

	LST_addElement(neighborhood(v), id);
	nbNeighbors(v)++;
}

void GPH_removeNeighbor(Vertex_t* v, unsigned id) {

	LST_removeElement(neighborhood(v), id);
	nbNeighbors(v)--;
}

void GPH_deleteVertex(Vertex_t* v) {

	LST_delete(neighborhood(v));
}

int GPH_nbVertex(Graph_t* g) {

	int i, cpt = 0;
	for (i=0; i<size(g); i++)
		if (id(vertex(g,i)) != -1)
			cpt++;

	return cpt;
}

void GPH_addAlloc(Graph_t* g, unsigned size) {

	int i;

	g->vertices = realloc(g->vertices, (size(g)+size)*sizeof(Vertex_t));

	for (i=0; i<size; i++)
		GPH_initVertex(vertex(g,size(g)+i));

	size(g) += size;
}

int GPH_getIndiceFree(Graph_t* g) {

	int i;

	for (i=0; i<size(g); i++)
		if (g->vertices[i].id == -1)
			return i;

	GPH_addAlloc(g, REALLOCSIZE);
	return i;
}

int GPH_getIndice(Graph_t* g, unsigned id) {

	int i;

	for (i=0; i<size(g); i++)
		if (g->vertices[i].id == id)
			return i;

	return -1;
}

unsigned GPH_addVertex(Graph_t* g, unsigned id) {

	unsigned indice = GPH_getIndice(g, id);

	if (indice == -1) {

		indice = GPH_getIndiceFree(g);
		id(vertex(g,indice)) = id;
	}

	return indice;
}

void GPH_removeVertex(Graph_t* g, unsigned id) {

	int i, indice = GPH_getIndice(g, id);

	if (indice != -1) {

		Vertex_t* v = vertex(g, indice);

		for (i=0; i<nbNeighbors(v); i++)
			GPH_removeNeighbor(vertex(g, GPH_getIndice(g, neighbor(v,i))), id);
	
		GPH_deleteVertex(v);
		GPH_initVertex(v);
	}
}

void GPH_addEdge(Graph_t* g, unsigned id1, unsigned id2) {

	int indice1 = GPH_getIndice(g, id1), indice2 = GPH_getIndice(g, id2);
	if (indice1 != -1 && indice2 != -1) {

		GPH_addNeighbor(vertex(g, indice1), id2);
		GPH_addNeighbor(vertex(g, indice2), id1);
	}
}

void GPH_removeEdge(Graph_t* g, unsigned id1, unsigned id2) {

	int indice1 = GPH_getIndice(g, id1), indice2 = GPH_getIndice(g, id2);
	if (indice1 != -1 && indice2 != -1) {

		GPH_removeNeighbor(vertex(g, indice1), id2);
		GPH_removeNeighbor(vertex(g, indice2), id1);
	}
}

int GPH_cycle(Graph_t* g, List_t* l, unsigned id, unsigned idP) {

	int i, tmp = 0;
	Vertex_t* v;

	if (id == elts(l,0))
		return 1;
	if (LST_nbElements(l) > 5)
		return 0;

	LST_addElement(l, id);
	v = vertex(g, GPH_getIndice(g, id));

	for (i=0; i<nbNeighbors(v) && tmp == 0; i++) {
		if (neighbor(v, i) != idP)
			tmp = GPH_cycle(g, l, neighbor(v, i), id);
	}

	LST_removeElement(l, id);
	return tmp;
}

List_t* GPH_seekCycle(Graph_t* g) {

	int i;
	List_t* out = LST_create();
	List_t* l = LST_create();
	Vertex_t* v, *n;

	for (i=0; i<size(g); i++) {
		if (nbNeighbors(vertex(g, i)) == 1)
			LST_addElement(l, id(vertex(g,i)));
	}

	while (LST_nbElements(l) != 0) {

		v = vertex(g, GPH_getIndice(g, elts(l,0)));
		n = vertex(g, GPH_getIndice(g, neighbor(v,0)));

		if (nbNeighbors(n) == 2)
			LST_addElement(l, id(n));

		GPH_removeVertex(g, elts(l,0));
		LST_removeElement(l, elts(l,0));
	}

	//chercher les cycles
	for (i=0; i<size(g); i++) {

		if (id(vertex(g,i))!=-1 && GPH_cycle(g, l, id(vertex(g,i)), -1))
			LST_addElement(out, id(vertex(g,i)));
	}

	LST_delete(l);

	return out;
}

unsigned GPH_checkVertex(Graph_t* g, unsigned id) {

	int i;

	for (i=0; i<size(g); i++)
		if (id(vertex(g,i)) == id)
			return 1;
	
	return 0;
}

unsigned GPH_checkBond(Graph_t* g, unsigned id1, unsigned id2) {

	return LST_check(neighborhood(vertex(g, GPH_getIndice(g, id1))), id2);
}

Graph_t* GPH_create() {

	Graph_t* g = malloc(sizeof(Graph_t));

	g->size = 0;
	g->vertices = NULL;

	return g;
}

void GPH_delete(Graph_t* g) {

	int i;

	for (i=0; i<size(g); i++)
		GPH_deleteVertex(vertex(g,i));

	free(g->vertices);
	free(g);
}

Graph_t* GPH_copy(Graph_t* g) {

	int i;
	Graph_t* copy = GPH_create();
	unsigned indice;

	GPH_addAlloc(copy, GPH_nbVertex(g));

	for (i=0; i<size(g); i++) {
		if (id(vertex(g,i)) != -1) {
			indice = GPH_addVertex(copy, id(vertex(g,i)));
			LST_delete(neighborhood(vertex(copy, indice)));
			neighborhood(vertex(copy, indice)) = LST_copy(neighborhood(vertex(g,i)));
			nbNeighbors(vertex(copy, indice)) = nbNeighbors(vertex(g,i));
		}
	}

	return copy;
}
