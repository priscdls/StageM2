#include "structure.h"

Graph_t* MolToGph(Molecule_t* m) {

	int i, j;
	Graph_t* g = GPH_create();
	GPH_addAlloc(g, size(m));

	for (i=0; i<size(m); i++) {

		GPH_addVertex(g, i);
		for (j=0; j<neighborhoodSize(atom(m,i)); j++) {
			if (neighbor(atom(m,i),j) != -1)
				GPH_addNeighbor(vertex(g,i), neighbor(atom(m,i),j));
		}
	}

	return g;
}

Graph_t* ShlToGph(Shell_t* s) {

	int i, j;
	Graph_t* g = GPH_create();

	for (i=0; i<size(s); i++) {

		if (flag(atom(s,i)) != -1) {

			GPH_addVertex(g, i);
			for (j=0; j<neighborhoodSize(atom(s,i)); j++) {
				if (neighbor(atom(s,i),j) != -1)
					GPH_addNeighbor(vertex(g,i), neighbor(atom(s,i),j));
			}
		}
	}

	return g;
}
