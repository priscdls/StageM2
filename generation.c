#include "generation.h"
#include "utile.h"
#include "output.h"

void generationDep(Main_t* m) {

	unsigned i, j, copy;
	Shell_t* mo;
	Vertex_t* v;

	//Création des différents moc en fonction des dépendances.
	//mocSize change au fur et à mesure des itérations
	for (i=0; i</*mocSize(m)*/1	; i++) { //Pour tous les mocs
		mo = moc(m,i);
		if (size(moc(m,i)) != 0) {
			for (j=0; j<size(bond(moc(m,i))); j++) {//Pour tous les sommets
				v = vertex(bond(moc(m,i)),j);

				if (id(v) != -1 && nbNeighbors(v) > 0) {

					while (nbNeighbors(v) != 0) {
						copy = MN_copyMoc(m, mo);
						//printf("copy %d \n", copy);
						//printf("number copy %d\n", copy);
						SHL_removeVertex(moc(m, copy), id(v));
						//printf("neighbor(v,0) = %d\n", neighbor(v,0));
						SHL_removeVertex(mo, neighbor(v,0));
					}
				}
			}
		}
		
	}
	printf("mocSize %d\n", mocSize(m));

	//Création des différents moc en fonction des types de D/A.
}

void checkInsertVertex(Shell_t* m, List_t* l, unsigned idv) {

	int i, indice = idv;
	float min, distance;
	AtomShl_t *v = atom(m,idv), *s;

	if (size(l) > 0 && elts(l,0) != -1) {
		min = dist(coords(v), coords(atom(m, elts(l,0))));
		indice = 0;
	}

	for (i=1; i<size(l) && elts(l,i) != -1; i++) {
		distance = dist(coords(v), coords(atom(m, elts(l,i))));
		if (min > distance) {
			min = distance;
			indice = i;
		}
	}

	if (indice != idv && min < MINDIS) {
		s = atom(m, indice);

		while (neighbor(s,0) != -1){
			LST_addElement(l,neighbor(s,0));
			SHL_removeEdge(m, idv, neighbor(s,0));
		}

		flag(v) = flag(s);
		SHL_removeVertex(m, idv);
	}

}

void insertDonor1(Shell_t* m, unsigned idv, Point_t normal, Point_t dir) {
	//La position de v est la position du premier.
	int indice;
	Point_t new_coords;
	AtomShl_t* v = atom(m,idv);
	List_t* l = LST_create();

	flag(v) = 3;
	dir = subPoint(initPoint(0), normalization(dir, SIMPLE));

	//Insérer tous les voisins de v dans la liste.
	while (neighbor(v,0) != -1){
		LST_addElement(l,neighbor(v,0));
		SHL_removeEdge(m, idv, neighbor(v,0));
	}

	//Position du deuxième : (rotation de normal, 120, -dir) + coords(v)
	new_coords = addPoint(coords(v), rotation(normal, 120, dir));
	indice = SHL_addAtom(m, new_coords, -1);

	//checkInsertVertex(m, l, indice);
	SHL_addEdge(m, idv, indice);
	if (flag(atom(m,indice)) < 2)
		flag(atom(m,indice)) = 1;
	
	v = atom(m,idv);
	
	//Position du troisième : (rotation de normal, -120, -dir) + coords(v)
	new_coords = addPoint(coords(v), rotation(normal, -120, dir));
	indice = SHL_addAtom(m, new_coords, -1);

	//checkInsertVertex(m, l, indice);
	SHL_addEdge(m, idv, indice);
	if (flag(atom(m,indice)) < 2)
		flag(atom(m,indice)) = 1;

	//Rattacher les nouveaux sommets à ceux de la liste l.
	SHL_linkBorder(m, idv, l);

	LST_delete(l);
}

void insertAcceptor1(Shell_t* m, unsigned idv, Point_t normal, Point_t dir) {
	//La position de v est la position du premier.
	int indice, idc;
	Point_t new_coords;
	AtomShl_t* v = atom(m, idv), *c;
	List_t* l = LST_create();

	flag(v) = 3;

	//Insérer tous les voisins de v dans la liste.
	while (neighbor(v,0) != -1){
		LST_addElement(l,neighbor(v,0));
		SHL_removeEdge(m, idv, neighbor(v,0));
	}


	//Position du deuxième sommet : centre du motif
	//Hydrogène+taille d'une liaison simple moyenne.
	dir = normalization(dir, (SIMPLE/2)+(MINDIS/2));
	new_coords = addPoint(coords(v), dir);
	idc = SHL_addAtom(m, new_coords, -1);

	//checkInsertVertex(m, l, idc);
	SHL_addEdge(m, idv, idc);
	c = atom(m,idc);
	flag(c) = 3;

	//Position du deuxième : (rotation de normal, 120, -dir) + coords(v)
	dir = subPoint(initPoint(0), normalization(dir, SIMPLE));
	new_coords = addPoint(coords(c), rotation(normal, 120, dir));
	indice = SHL_addAtom(m, new_coords, -1);

	//checkInsertVertex(m, l, indice);
	SHL_addEdge(m, idc, indice);
	if (flag(atom(m,indice)) < 2)
		flag(atom(m,indice)) = 1;

	//Position du troisième : (rotation de normal, 120, -dir) + coords(v)
	new_coords = addPoint(coords(c), rotation(normal, -120, dir));
	indice = SHL_addAtom(m, new_coords, -1);

	//checkInsertVertex(m, l, indice);
	SHL_addEdge(m, idc, indice);
	if (flag(atom(m,indice)) < 2)
		flag(atom(m,indice)) = 1;

	//Rattacher les nouvaux sommets à ceux de la liste l.
	SHL_linkBorder(m, idc, l);
	
	LST_delete(l);
}

void insertAcceptor2(Shell_t* m, unsigned idv, Point_t normal, Point_t dir) {
	//La position de v est la position du premier.
	int indice, idc;
	Point_t x1, x2, x3, x4;
	AtomShl_t* v = atom(m, idv), *c;
	List_t* l = LST_create();

	flag(v) = 3;

	//Insérer tous les voisins de v dans la liste.
	while (neighbor(v,0) != -1){
		LST_addElement(l,neighbor(v,0));
		SHL_removeEdge(m, idv, neighbor(v,0));
	}

	x1 = coords(v);

	//Position du deuxième sommet : centre du motif
	//Hydrogène+taille d'une liaison simple moyenne.
	dir = normalization(dir, (SIMPLE/2)+(MINDIS/2));
	x2 = addPoint(coords(v), dir);
	idc = SHL_addAtom(m, x2, -1);

	//checkInsertVertex(m, l, idc);
	SHL_addEdge(m, idv, idc);
	c = atom(m,idc);
	flag(c) = 3;

	//Deuxième sommet du tétraèdre
	//x2 = AX1E3(coords(c), x1, normal, SIMPLE);
	x2 = normalization(vector(x1, coords(c)), 1);
	x2 = rotation(normal, 109.47, x2);
	x2 = normalization(x2, SIMPLE);
	x2 = addPoint(coords(c), x2);
	indice = SHL_addAtom(m, x2, -1);

	//checkInsertVertex(m, l, indice);
	SHL_addEdge(m, idc, indice);
	flag(atom(m,indice)) = 1;

	//Troisième sommet du tétraèdre
	x3 = AX2E2(coords(c), x1, x2, SIMPLE);
	indice = SHL_addAtom(m, x3, -1);

	//checkInsertVertex(m, l, indice);
	SHL_addEdge(m, idc, indice);
	flag(atom(m,indice)) = 1;

	//Troisième sommet du tétraèdre
	x4 = AX3E1(coords(c), x1, x2, x3, SIMPLE);
	indice = SHL_addAtom(m, x4, -1);

	//checkInsertVertex(m, l, indice);
	SHL_addEdge(m, idc, indice);
	flag(atom(m,indice)) = 1;

	//Rattacher les nouvaux sommets à ceux de la liste l.
	SHL_linkBorder(m, idc, l);
	
	LST_delete(l);
}

void generationHydro(Main_t* m) {

	int i, j, idv;
	AtomShl_t *v;
	Atom_t *parent;

	for (i=0; i<1; i++) {
		for (j=0; j<size(bond(moc(m,i))); j++) {
			idv = id(vertex(bond(moc(m,i)),j));

			if (idv != -1) {
				//printf("idv = %d\n", idv);
				v = atom(moc(m,i), idv);
				parent = atom(substrat(m), parentAtom(v));
				if (!strcmp(symbol(parent), "H")) {
					insertDonor1(moc(m,i), idv, MOL_seekNormal(substrat(m), parentAtom(v), -1), vector(coords(parent), coords(v)));
				}
				else {
					if (steric(parent) == 3)
						insertAcceptor1(moc(m,i), idv, MOL_seekNormal(substrat(m), parentAtom(v), -1), vector(coords(parent), coords(v)));
					else
						insertAcceptor2(moc(m,i), idv, MOL_seekNormal(substrat(m), parentAtom(v), -1), vector(coords(parent), coords(v)));
				}
			}
		}

		SHL_testDis(moc(m,i));
	}
}


void generationCycle(Shell_t* s) {
	int i, j;
	AtomShl_t* a;
	List_t* nei;
	List_t* atomT = LST_create();

	for (i=0; i<size(s); i++) {
		if (flag(atom(s,i)) != -1)
			if (cycle(s, i))				
				LST_addElement(atomT, i);
	}

	for (i=0; i<size(atomT) && elts(atomT,i) != -1; i++) {

		a = atom(s, elts(atomT,i));
		nei = LST_create();

		//pour tous les voisins de a
		for (j=0; j<neighborhoodSize(a) && neighbor(a,j) != -1; j++) {

			if (!LST_check(atomT, neighbor(a,j)) ||
			dist(coords(a), coords(atom(s,neighbor(a,j)))) > 1.7) {

				LST_addElement(nei,neighbor(a,j));

			}
		}
		
		//S'il existe au moins deux voisins de a pouvant participés au motif.
		if (SHL_nbNeighborhood(a) - LST_nbElements(nei) > 1) {

			//Retirer les anciens liens entre l'atome et ses voisins
			for (j=0; j<size(nei) && elts(nei,j) != -1; j++)
				SHL_removeEdge(s, elts(atomT,i), elts(nei,j));

			flag(a) = 2;
			if (SHL_nbNeighborhood(a) == 2) {
				int id = -1;
				Point_t newCoords = autre(coords(a), coords(atom(s,neighbor(a,0))),
							coords(atom(s,neighbor(a,1))), 1.4);
				id = SHL_addAtom(s, newCoords, -1);

				for (j=0; j<size(s); j++)
					if (flag(atom(s,j)) != -1 &&
					dist(newCoords, coords(atom(s,j))) < 0.7){
						if (LST_check(nei, j)) LST_removeElement(nei, j);
						if (cycle(s,j)) LST_addElement(atomT, id);
						SHL_mergeAtom(s, id, j);
						//S'il le sommet appartenant à la liste de cycle il faut le rajouter dans atomT
					}
				if (flag(atom(s,id)) < 2) flag(atom(s,id)) = 1;

				//si autre sommet trop proche newCoord == coords(sommet)

				//Ajoute les arêtes entre le nouveau sommet et i	
				SHL_addEdge(s, elts(atomT,i), id);
			}

			SHL_linkBorder(s, elts(atomT,i), nei);
		}
		LST_delete(nei);
	}
	LST_delete(atomT);
}

void generationMoc(Main_t* m) {

	//unsigned idMoc = MN_copyMoc(m, envelope(m));
	//printf("idMoc %d\n", idMoc);

	//SHL_write(moc(m, idMoc));
	envarom(m) = SHL_copy(envelope(m));
	//generationDep(m);
	printf("Génération des motifs aromatiques.\n");
	generationCycle(envarom(m));
	//SHL_write(envarom(m));

	printf("Génération des motifs hydrogènes.\n");
	MN_copyMoc(m, envarom(m));
	generationDep(m);
	generationHydro(m);

}
