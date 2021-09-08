#include "expansion.h"
#include "utile.h"
#include "interface.h"

#include "output.h"

/**
* Expension de groupement atomique sterique 2.
* Ajout d'un point dans l'enveloppe.
* 
* @param m Molécule (entrée)
* @param s Enveloppe (sortie)
* @param id Identifiant du sommet traité de la molécule.
*/
void expansion_AX1E1(Molecule_t* m, Shell_t* s, unsigned id) {

	unsigned indice;
	Point_t newCoords;
	Atom_t* a = atom(m,id), *x1 = atom(m,neighbor(a,0));

  newCoords = AX1E1(coords(a), coords(x1), HYDRO);

  indice = SHL_addAtom(s, newCoords, id);

  if (checkVertex(m,id))
			SHL_addVertex(s, indice);
}

/**
* Expension de groupement atomique sterique 3.
* Ajout de un à quatre points dans l'enveloppe.
* 
* @param m Molécule (entrée)
* @param s Enveloppe (sortie)
* @param id Identifiant du sommet traité de la molécule.
*/
void expansion_steric3(Molecule_t* m, Shell_t* s, unsigned id) {

	unsigned indice;
	Point_t a, x1, x2, x3, normal;
	a = coords(atom(m,id));
	x1 = coords(atom(m,neighbor(atom(m,id),0)));

	//Ajout d'un point si le groupement possède moins de deux doublets.
	//X2 devient soit le point ajouté à l'enveloppe, soit le second voisin dans la molécule.
	if (ligands(atom(m,id)) < 2) {

		if (neighbor(atom(m,neighbor(atom(m,id),0)),0) == id)
			x2 = coords(atom(m,neighbor(atom(m,neighbor(atom(m,id),0)),1)));
		else
			x2 = coords(atom(m,neighbor(atom(m,neighbor(atom(m,id),0)),0)));

		normal = planNormal(a, x1, x2);
		x2 = AX1E2(a, x1, normal, HYDRO);
		indice = SHL_addAtom(s, x2, id);

		if (checkVertex(m,id))
			SHL_addVertex(s, indice);

	}
	else
		x2 = coords(atom(m,neighbor(atom(m,id),1)));

	//Ajout d'un point si le groupement possède moins de trois doublets.
	//X3 devient soit le point ajouté à l'enveloppe, soit le troisième voisin dans la molécule.
	if (ligands(atom(m,id)) < 3) {
		x3 = AX2E1(a, x1, x2, HYDRO);
		indice = SHL_addAtom(s, x3, id);

		if (checkVertex(m,id))
			SHL_addVertex(s, indice);
	}
	else
		x3 = coords(atom(m,neighbor(atom(m,id),2)));

	//Ajout des deux extensions perpendiculaires.
	normal = normalization(planNormal(x1, x2, x3), HYDRO);

	if (cycle(m,id)) {
		SHL_addCycle(s, SHL_addAtom(s, addPoint(a, normal), id));
		SHL_addCycle(s, SHL_addAtom(s, subPoint(a, normal), id));
	}
	else {
		SHL_addAtom(s, addPoint(a, normal), id);
		SHL_addAtom(s, subPoint(a, normal), id);
	}
}

/**
* Expension de groupement atomique sterique 4.
* Ajout de zéro à trois points dans l'enveloppe.
* 
* @param m Molécule (entrée)
* @param s Enveloppe (sortie)
* @param id Identifiant du sommet traité de la molécule.
*/
void expansion_steric4(Molecule_t* m, Shell_t* s, unsigned id) {

	unsigned indice;
	Point_t a, x1, x2, x3, normal;

	a = coords(atom(m,id));
	x1  = coords(atom(m, neighbor(atom(m,id), 0)));

	//Ajout d'un point si le groupement possède moins de deux doublets.
	//X2 devient soit le point ajouté à l'enveloppe, soit le second voisin dans la molécule.
	if (ligands(atom(m,id)) < 2) {

		if (neighbor(atom(m,neighbor(atom(m,id),0)),0) == id)
			x2 = coords(atom(m,neighbor(atom(m,neighbor(atom(m,id),0)),1)));
		else
			x2 = coords(atom(m,neighbor(atom(m,neighbor(atom(m,id),0)),0)));

		//à voir si on le garde
		x2 = addPoint(a, planNormal(a, x1, x2));
		normal = planNormal(a, x1, x2);
		x2 = AX1E3(a, x1, normal, HYDRO);
		indice = SHL_addAtom(s, x2, id);

		if (checkVertex(m,id))
			SHL_addVertex(s, indice);
	}
	else
		x2 = coords(atom(m, neighbor(atom(m,id), 1)));

	//Ajout d'un point si le groupement possède moins de trois doublets.
	//X3 devient soit le point ajouté à l'enveloppe, soit le troisième voisin dans la molécule.
	if (ligands(atom(m,id)) < 3) {

		x3 = AX2E2(a, x1, x2, HYDRO);
		indice = SHL_addAtom(s, x3, id);

		if (checkVertex(m,id))
			SHL_addVertex(s, indice);
	}
	else
		x3 = coords(atom(m, neighbor(atom(m,id), 2)));

	if (ligands(atom(m,id)) < 4) {

		indice = SHL_addAtom(s, AX3E1(a, x1, x2, x3, HYDRO), id);

		if (checkVertex(m,id))
			SHL_addVertex(s, indice);
	}
}

/********************* A refaire ************************/
/**
* Expension de groupement atomique avec un angle à 180°.
* Ajout de quatre points dans l'enveloppe.
* 
* @param m Molécule (entrée)
* @param s Enveloppe (sortie)
* @param id Identifiant du sommet traité de la molécule.
*/
void expansion_AX2E0(Molecule_t* m, Shell_t* s, unsigned id) {
	Point_t normal, newCoords;
	Atom_t* a = atom(m,id),	*x1 = atom(m,neighbor(a,0)),
	*x2 = atom(m, neighbor(a,1));

	normal = normalization(vector(coords(a), coords(x1)),1);

	newCoords = normalization(
		planNormal(coords(a), coords(x1), coords(x2)), HYDRO);

	SHL_addAtom(s, addPoint(coords(a), newCoords), id);

	SHL_addAtom(s, addPoint(coords(a), normalization(
		rotation(normal, 90, newCoords), HYDRO)), id);

	SHL_addAtom(s, addPoint(coords(a), normalization(
		rotation(normal, 180, newCoords), HYDRO)), id);

	SHL_addAtom(s, addPoint(coords(a), normalization(
		rotation(normal, -90, newCoords), HYDRO)), id);
}

/**
* Expansion de l'ensemble de la molécule.
* Ajout de zéro à trois points dans l'enveloppe.
* 
* @param m Molécule (entrée) Déjà instanciée.
* @param s Enveloppe (sortie) créée.
*/
void expansion(Molecule_t* m, Shell_t* s) {
	int i, j;

	//Création du nuage de points constituant l'enveloppe.
	/*rajouter une liste pour être sur qu'ils soient fait dans le bon ordre*/
	/*vérifier que toutes les donnés ont été calculé*/
	for (i=0; i<size(m); i++) {
		if (ligands(atom(m,i)) == 1 && lonePairs(atom(m,i)) == 1)
			expansion_AX1E1(m, s, i);
		else if (steric(atom(m,i)) == 3)
			expansion_steric3(m, s, i);
		else if (steric(atom(m,i)) == 4)
			expansion_steric4(m, s, i);
		else if (ligands(atom(m,i)) == 2 && lonePairs(atom(m,i)) == 0)
			expansion_AX2E0(m, s, i);
	}

	//Construction d'arrêtes du graphe.
	for (i=0; i<size(bond(s)) && id(vertex(bond(s),i)) != -1; i++)
		for (j=i+1; j<size(bond(s)) && id(vertex(bond(s),j)) != -1; j++) {
			if (parentAtom(atom(s,id(vertex(bond(s),i)))) == parentAtom(atom(s,id(vertex(bond(s),j))))
					|| checkBond(m,
						parentAtom(atom(s,id(vertex(bond(s),i)))), 
						parentAtom(atom(s,id(vertex(bond(s),j)))))
					)
				SHL_addBond(s, id(vertex(bond(s),i)), id(vertex(bond(s),j)));
		}

}

/**
* Construction des arêtes de l'enveloppe 
* Appel à R.
* 
* @param s			Enveloppe possédant déjà un nuage de points.
* @param alpha 	Paramètre de la sphère. (2 est souvent une bonne mesure, 3 sinon)
*/
void alphaShape(Shell_t* s, double alpha) {
	int i;
	Ashape_t* as3d = Cashape3d(s, alpha);

	for (i=0; i<(as3d->nb_edge/2); i++) {
		SHL_addEdge(s, as3d->edge[i]-1, as3d->edge[i+as3d->nb_edge/2]-1);
	}

	for (i=0; i<size(s); i++) {
		if (neighborhoodSize(atom(s,i)) == 0) {
			SHL_removeAtom(s, i);
			GPH_removeVertex(bond(s),i);
		}
	}
	
	ASP_delete(as3d);
}

Ashape_t* alphaShape2(Shell_t* s, double alpha) {
	int i;
	Ashape_t* as3d = Cashape3d(s, alpha);

	for (i=0; i<(as3d->nb_edge/2); i++) {
		SHL_addEdge(s, as3d->edge[i]-1, as3d->edge[i+as3d->nb_edge/2]-1);
	}

	for (i=0; i<size(s); i++) {
		if (neighborhoodSize(atom(s,i)) == 0) {
			SHL_removeAtom(s, i);
			GPH_removeVertex(bond(s),i);
		}
	}
	return as3d;
}

/**
* Création complète de l'enveloppe à partir de la molécule.
* 
* @param m Molécule (entrée) Déjà instanciée.
* @param s Enveloppe (sortie) créée.
*/
Shell_t* createShell(Molecule_t* m, double alpha) {
	Shell_t* s = SHL_create();
	
	expansion(m, s);
	
	SHL_writeMol2("Results/vec.mol2", s);
	alphaShape(s, alpha);
	printf("Graphe de dépendance de l'enveloppe.\n");
	GPH_write(bond(s));	

	return s;
}

Shell_t* createShell2(Molecule_t* m, double alpha, Ashape_t** as3d) {
	Shell_t* s = SHL_create();
	
	expansion(m, s);
	
	SHL_writeMol2("Results/vec.mol2", s);
	*as3d = alphaShape2(s, alpha);
	printf("Graphe de dépendance de l'enveloppe.\n");
	GPH_write(bond(s));	

	return s;
}
