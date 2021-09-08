#include "structure.h"
#include "utile.h"

/**************************************/
/* MOLECULE ***************************/
/**************************************/

/*
*	Fonction qui ajoute l'identifiant d'un voisin à un atome
*
*	@param 	a 	Pointeur de l'atome
* @param	id 	Identifiant à rajouter
*
*/
void MOL_addNeighbor(Atom_t* a, unsigned id) {

	if (LST_getIndice(neighborhood(a), id) == -1) {
		int indice = LST_getIndiceFree(neighborhood(a));

		if (indice == -1) {
			printf("MOL_addNeighbor : -1 \n");
			exit(1);
		}

		neighbor(a, indice) = id;
	}
}

/*
*	Fonction qui supprime un voisin d'un atome.
*
*	@param 	a 	Pointeur de l'atome
* @param	id 	Identifiant à retirer
*
*/
void MOL_removeNeighbor(Atom_t* a, unsigned id) {
	int indice = LST_getIndice(neighborhood(a), id);

	if (indice != -1)
		neighbor(a, indice) = -1;
}

/*
*	Fonction qui calcule le nombre de doublets liants d'un atome.
* Nombre de doublets liants = Nombre de voisins dans la molécule.
*
*	@param 	a 	Pointeur de l'atome.
*
*/
void MOL_nbLigands(Atom_t* a) {
	int i, cpt = 0;

	for (i=0; i<4; i++)
		if (neighbor(a,i) != -1)
			cpt++;

	ligands(a) = cpt;
}

/*
*	Fonction qui calcule le nombre de doublets non liants d'un atome.
*
*	@param 	a 			Pointeur de l'atome.
* @param	alpha 	Moyenne dans angles formés par les voisins de l'atome.
* @param	stericNeighbor	Nombre de doublets du voisin de l'atome (pour le cas où il n'en possède qu'un seul.
*	@param	cycle 	Booléen qui indique si l'atome appartient à un cycle.
*/
void MOL_nbLonePairs(Atom_t* a, float alpha, int stericNeighbor, unsigned cycle) {

	if (ligands(a) == 1)
	{
		if (!strcmp(symbol(a), "H")) {
			lonePairs(a) = 1;
		}
		else if(!strcmp(symbol(a), "Cl") || !strcmp(symbol(a), "Br")
			|| !strcmp(symbol(a), "F") || !strcmp(symbol(a), "I")) {
			lonePairs(a) = 3;
		}
		else {
				lonePairs(a) = stericNeighbor - 1;
		}
	}
	else {
		if (ligands(a) == 4)
			lonePairs(a) = 0;
		
		else if (abs(120-alpha) < 4)
			lonePairs(a) = 3 - ligands(a);

		else if (alpha-109 < 7)
			if (cycle == 1 && stericNeighbor == 3)
				lonePairs(a) = 3 - ligands(a);
			else if (cycle == 1 && stericNeighbor == -1)
				lonePairs(a) = -1;
			else
				lonePairs(a) = 4 - ligands(a);
		else 
			lonePairs(a) = 0;
	}
}

/*
* Fonction qui trouve tous les sommets appartenant à un cycle d'une molécule.
*
* @param	m 	Adresse de la molécule.
*/
void MOL_seekCycle(Molecule_t* m) {

	Graph_t* g = MolToGph(m);

	m->cycle = GPH_seekCycle(g);

	GPH_delete(g);
}

/*
*	Fonction qui calcule ne nombre d'arêtes (de liaisons) d'un molécule.
*
*	@param 	m 	Pointeur de la molécule.
*
*/
int MOL_nbEdges(Molecule_t* m) {
	int i, cpt = 0;

	for (i=0; i<size(m); i++)
		cpt += ligands(atom(m,i));

	return cpt/2;
}

/*
*	Fonction qui ajoute une arête entre deux sommets d'une molécule.
*
*	@param 	m 	Pointeur de la molécule.
* @param	id1	Identifiant du premier sommet.
* @param	id2	Identifiant du deuxième sommet.
*
*/
void MOL_addEdge(Molecule_t* m, unsigned id1, unsigned id2) {
	if (id1 != id2) {
		MOL_addNeighbor(atom(m,id1), id2);
		MOL_addNeighbor(atom(m,id2), id1);
	}
}

/*
*	Fonction qui supprime une arête entre deux sommets d'une molécule.
*
*	@param 	m 	Pointeur de la molécule.
* @param	id1	Identifiant du premier sommet.
* @param	id2	Identifiant du deuxième sommet.
*
*/
void MOL_removeEdge(Molecule_t* m, unsigned id1, unsigned id2) {
	
	MOL_removeNeighbor(atom(m,id1), id2);
	MOL_removeNeighbor(atom(m,id2), id1);
}

/*
*	Fonction qui initialise un nouvel atome.
*
*	@param 	a 	Pointeur de l'atome.
*/
void MOL_createAtom(Atom_t* a) {

	symbol(a)[0] = '\0';
	radius(a) = -1;
	ligands(a) = -1;
	lonePairs(a) = -1;

	atomX(a) = 0;
	atomY(a) = 0;
	atomZ(a) = 0;

	neighborhood(a) = LST_create();
}

/*
* Fonction qui crée le graphe de dépendance d'une molécule.
* Le graphe de dépendance regroupe les atomes de la molécule qui peut participer à une liaison hydrogène.
* Il y a une arête entre deux sommets s'ils ne peuvent pas participer à des liaisons hydrogènes en même temps.
*
*	@param	m 	Molécule.
*
*/
void MOL_createBond(Molecule_t* m) {

	int i, j, k;
	Atom_t* a;
	List_t* lh;

	for (i=0; i<size(m); i++) {
		a = atom(m, i);

		if (!strcmp(symbol(a), "O") || !strcmp(symbol(a), "N")
			||  !strcmp(symbol(a), "F")
			/*&& (ligands(a) != 4 && lonePairs(a) != 0)*/) {
			lh = LST_create();
			if (lonePairs(a) > 0 /*sauf (2,2) && (1,3)*/) {
				LST_addElement(lh, i);
				GPH_addVertex(bond(m), i);
			}

			for (j=0; j<neighborhoodSize(a) && neighbor(a,j) != -1; j++)
				if (!strcmp(symbol(atom(m, neighbor(a,j))), "H")) {
					LST_addElement(lh, neighbor(a,j));
					GPH_addVertex(bond(m), neighbor(a,j));
				}

			for (j=0; j<size(lh) && elts(lh,j)!=-1; j++)
				for (k=j+1; k<size(lh) && elts(lh,k)!=-1; k++)
					GPH_addEdge(bond(m), elts(lh,j), elts(lh,k));

			LST_delete(lh);
		}
	}
}

Point_t MOL_seekNormal(Molecule_t* m, unsigned ida, unsigned dad){

	Atom_t* a = atom(m,ida);

	if (ligands(a) == 1)
		return MOL_seekNormal(m, neighbor(a,0), ida);


	if (steric(a) == 2) {
		if (neighbor(a,0) != dad)
			return MOL_seekNormal(m, neighbor(a,0), ida);
		else
			return MOL_seekNormal(m, neighbor(a,1), ida);
	}

	return planNormal(coords(a),
					coords(atom(m,neighbor(a,0))),
					coords(atom(m,neighbor(a,1))));
}

/*
*	Fonction qui crée et initialise une nouvelle molécule.
*
*	@param 	size 	Nombre de sommets de la molécule.
* @return				Nouvelle molécule.
*/
Molecule_t* MOL_create(unsigned size) {
	int i;
	Molecule_t *m = malloc(sizeof(Molecule_t));

	size(m) = size;
	m->atoms = malloc(size*sizeof(Atom_t));

	for (i=0; i<size; i++)
		MOL_createAtom(atom(m,i));

	m->cycle = NULL;
	m->bond = GPH_create();

	return m;
}

/*
*	Fonction qui supprime un atome d'une molécule.
* Supprimer un atome revient à supprimer la liste de ses voisins.
*
*	@param 	a 	Atome à supprimer.
*/
void MOC_deleteAtom(Atom_t* a) {

	LST_delete(neighborhood(a));
}

/*
*	Fonction qui supprime une molécule.
*
*	@param 	m 	Molécule à supprimer.
*/
void MOL_delete(Molecule_t* m) {

	int i;
	for (i=0; i<size(m); i++)
		MOC_deleteAtom(atom(m,i));

	free(m->atoms);
	LST_delete(m->cycle);
	GPH_delete(m->bond);
	free(m);
}