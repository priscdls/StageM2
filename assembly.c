#include "assembly.h"
#include "interface.h"
#include "expansion.h"
#include "output.h"
#include <float.h>

#define NB_MOTIF 5
#define MIN_DIST 2

void affichage(Shell_t* s) {
	
	for (int i = 0; i < SHL_nbAtom(s); i++)
	{
		printf("Atome %d : ", i/*+1*/);
		for (int j = 0; j < neighborhoodSize(atom(s,i)); j++)
		{
			printf("%d ",neighbor(atom(s,i),j)/*+1*/);
		}
		printf("flag %d\n",flag(atom(s,i)));
	}
}

// Vérifie si le point passé en argument est dans l'enveloppe
// Retourne 1 si dans l'enveloppe et 0 si en dehors
int inAShape(Ashape_t* as3d, Point_t p) {
	
	double* point = malloc(3* sizeof(double));
	
	point[0] = p.x;
	point[1] = p.y;
	point[2] = p.z;
	
	int* res = Cinashape3d(as3d, point, 3);
	//printf("InAShape : %d\n", res[0]);
	
	free(point);
	return res[0];
}

/**************************************/
/* Ajout des motifs *******************/
/**************************************/

// Donne le type de l'atome inseré
int typeInsert(int numMotif){
	if (numMotif == 0) // Oxygene
	{
		return 1;
	}
	else if (numMotif == 1) // Azote
	{
		return 3;
	}
	else //Carbone
	{
		return 4;
	}
}

// Ajout du motif 4 (cycle aromatique) perpendiculaire au plan avec son voisin
void ajoutMotif4(Shell_t* mocTraite, List_m* mocAtt, int depart, List_d* nvDepart, Point_t posNvDprt, Ashape_t* as3d) {
	
	for (int i = 0; i < neighborhoodSize(atom(mocTraite,depart)); i++) // Pour tous les plans possibles avec les voisins
	{
		if (neighbor(atom(mocTraite, depart), i) != -1)
		{
			int idSuiv, idCycle;
			//int inASh = 0;
			Point_t posSuiv;
			Point_t positionNvDprt = posNvDprt;
			Shell_t* moc = SHL_copy(mocTraite);
			
			Point_t v1 = coords(atom(moc, neighbor(atom(moc, depart), i)));
			Point_t dpt = coords(atom(moc, depart));
			
			// Ajouter le premier atome du cycle dont on a déjà la position
			int id = SHL_addAtom(moc, positionNvDprt, -1);
			flag(atom(moc, id)) = 4;
			SHL_addEdge(moc, depart, id);
			
			// Cherche la normal pour positionné le cycle
			Point_t normal = planNormal(positionNvDprt, dpt, v1);
			normal = rotation(normalization(vector(positionNvDprt, dpt), 1),  90, normal); // Perpendiculaire
			
			// Positionner les autres atomes du cycle
			v1 = AX1E2(positionNvDprt, coords(atom(moc, depart)), normal, SIMPLE); // Voisin
			positionNvDprt = AX2E1(positionNvDprt, coords(atom(moc, depart)), v1, SIMPLE); 
			/*if (inAShape(as3d, positionNvDprt) == 1) // Si le point est dans l'enveloppe
			{
				inASh++;
			}*/

			idCycle = SHL_addAtom(moc, positionNvDprt, -1);
			flag(atom(moc, idCycle)) = 4;
			SHL_addEdge(moc, id, idCycle);
			
			int idVoisin = idCycle;
			for (int i = 0; i < 4; i++) 
			{
				v1 = coords(atom(moc, neighbor(atom(moc, idVoisin), 0)));
				positionNvDprt = AX1E2(positionNvDprt, v1, normal, SIMPLE);
				
				/*if (inAShape(as3d, positionNvDprt) == 1) // Si le point est dans l'enveloppe
				{
					inASh++;
				}*/
				idCycle = SHL_addAtom(moc, positionNvDprt, -1);
				flag(atom(moc, idCycle)) = 4;
				SHL_addEdge(moc, idVoisin, idCycle);
				
				if ( i == 1 ) // Position du prochain depart pour continuer le chemin
				{
					posSuiv = positionNvDprt;
					idSuiv = idCycle;
				}
				
				idVoisin = idCycle;
			}
			
			SHL_addEdge(moc, id, idCycle);
			
			// Positionner atome qui suit le cycle
			v1 = coords(atom(moc, neighbor(atom(moc, idSuiv), 0)));
			Point_t v2 = coords(atom(moc, neighbor(atom(moc, idSuiv), 1)));
			positionNvDprt = AX2E1(posSuiv, v1, v2, SIMPLE); 
			/*if (inAShape(as3d, positionNvDprt) == 1) // Si le point est dans l'enveloppe
			{
				inASh++;
			}*/
			int idSuiv2 = SHL_addAtom(moc, positionNvDprt, -1);
			flag(atom(moc, idCycle)) = 4;
			SHL_addEdge(moc, idSuiv, idSuiv2);
			
			//if (inASh == 0) // Si tous les atomes du motif sont hors de l'enveloppe
			//{
				LSTm_addElement(mocAtt, moc);
				LSTd_addElement(nvDepart, idSuiv2);
			/*}
			else
			{
				SHL_delete(moc);
			}*/
			
		}
	}
}

/*void aff(Shell_t* mocTraite, Point_t point) {
	
	Point_t a = point;
	a.x += 0.5;
	int id = SHL_addAtom(mocTraite, a, -1);
	flag(atom(mocTraite,id)) = 2;
}*/

// Ajout seulement du 0 sur un C d'un motif 3
List_m* ajoutOMotif3(Shell_t* mocTraite, int depart, Ashape_t* as3d) {
	
	List_m* mAtt = LSTm_init();
	Point_t dpt = coords(atom(mocTraite, depart));
	int voisin1 = neighbor(atom(mocTraite, depart), 0); // Voisin
	Point_t v1 = coords(atom(mocTraite, voisin1));
	
	for (int i = 0; i < neighborhoodSize(atom(mocTraite,voisin1)); i++) // Pour tous les plans possibles avec les voisins
	{
		if (neighbor(atom(mocTraite, voisin1), i) != -1 && neighbor(atom(mocTraite, voisin1), i) != depart)
		{
			Shell_t* moc = SHL_copy(mocTraite);
			Shell_t* moc2 = SHL_copy(mocTraite);
			
			Point_t v2 = coords(atom(moc, neighbor(atom(moc, voisin1), i)));
						
			// Cherche la normal pour positionné le O
			Point_t normal = planNormal(dpt, v1, v2);
			
			// Atome O
			//Premier position
			Point_t positionO = AX1E2(dpt, v1, normal, SIMPLE);
			
			if(inAShape(as3d, positionO) == 0) // Si le point n'est pas dans l'enveloppe
			{
				int id3 = SHL_addAtom(moc, positionO, -1);
				flag(atom(moc, id3)) = 1;
				SHL_addEdge(moc, depart, id3);
				
				LSTm_addElement(mAtt, moc);
			}
			else
			{
				SHL_delete(moc);
			}
						
			//Seconde position
			positionO = AX2E1(dpt, v1, positionO, SIMPLE);
			
			//if (inAShape(as3d, positionO) == 0) // Si le point n'est pas dans l'enveloppe
			//{
				int id4 = SHL_addAtom(moc2, positionO, -1);
				flag(atom(moc2, id4)) = 1;
				SHL_addEdge(moc2, depart, id4);
				
				LSTm_addElement(mAtt, moc2);
			/*}
			else
			{
				SHL_delete(moc2);
			}*/
			
		}
	}
	
	return mAtt;
}

// Ajout du motif 3 (C = O) dans le plan avec son voisin
void ajoutMotif3(Shell_t* mocTraite, List_m* mocAtt, int depart, List_d* nvDepart, Point_t positionNvDprt, Ashape_t* as3d) {
	
	for (int i = 0; i < neighborhoodSize(atom(mocTraite,depart)); i++) // Pour tous les plans possibles avec les voisins
	{
		if (neighbor(atom(mocTraite, depart), i) != -1)
		{
			Shell_t* moc = SHL_copy(mocTraite);
			Shell_t* moc2 = SHL_copy(mocTraite);
			
			Point_t v1 = coords(atom(moc, neighbor(atom(moc, depart), i)));
			Point_t dpt = coords(atom(moc, depart));
			
			// Atome C
			int id = SHL_addAtom(moc, positionNvDprt, -1);
			flag(atom(moc, id)) = 4;
			SHL_addEdge(moc, depart, id);
			
			int id2 = SHL_addAtom(moc2, positionNvDprt, -1);
			flag(atom(moc2, id2)) = 4;
			SHL_addEdge(moc2, depart, id2);
			
			// Cherche la normal pour positionné le O
			Point_t normal = planNormal(positionNvDprt, dpt, v1);
			
			// Atome O
			//Premier position
			Point_t positionO = AX1E2(positionNvDprt, dpt, normal, SIMPLE);
			
			//if (inAShape(as3d, positionO) == 0) // Si le point n'est pas dans l'enveloppe
			//{
				int id3 = SHL_addAtom(moc, positionO, -1);
				flag(atom(moc, id3)) = 1;
				SHL_addEdge(moc, id, id3);
				
				LSTm_addElement(mocAtt, moc);
				LSTd_addElement(nvDepart, id);
			/*}
			else
			{
				SHL_delete(moc);
			}*/
			
			//Seconde position
			positionO = AX2E1(positionNvDprt, dpt, positionO, SIMPLE);
			
			//if (inAShape(as3d, positionO) == 0) // Si le point n'est pas dans l'enveloppe
			//{
				int id4 = SHL_addAtom(moc2, positionO, -1);
				flag(atom(moc2, id4)) = 1;
				SHL_addEdge(moc2, id2, id4);
				
				LSTm_addElement(mocAtt, moc2);
				LSTd_addElement(nvDepart, id2);
			/*}
			else
			{
				SHL_delete(moc2);
			}*/
			
		}
	}
}

// Ajout l'atome projeté a l'enveloppe
void ajoutProjection(Shell_t* mocTraite, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Point_t positionNvDprt, Ashape_t* as3d) {
	
	if (numMotif == 3) // C = 0
	{
		ajoutMotif3(mocTraite, mocAtt, depart, nvDepart, positionNvDprt, as3d);
	}
	else if (numMotif == 4)
	{
		ajoutMotif4(mocTraite, mocAtt, depart, nvDepart, positionNvDprt, as3d);
	}
	else
	{
		Shell_t* moc = SHL_copy(mocTraite);
	
		int id = SHL_addAtom(moc, positionNvDprt, -1);
		flag(atom(moc, id)) = typeInsert(numMotif);
		SHL_addEdge(moc, depart, id);
		
		LSTm_addElement(mocAtt, moc);
		LSTd_addElement(nvDepart, id);
		
		// Visualisation
		/*Point_t v = positionNvDprt;
		v.z += 0.5;
		int id2 = SHL_addAtom(moc, v, -1);
		flag(atom(moc, id2)) = 2;
		* */
	}
	
}

/**************************************/
/* Projection emplacement *************/
/**************************************/

// Projection pour un atome avec 1 voisin
void projectionOCN_AX1E3(Shell_t* moc, List_m* mocAtt, int depart, int arrivee, List_d* nvDepart, int numMotif, Ashape_t* as3d) {
	
	List_s* positions = LSTs_init();
	Point_t dpt = coords(atom(moc, depart));
	Point_t arv = coords(atom(moc, arrivee));
	Point_t v1 = coords(atom(moc, neighbor(atom(moc, depart), 0)));
	Point_t x2; // Voisin du voisin
	
	if (neighbor(atom(moc,neighbor(atom(moc,depart),0)),0) == depart)
		x2 = coords(atom(moc,neighbor(atom(moc,neighbor(atom(moc,depart),0)),1)));
	else
		x2 = coords(atom(moc,neighbor(atom(moc,neighbor(atom(moc,depart),0)),0)));

	Point_t normal = planNormal(dpt, v1, x2);
	
	Point_t positionNvDprt = AX1E3(dpt, v1, normal, SIMPLE);
	LSTs_addElement(positions, positionNvDprt);
	
	for (int i = 0; i < 11; i++) // Rotation a 360°
	{
		normal = rotation(normalization(vector(dpt, v1), 1),  30, normal); // Rotation de 30° de la normal
		positionNvDprt = AX1E3(dpt, v1, normal, SIMPLE);
		
		//if (inAShape(as3d, positionNvDprt) == 0) // Si le point n'est pas dans l'enveloppe
		//{
			LSTs_addElement(positions, positionNvDprt);
		//}
	}
	
	for (int i = 0; i < 3 && positions->premier; i++) // 3 positions les mieux placés (distance min avec arrivée)
	{
		positionNvDprt = distMin(positions, arv); 
		LSTs_removeElement(positions, positionNvDprt);
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt, as3d); // Ajout a l'enveloppe
	}
	
	LSTs_delete(positions);
}

// Projection pour un azote avec 2 voisins
void projectionN_AX2E2(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Ashape_t* as3d) {
	
	Point_t positionNvDprt = AX2E2(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), SIMPLE);
	
	//if (inAShape(as3d, positionNvDprt) == 0) // Si le point n'est pas dans l'enveloppe
	//{
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt, as3d); // Ajout a l'enveloppe
	//}
}

// Projection pour un carbone avec 2 voisins dont un oxygene
void projectionC_AX2E1(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Ashape_t* as3d) {
	
	Point_t positionNvDprt = AX2E1(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), SIMPLE);
	
	//if (inAShape(as3d, positionNvDprt) == 0) // Si le point n'est pas dans l'enveloppe
	//{
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt, as3d); // Ajout a l'enveloppe
	//}
}

// Projection pour un carbone avec 2 voisins
void projectionC_AX2E2(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Ashape_t* as3d) {
	
	Point_t positionNvDprt = AX2E2(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), SIMPLE);
	
	//if (inAShape(as3d, positionNvDprt) == 0) // Si le point n'est pas dans l'enveloppe
	//{
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt, as3d); // Ajout a l'enveloppe
	//}
	
	Point_t positionNvDprt2 = AX3E1(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), positionNvDprt, SIMPLE);
	
	//if (inAShape(as3d, positionNvDprt2) == 0) // Si le point n'est pas dans l'enveloppe
	//{
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt2, as3d); // Ajout a l'enveloppe
	//}
}

// Projection pour un carbone avec 3 voisins
void projectionC_AX3E1(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Ashape_t* as3d) {
	
	Point_t positionNvDprt = AX3E1(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), coords(atom(moc, neighbor(atom(moc, depart), 2))), SIMPLE);
	
	//if (inAShape(as3d, positionNvDprt) == 0) // Si le point n'est pas dans l'enveloppe
	//{
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt, as3d); // Ajout a l'enveloppe
	//}
}

/**************************************/
/* Générer chemin *********************/
/**************************************/

// Insertion du motif passé en argument
void insererMotif(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, int arrivee, Ashape_t* as3d){
	
	if ( LST_nbElements(neighborhood(atom(moc, depart))) == 1 ) // Oxygene ou Azote ou Carbone avec 1 voisin
	{
		//Projections
		//Diff rotations
		projectionOCN_AX1E3(moc, mocAtt, depart, arrivee, nvDepart, numMotif, as3d);
	}
	else if (flag(atom(moc, depart)) == 3 && LST_nbElements(neighborhood(atom(moc, depart))) == 2) // Azote avec 2 voisins
	{
		//Projection
		projectionN_AX2E2(moc, mocAtt, depart, nvDepart, numMotif, as3d);
	}
	else if (flag(atom(moc, depart)) == 4) // Carbone
	{
		if (LST_nbElements(neighborhood(atom(moc, depart))) == 2) // 2 voisins
		{
			if (flag(atom(moc, neighbor(atom(moc, depart), 0))) == 1 || flag(atom(moc, neighbor(atom(moc, depart), 1))) == 1) // Si 1 des 2 voisins est un oxygene
			{
				//Projection
				projectionC_AX2E1(moc, mocAtt, depart, nvDepart, numMotif, as3d);
			}
			else
			{
				// 2 Projections
				projectionC_AX2E2(moc, mocAtt, depart, nvDepart, numMotif, as3d);
			}
		}
		else // 3 voisins
		{
			//Projection
			projectionC_AX3E1(moc, mocAtt, depart, nvDepart, numMotif, as3d);
		}
		
	}
}

// Calcule si le nouveau depart est plus loin de l'arrivée que l'ancien départ
int eloigne(Point_t depart, Point_t nvDepart, Point_t arrivee){
	
	float d1 = dist(depart, arrivee);
	float d2 = dist(nvDepart, arrivee);
	
	if (d1 > d2)
	{
		return 0; // L'ancien est plus eloigné : On s'est rapprocher
	}
	else
	{
		return 1; // Le nouveau est plus eloigné
	}
	
}

// Génère le chemin entre 2 groupements de motifs
/*void genererChemin(Molecule_t* mol, List_m* mocAtt, Shell_t* mocTraite, int depart, int arrivee, Elem_d* sommetInter){
	
	for (int i = 0; i < 1; i++) // NB_MOTIF
	{
		List_m* moc = LSTm_init();
		List_d* nvDepart = LSTd_init();
		
		insererMotif(mocTraite, moc, depart, nvDepart, i, sommetInter->sommet);
		
		while (moc->premier)
		{
			if (eloigne( coords(atom(mocTraite, depart)), coords(atom(moc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, arrivee)) )) // Si le nv depart est plus éloigné 
			{
				if (sommetInter->suivant != NULL) // Si ce n'est pas la derniere arrivee
				{
					genererChemin(mol, mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, sommetInter->suivant);
				}
				else // C'est la derniere arrivee
				{
					if (dist( coords(atom(moc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, arrivee)) ) < MIN_DIST) // Proche de l'arrivée 
					{
						SHL_addEdge(moc->premier->moc, nvDepart->premier->sommet, arrivee);// Ajout lien entre dernier sommet du chemin et arrivee
						LSTm_addElement(mocAtt, SHL_copy(moc->premier->moc)); // Ajout dans la liste a traiter
					}
					else
					{
						printf("Modifier angles\n");
						// Modifier angles
						SHL_addEdge(moc->premier->moc, nvDepart->premier->sommet, arrivee);
						LSTm_addElement(mocAtt, SHL_copy(moc->premier->moc)); // Ajout dans la liste a traiter
					}
				}
			}
			else
			{
				genererChemin(mol, mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, sommetInter);
			}
			LSTm_removeFirst(moc);
			LSTd_removeFirst(nvDepart);
		}
	}
}*/

// Génère le chemin entre 2 groupements de motifs
/*void genererChemin2(List_m* mocAtt, Shell_t* mocTraite, int depart, int arrivee, Elem_d* sommetInter, char* InputFile){
	
	for (int i = 0; i < 1; i++) // NB_MOTIF
	{
		List_m* moc = LSTm_init();
		List_d* nvDepart = LSTd_init();
		
		insererMotif(mocTraite, moc, depart, nvDepart, i, sommetInter->sommet);
		
		if (moc->premier) // Pour toutes les solutions générées en generant le chemin
		{
			if (eloigne( coords(atom(mocTraite, depart)), coords(atom(moc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, sommetInter->sommet)) )) // Si le nv depart est plus éloigné 
			{
				if (sommetInter->suivant != NULL) // Si ce n'est pas la derniere arrivee
				{
					genererChemin2(mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, sommetInter->suivant, InputFile);
				}
				else // C'est la derniere arrivee
				{
					if (dist( coords(atom(moc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, arrivee)) ) < MIN_DIST) // Proche de l'arrivée 
					{
						//printf("Modifier angles\n");
						// Modifier angles
						SHL_addEdge(moc->premier->moc, nvDepart->premier->sommet, arrivee);// Ajout lien entre dernier sommet du chemin et arrivee
						LSTm_addElement(mocAtt, SHL_copy(moc->premier->moc)); // Ajout dans la liste a traiter
						printf("AJOUT\n");
						outputShell(InputFile, moc->premier->moc); // A RETIRER
					}*/
					/*else
					{
						//printf("Modifier angles\n");
						// Modifier angles
						SHL_addEdge(moc->premier->moc, nvDepart->premier->sommet, arrivee);
						LSTm_addElement(mocAtt, SHL_copy(moc->premier->moc)); // Ajout dans la liste a traiter
						printf("AJOUT2\n");
					}*/
				/*}
			}
			else
			{
				genererChemin2(mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, sommetInter, InputFile);
			}
			LSTm_removeFirst(moc);
			LSTd_removeFirst(nvDepart);
		}
		LSTm_delete(moc);
		LSTd_delete(nvDepart);
	}
}*/

// Génère le chemin entre 2 groupements de motifs
// Sans sommets intermediaires
void genererChemin3(Main_t* m, List_m* mocAtt, Shell_t* mocTraite, int depart, int arrivee, int nbMotif3, int nbMotif4, char* InputFile, Ashape_t* as3d){
	
	// Visualisation
	/*Point_t v = coords(atom(mocTraite, depart));
	v.z += 0.5;
	int id2 = SHL_addAtom(mocTraite, v, -1);
	flag(atom(mocTraite, id2)) = 2;*/
	
	for (int i = 0; i < NB_MOTIF; i++)
	{
		List_m* moc = LSTm_init();
		List_d* nvDepart = LSTd_init();
		
		insererMotif(mocTraite, moc, depart, nvDepart, i, arrivee, as3d);
		
		if (moc->premier) // Pour toutes les solutions générées en générant le chemin / Diff rotations
		{
			// Compte le nombre de motif 3 d'affilée (C = O)
			if (i == 3)
			{
				nbMotif3++;
			}
			else
			{
				nbMotif3 = 0;
			}
			
			// Compte le nombre de motif 4 (Cycle) 
			if (i == 4) 
			{
				nbMotif4++;
			}
			
			if (dist( coords(atom(moc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, arrivee)) ) < MIN_DIST) // Proche de l'arrivée 
			{
				//printf("Modifier angles\n");
				// Modifier angles
				SHL_addEdge(moc->premier->moc, nvDepart->premier->sommet, arrivee); // Ajout lien entre dernier sommet du chemin et arrivee
				LSTm_addElement(mocAtt, SHL_copy(moc->premier->moc)); // Ajout dans la liste a traiter
				
				//printf("AJOUT\n");
				//outputShell(InputFile, moc->premier->moc); // A RETIRER
			}
			else if (nbMotif3 < 4 && nbMotif4 < 2) // Maximum 4 motifs 3 d'affilée et 2 motifs 4 en tout
			{
				genererChemin3(m, mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, nbMotif3, nbMotif4, InputFile, as3d);
			}
			else // A retirer
			{
				//printf("SORTIE\n");
				//outputShell(InputFile, moc->premier->moc);
			}
			
			LSTm_removeFirst(moc);
			LSTd_removeFirst(nvDepart);
		}
		LSTm_delete(moc);
		LSTd_delete(nvDepart);
	}
}

/***********************************************/
/* Chemins avec sommets intermédiaires *********/
/***********************************************/

void initDijkstra(Shell_t* s, int depart, int arrivee, float** dist, int** predecesseur, List_d** Q) {
	
	*dist = malloc(sizeof(float) * SHL_nbAtom(s));
	
	for (int i = 0; i < SHL_nbAtom(s); i++)
	{
		(*dist)[i] = FLT_MAX;
	}
	(*dist)[depart] = 0;
	
	*predecesseur = malloc(sizeof(int) * SHL_nbAtom(s));
	for (int i = 0; i < SHL_nbAtom(s); i++)
	{
		(*predecesseur)[i] = -1;
	}
	
	*Q = LSTd_init();
	for (int i = 0; i < SHL_nbAtom(s) ; i++)
	{
		if (flag(atom(s, i)) == 0 || depart == i || arrivee == i)
		{
			LSTd_addElement(*Q, i);
		}
	}
}

int trouveMin(float* dist, List_d* Q) {
	float mini = FLT_MAX;
	int sommet = -5;
	Elem_d* cursor = Q->premier;
	while (cursor)
	{		
		if (dist[cursor->sommet] <= mini)
		{
			mini = dist[cursor->sommet];
			sommet = cursor->sommet;
		}
		cursor = cursor->suivant;
	}
	
	return sommet;
}

// Verifie si le point passer en argument est trop interne
// dans ce cas a enlever de la liste
int ptInterne(Shell_t* envelope2, Point_t p) {
	
	for (int i = 0; i < SHL_nbAtom(envelope2) ; i++)
	{
		if (p.x == atomX(atom(envelope2, i)) && p.y == atomY(atom(envelope2, i)) && p.z == atomZ(atom(envelope2, i)))
		{
			return 0; // Il est sur l'enveloppe grossiere donc pas interne
		}
	}
	return 1; // Il n'est pas sur l'enveloppe grossiere donc trop interne
} 

void majDistances(Shell_t* envelope2, Shell_t* s, float* dist, int* predecesseur, int s1, int s2, int arrivee) {
	float poids = PT_distance(coords(atom(s, s1)), coords(atom(s, s2)));
	if ( ptInterne(envelope2, coords(atom(s,s2))) && s1!=arrivee && s2!=arrivee ) // Si le point est trop interne le passer
	{
		poids += 1000;
	}
	
	if (dist[s2] > dist[s1] + poids)
	{
		dist[s2] = dist[s1] + poids;
		predecesseur[s2] = s1;
	}
}

// Plus court chemin
int* dijkstra(Shell_t* envelope2, Shell_t* s, int depart, int arrivee) {
	float* dist = NULL;
	int* predecesseur = NULL;
	List_d* Q = NULL;
	
	initDijkstra(s, depart, arrivee, &dist, &predecesseur, &Q);
		
	while (Q->premier) 
	{
		int s1 = trouveMin(dist, Q);
		LSTd_removeSommet(Q, s1);
		
		for (int i = 0; i < LST_nbElements(neighborhood(atom(s,s1))); i++)
		{
			int s2 = neighbor(atom(s,s1), i);
			majDistances(envelope2, s, dist, predecesseur, s1, s2, arrivee);
		}
	}
	
	free(dist);
	LSTd_delete(Q);
	
	return predecesseur;
}

// Determine les sommets intermédiaires du chemin 
List_d* sommetIntermediaire(Main_t* m, Shell_t* s, int depart, int arrivee) {
	
	double alpha2 = 20.0;
	Shell_t* envelope2 = createShell(substrat(m), alpha2); // Enveloppe grossiere
	
	int* predecesseur = dijkstra(envelope2, s, depart, arrivee);
	
	List_d* sommets = LSTd_init();
	
	int si = arrivee;
	while (si != depart)
	{
		//flag(atom(s, si)) = 2; // Visualisation 
		LSTd_addElement(sommets, si);
		si = predecesseur[si];
	}
	
	free(predecesseur);
	SHL_delete(envelope2);
	
	return sommets;
}

/*************************************************/
/* Sommets départ et arrivée du chemin ***********/
/*************************************************/

// Vérifie si le sommet se situe en bordure de motif
int bordureCheck(Shell_t* s, AtomShl_t* sommet) {
	
	for (int i = 0; i < LST_nbElements(neighborhood(sommet)); i++) // Pour tous les voisins du sommet
	{
		// Si ce sommet est dans un motif et qu'un de ses voisins est de l'enveloppe
		if (flag(atom(s, neighbor(sommet, i))) == 0 && flag(sommet) != 0 )
		{
			return 1; 
		}
	}
	
	return 0;
}

// Parcours en profondeur en fonction des indices sur les sommets des motifs uniquement
int parcours(Shell_t* s, List_t* marquer, int indice1, int indice2) {
	
	AtomShl_t* a = atom(s, indice1);
	LST_addElement(marquer, indice1);
	
	if (neighborhoodSize(a) == 0)
	{
		return 0;
	}
	else
	{
		for (int i = 0; i < neighborhoodSize(a) && neighbor(a, i) != -1; i++) // Pour tous les voisins de a
		{
			if (flag(atom(s, neighbor(a, i))) != 0) // Si le sommet est dans un motif donc de priorité != 0
			{
				if (neighbor(a, i) == indice2) // Si l'identifiant recherché est trouvé
				{
					return 1;
				}
				else
				{
					if (!LST_check(marquer, neighbor(a, i))) // Si l'identifiant de ce sommet n'est pas déjà marqué
					{
						int valide = parcours(s, marquer, neighbor(a, i), indice2);
						if (valide)
						{
							return 1;
						}
					}
				}
			}
			
			
		}
	}
	return 0;
}

// Vérifie s'il existe un chemin passant seulement par des sommets qui appartiennent aux motifs donc de priorité != 0
// donc si les 2 sommets sont du même groupement
int existeChemin(Shell_t* s, int indice1, int indice2){
	
	List_t* marquer = NULL;
	marquer = LST_create();
	
	int existe = parcours(s, marquer, indice1, indice2);
	
	LST_delete(marquer);
	
	return existe;
}

// Génère tous les couples de sommets à relier possible entre des groupements
List_p* choixSommets(Shell_t* s){
	
	List_p* sommets = LST2_init();
	
	for (int i = 0; i < SHL_nbAtom(s) - 1; i++) //Pour tous les sommets en bordure
	{
		if ( bordureCheck(s, atom(s, i)) )
		{
			for (int j = i+1; j < SHL_nbAtom(s); j++)
			{
				if ( bordureCheck(s, atom(s, j)) )
				{
					if (!existeChemin(s, i, j)) // Si i et j sont de groupements différents
					{
						LST2_addElement(sommets, i, j);
					}
				}
			}
		}
	}
	
	return sommets;
}

// Crée une liste de moc a traiter et vide le tableau de solutions finales
List_m* initMocAtt(Main_t* m){
	List_m* mocAtt = LSTm_init();
	
	for (int i=0; i<mocSize(m); i++) //Pour tous les mocs
	{
		if (moc(m,i) != NULL)
		{
			if (i == 0) // Traite juste le premier mocs
			{
				LSTm_addElement(mocAtt, moc(m, i)); // Les mettre dans la liste a traiter
			}
			if (i != 0)
			{
				SHL_delete(moc(m,i)); // Les supprimer du tableau de solutions finales
			}
		}
	}
	
	free(m->mocs);
	m->mocs = NULL;
	mocSize(m) = 0;
	
	return mocAtt;
}

/**************************************/
/* Fonction principale ****************/
/**************************************/

// Fonction principale
/*void assemblage(Main_t* m){
	List_m* mocAtt = initMocAtt(m);
	
	while (mocAtt->premier) // Tant qu'il existe un moc a traiter
	{
		List_p* sommets = choixSommets(mocAtt->premier->moc);
		
		if (!sommets->premier) // S'il n'y qu'un groupement de motifs
		{
			m->mocs[MN_getIndiceFree(m)] = mocAtt->premier->moc; // Ajout au tableau des solutions finales
			LSTm_removeFirst(mocAtt); // Suppression dans la liste a traiter
		}
		else // S'il y a au moins 2 groupements de motifs
		{
			Shell_t* mocTraite = mocAtt->premier->moc; // Copie le moc a traiter
			mocAtt->premier = mocAtt->premier->suivant; // Supprime de la liste à traiter
			
			while (sommets->premier) // Pour tous les couples de sommets à relier
			{
				Shell_t* mocTraite2 = SHL_copy(mocTraite); // Cree un nouveau moc dans la liste a traiter
				printf("33333333");
				int depart = sommets->premier->depart;
				int arrivee = sommets->premier->arrivee;
				printf("444444444");
				LST2_removeFirst(sommets);
				printf("555555555");
					
				// Choix sommets intermédiaire et generer directions
				// genererChemin();
				
				LSTm_addElement(mocAtt, mocTraite2); // Ajout dans la liste a traiter
			}
			SHL_delete(mocTraite);
		}
	}
	free(mocAtt);
}*/

void assemblage2(char* InputFile, Main_t* m, double alpha, Ashape_t* as3d){
	List_m* mocAtt = initMocAtt(m); // ! Prend le premier moc seulement
	
	while (mocAtt->premier) // Tant qu'il existe un moc a traiter
	{
		List_p* sommets = choixSommets(mocAtt->premier->moc);
		
		printf("111111");
		if (!sommets->premier) // S'il n'y a plus qu'un groupement de motifs
		{
			printf("222222");
			// Ecrire directement les solutions
			
			outputShell(InputFile, mocAtt->premier->moc); // Ecriture de la sortie
			LSTm_removeFirst(mocAtt); // Suppression dans la liste a traiter
			
			//mocAtt->premier = mocAtt->premier->suivant; // Suppression dans la liste a traiter
			//m->mocs[MN_getIndiceFree2(m)] = mocAtt->premier->moc; // Ajout au tableau des solutions finales
			
		}
		else // S'il y a au moins 2 groupements de motifs
		{
			Shell_t* mocTraite = mocAtt->premier->moc; // Copie le moc a traiter
			mocAtt->premier = mocAtt->premier->suivant; // Supprime de la liste à traiter
			
			while (sommets->premier) // Pour tous les couples de sommets à relier
			{
				Shell_t* mocTraite2 = SHL_copy(mocTraite); // Cree un nouveau moc dans la liste a traiter
				printf("33333333");
				int depart = sommets->premier->depart;
				int arrivee = sommets->premier->arrivee;
				printf("444444444");
				LST2_removeFirst(sommets);
				printf("555555555");
					
				//List_d* sommetInter = sommetIntermediaire(m, mocTraite2, depart, arrivee); // Choix sommets intermédiaires
				
				printf("666666");

				for (int i = 0; i < LST_nbElements(neighborhood(atom(mocTraite2, depart))); i++) // Retire les voisins enveloppe de l'atome de départ (bordure)
				{
					if (flag(atom(mocTraite2, neighbor(atom(mocTraite2, depart), i))) == 0)
					{
						SHL_removeEdge(mocTraite2, depart, neighbor(atom(mocTraite2, depart), i));
						i--;
					}
				}
								
				for (int i = 0; i < 4 ; i++) // Attribution de tous les types a l'atome de départ (sommet en bordure)
				{
					flag(atom(mocTraite2, depart)) = typeInsert(i);
					
					if (i == 3) 
					{
						if (LST_nbElements(neighborhood(atom(mocTraite2, depart))) == 1) // Motif possible que si le depart n'a qu'un voisin
						{
							List_m* mAtt = ajoutOMotif3(mocTraite2, depart, as3d); // Ajout du O du motif 3 ( C = O )
							
							while (mAtt->premier) // Traiter tous les mocs générés par cet ajout
							{
								//genererChemin2(mocAtt, mAtt->premier->moc, depart, arrivee, sommetInter->premier, InputFile);
								genererChemin3(m, mocAtt, mAtt->premier->moc, depart, arrivee, 0, 0, InputFile, as3d);
								LSTm_removeFirst(mAtt);
							}
							LSTm_delete(mAtt);
						}
					}
					else
					{
						//genererChemin2(mocAtt, mocTraite2, depart, arrivee, sommetInter->premier, InputFile);
						genererChemin3(m, mocAtt, mocTraite2, depart, arrivee, 0, 0, InputFile, as3d);

					}
					
				}
				
				printf("777777777\n");
				
				//LSTd_delete(sommetInter); // Supprime la liste des sommets intermediaires
				SHL_delete(mocTraite2);
			}
			SHL_delete(mocTraite);
		}
		LST2_delete(sommets);
	}
	
	free(mocAtt);
}

// Fonction pour tester la fonction inashape3d
void testEnveloppe2(Main_t* m, double alpha) {
	Shell_t* sh = SHL_create();
	Point_t p = PT_init();
	p.x = -1;
	p.y = -1;
	p.z = -1;
	int id = SHL_addAtom(sh, p, -1);
	
	p.x = 1;
	int id2 = SHL_addAtom(sh, p, -1);
	
	p.y = 1;
	int id3 = SHL_addAtom(sh, p, -1);
	
	p.x = -1;
	int id4 = SHL_addAtom(sh, p, -1);

	p.z = 1;
	int id5 = SHL_addAtom(sh, p, -1);

	p.y = -1;
	int id6 = SHL_addAtom(sh, p, -1);
	
	p.x = 1;
	int id7 = SHL_addAtom(sh, p, -1);
	
	p.y = 1;
	int id8 = SHL_addAtom(sh, p, -1);
	
	p.x = 0;
	p.y = -0.5;
	p.z = 0;
	int id9 = SHL_addAtom(sh, p, -1);

	SHL_addEdge(sh, id, id2);
	SHL_addEdge(sh, id, id4);
	SHL_addEdge(sh, id, id6);
	SHL_addEdge(sh, id2, id3);
	SHL_addEdge(sh, id2, id7);
	SHL_addEdge(sh, id3, id4);
	SHL_addEdge(sh, id3, id8);
	SHL_addEdge(sh, id4, id5);
	SHL_addEdge(sh, id6, id5);
	SHL_addEdge(sh, id6, id7);
	SHL_addEdge(sh, id8, id5);
	SHL_addEdge(sh, id8, id7);
	SHL_addEdge(sh, id9, id7);
	SHL_addEdge(sh, id9, id6);
	SHL_addEdge(sh, id9, id);
	SHL_addEdge(sh, id9, id2);


	Ashape_t* as3d = Cashape3d(sh, alpha); //envelope(m) 1
	
	printf("TRIANGLE %d\n", as3d->nb_triang);
	for (int i = 0; i < as3d->nb_triang/9 ; i++)
	{
		printf("%lf, ", as3d->triang[i]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*2]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*3]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*4]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*5]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*6]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*7]);
		printf("%lf\n", as3d->triang[i+(as3d->nb_triang/9)*8]);
	}
	
	double* point = malloc(2 * 3* sizeof(double));
	// Point dans l'enveloppe
	point[0] = 0;
	point[2] = 0;
	point[4] = 0;
	
	// Point hors de l'enveloppe
	point[1] = 7;
	point[3] = -4;
	point[5] = 2;
	
	p.x = point[0];
	p.y = point[2];
	p.z = point[4];
	Point_t p2 = PT_init();
	p2.x = point[1];
	p2.y = point[3];
	p2.z = point[5];
	
	//Cinashape3d(as3d, point, 6);
	
	int a = inAShape(as3d, p);
	int b = inAShape(as3d, p2);
	
	printf("A : %d B : %d\n", a, b);
	
	free(point);
	ASP_delete(as3d);
	SHL_delete(sh);
}

void testEnveloppe3(Main_t* m, double alpha, Ashape_t* as3d) {
	/*printf("TRIANGLE %d\n", as3d->nb_triang);
	for (int i = 0; i < as3d->nb_triang/9 ; i++)
	{
		printf("%lf, ", as3d->triang[i]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*2]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*3]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*4]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*5]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*6]);
		printf("%lf, ", as3d->triang[i+(as3d->nb_triang/9)*7]);
		printf("%lf\n", as3d->triang[i+(as3d->nb_triang/9)*8]);
	}*/
	
	double* point = malloc(2 * 3* sizeof(double));
	// Point dans l'enveloppe
	/*point[0] = 3.6434;
	point[2] = -2.5436;
	point[4] = 2.6012;
	*/
	point[0] = -3.8818;
	point[2] = 2.8605;
	point[4] = 3.4921;
	// Point hors de l'enveloppe
	point[1] = 7;
	point[3] = -4;
	point[5] = 2;
	
	Point_t p = PT_init();
	p.x = point[0];
	p.y = point[2];
	p.z = point[4];
	Point_t p2 = PT_init();
	p2.x = point[1];
	p2.y = point[3];
	p2.z = point[5];
	
	/*int id = SHL_addAtom(moc(m,0), p, -1);
	flag(atom(moc(m,0),id)) = 3; // Bleue
	id = SHL_addAtom(moc(m,0), p2, -1);
	flag(atom(moc(m,0),id)) = 2; // Jaune
	*/
	//Cinashape3d(as3d, point, 6);
	
	printf("POUR A :\n");
	int a = inAShape(as3d, p);
	
	printf("POUR B :\n");
	int b = inAShape(as3d, p2);
	
	printf("A : %d B : %d\n", a, b);
	
	free(point);
}
