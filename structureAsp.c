#include "structure.h"

/**************************************/
/* ASHAPE3D ***************************/
/**************************************/

/**
* Fonction qui alloue un nouveau Ashape.
*
*/
Ashape_t* ASP_create() {
	Ashape_t* as3d = malloc(sizeof(Ashape_t));

	as3d->nb_triang = 0;
	as3d->nb_edge = 0;
	as3d->nb_vertex = 0;
	as3d->nb_x = 0;
	as3d->nb_alpha = 0;
	
	as3d->triang = NULL;
	as3d->edge = NULL;
	as3d->vertex = NULL;
	as3d->x = NULL;
	as3d->alpha = NULL;

	return as3d;
}

/**
* Fonction qui libère la mémoire utilisé par l'ashape.
*
* @param		as3d 		Adresse de l'ashape à détruire.
*/
void ASP_delete(Ashape_t* as3d) {

	if (as3d != NULL) {
		free(as3d->triang);
		free(as3d->edge);
		free(as3d->vertex);
		free(as3d->x);
		free(as3d->alpha);
	}

	free(as3d);
}