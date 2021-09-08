#include "interface.h"
#include <string.h>

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>

/*
* Appel à R.
* @param s 		 Enveloppe possédant déjà un nuage de points.
* @param alpha Paramètre de l'aphashape.
*/
Ashape_t* Cashape3d(Shell_t* s, double alpha) {
	Ashape_t* as3d = ASP_create();
	int i, j, size = SHL_nbAtom(s);
	double* data = malloc(size*3*sizeof(double));

	for (i=0, j=0; i<size(s); i++) {
		if (flag(atom(s,i)) != -1) {
			data[j] = atomX(atom(s,i));
			data[j+size] = atomY(atom(s,i));
			data[j+2*size] = atomZ(atom(s,i));
			j++;		
		}
	}

	//Allocation d'un vecteur R avec la copie du vecteur C.
	SEXP arg;
	PROTECT(arg = allocVector(REALSXP, size*3));
	memcpy (REAL(arg), data, 3*size*sizeof(double));

	SEXP alp2;
	PROTECT (alp2 = allocVector (REALSXP, 1) );
	memcpy (REAL(alp2), &alpha, sizeof (double));

	//Configuration de la fonction R pour l'appel.
	SEXP Rashape3d_call;
	PROTECT (Rashape3d_call = lang3(install("Rashape3d"), arg, alp2));
	//Execution de la fonction R
	int errorOccurred;
	SEXP ret = R_tryEval(Rashape3d_call, R_GlobalEnv, &errorOccurred);

	if (!errorOccurred && length(ret)==5) {
		as3d->nb_triang = length(VECTOR_ELT(ret, 0));
		as3d->triang = malloc(as3d->nb_triang*sizeof(double));
		memcpy(as3d->triang, REAL(VECTOR_ELT(ret,0)), as3d->nb_triang*sizeof(double));

		as3d->nb_edge = length(VECTOR_ELT(ret, 1));
		as3d->edge = malloc(as3d->nb_edge*sizeof(double));
		memcpy(as3d->edge, REAL(VECTOR_ELT(ret,1)), as3d->nb_edge*sizeof(double));

		as3d->nb_vertex = length(VECTOR_ELT(ret, 2));
		as3d->vertex = malloc(as3d->nb_vertex*sizeof(double));
		memcpy(as3d->vertex, REAL(VECTOR_ELT(ret,2)), as3d->nb_vertex*sizeof(double));
	
		as3d->nb_x = length(VECTOR_ELT(ret, 3));
		as3d->x = malloc(as3d->nb_x*sizeof(double));
		memcpy(as3d->x, REAL(VECTOR_ELT(ret,3)), as3d->nb_x*sizeof(double));

		as3d->nb_alpha = length(VECTOR_ELT(ret, 4));
		as3d->alpha = malloc(as3d->nb_alpha*sizeof(double));
		memcpy(as3d->alpha, REAL(VECTOR_ELT(ret,4)), as3d->nb_alpha*sizeof(double));
	}

	UNPROTECT(3);

	free(data);
	return as3d;
}

int* Cinashape3d(Ashape_t* as3d, double* points, int nb_points) {
	SEXP triang;
	PROTECT(triang = allocVector(REALSXP, as3d->nb_triang));
	memcpy (REAL(triang), as3d->triang, as3d->nb_triang*sizeof(double));

	SEXP x;
	PROTECT(x = allocVector(REALSXP, as3d->nb_x));
	memcpy (REAL(x), as3d->x, as3d->nb_x*sizeof(double));

	SEXP alpha;
	PROTECT(alpha = allocVector(REALSXP, as3d->nb_alpha));
	memcpy (REAL(alpha), as3d->alpha, as3d->nb_alpha*sizeof(double));

	SEXP point;
	PROTECT(point = allocVector(REALSXP, nb_points));
	memcpy (REAL(point), points, nb_points*sizeof(double));

	//Configuration de la fonction R pour l'appel.
	SEXP Rinashape3d_call;
	PROTECT (Rinashape3d_call = lang5(install("Rinashape3d"), triang, x, alpha, point));
	//Execution de la fonction R
	int errorOccurred;
	SEXP ret = R_tryEval(Rinashape3d_call, R_GlobalEnv, &errorOccurred);

	if (!errorOccurred) {
		/*for (int i = 0; i < length(ret) ; i++)
		{
			printf("SHAPE : %d %d\n", LOGICAL(ret)[i], length(ret));
		}*/
	}

	UNPROTECT(5);
	return LOGICAL(ret);
}
