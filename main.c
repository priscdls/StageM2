#include "structure.h"
#include "initialization.h"
#include "expansion.h"
#include "generation.h"
#include "output.h"
#include "utile.h"
#include "assembly.h"

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>

#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define PATHNAME "alphashape.R"

void source (const char* name) {
	SEXP e;

	PROTECT(e = lang2(install("source"), mkString(name)));
  R_tryEval(e, R_GlobalEnv, NULL);
  UNPROTECT(1);
}

int main(int argc, char** argv) {
	
	time_t debut = time(NULL);
	
	/************ Initialisation de l'environnement R ************/
	//char * oldPath;
	int r_argc = 2;

	//oldPath = getenv ("R_HOME");
	setenv ("R_HOME", "/usr/lib/R", 1);

	char *r_argv [] = {"R", "--silent"};
	Rf_initEmbeddedR (r_argc, r_argv);

	source (PATHNAME);

	//Vérification qu'une entrée est passée en paramètre
	if (argc < 2) {
		printf("Veuillez rentrer le nom du fichier de la molécule et/ou l'alpha\n");
		exit(1);
	}

	char* name = argv[1];
	double alpha = atof(argv[2]);

	Main_t* m = MN_create();
	Ashape_t* as3d = NULL;
	
	substrat(m) = initMolecule(name);
	MOL_write(substrat(m));
	//envelope(m) = createShell(substrat(m), alpha);
	envelope(m) = createShell2(substrat(m), alpha, &as3d);

	SHL_write(envelope(m));
	printf("alpha = %0.1f, Nb sommets env = %d\n", alpha, SHL_nbAtom(envelope(m)));
	generationMoc(m);

	//SHL_write(envelope);
	//SHL_write(moc(m,0));
	
	
	/********** Assemblage des motifs **********/
	
	//printf("PRISCILLE\n");
	//alpha = 10;
	assemblage2(name, m, alpha, as3d);
	//testEnveloppe3(m, alpha, as3d);
	//testEnveloppe2(m, alpha);
	
	ASP_delete(as3d);

	/********** Écriture des résultats dans des fichiers **********/

	//Création du dossier de sortie.
	/*char* name = getBasename (name);
  char* dirName = createDir(name);
  char outputname[512];

  printf("Écriture de la molécule %s.\n", name);
  sprintf(outputname, "%s/%s.mol2", dirName, name);
  MOL_writeMol2(outputname, substrat(m));

  printf("Écriture de l'enveloppe.\n");
  sprintf(outputname, "%s/%s_shell.mol2", dirName, name);
  SHL_writeMol2(outputname, envelope(m));

  printf("Écriture de l'enveloppe après insertion des motifs aromatique.\n");
  sprintf(outputname, "%s/%s_aro.mol2", dirName, name);
  SHL_writeMol2(outputname, envarom(m));

  //Suppression des stuctures.
	free(name);
	free(dirName);*/

	output(name, m);
	
	MN_delete(m);

	/************** Fermeture de l'environnement R ***************/
	Rf_endEmbeddedR (0);
	//setenv ("R_HOME", oldPath, 1);
		
	time_t fin = time(NULL);
	long secondes = (long) difftime(fin, debut);
		
	int heure = secondes / 3600;
	secondes -= heure * 3600;
	int minute = secondes / 60;
	secondes -= minute * 60;
	printf("Temps d'execution : %d heure(s) %d minute(s) %ld seconde(s)\n", heure, minute, secondes);
	
	return 0;
}
