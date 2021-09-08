#include "input.h"
#include <string.h>

/**************************************/
/* INITIALISATION MOLÉCULE ************/
/**************************************/

/**
* Fonction qui récupère le fichier contenant les données de la molécule.
* Doit avoir une extention .xyz.
*
* @param 		inputname		Nom du fichier contenant les données de la molécule.
* @param 		m 					Adresse de stockage de la molécule.
*/
Molecule_t* readInput_xyz(char* inputname) {
	FILE* filestream = NULL;
	int i, size, ret;
	Molecule_t* m;
	
	filestream = fopen(inputname, "r");
	
	ret = fscanf(filestream, "%d", &size);
	m = MOL_create(size);
	
	for (i=0; i<size(m); ++i) {
		ret = fscanf(filestream, "%s %f %f %f", symbol(atom(m,i)), 
			&atomX(atom(m,i)), &atomY(atom(m,i)),	&atomZ(atom(m,i)));
	}
	
	fclose(filestream);
	
	if (ret<0) {
		printf("La lecture du fichier ne s'est pas bien passé.\n");
		exit(2);
	}

	return m;
}

/**
* Fonction qui récupère les rayons de covalence des atomes.
* Les rayons de covalence des atomes sont contenus dans le fichier rdc.dat.
*
* @param		m 					Adresse de la molécule.
*/
void readCovalence(Molecule_t* m) {
  FILE* filestream = NULL;
  int i, j, number, ret;
  Atom_t* atoms;
  
  filestream = fopen("rdc.dat", "r");
  
  ret = fscanf(filestream, "%d", &number);
  atoms = malloc(number*sizeof(Atom_t));

  //Récupére les rayons de covalence des atomes inscrits dans le fichier.
  for (i=0; i<number; ++i) {
    ret = fscanf(filestream, "%s %u", symbol(atoms+i), &radius(atoms+i));
  }
  
  for (i=0; i<size(m); ++i) {
    //Trouve le symbole qui lui correspond dans les covalence
    for (j=0; strcmp(symbol(atom(m,i)), symbol(atoms+j)) && j<number; ++j);

    //Vérifie si son symbole est référencé dans les covalence du fichier
    if (number <= j) {
      printf("%s n'est pas référencé.\n", symbol(atom(m,i)));
      exit(3);
    }
    else {
      radius(atom(m,i)) = radius(atoms+j);
    }
  }

  fclose(filestream);

  free(atoms);

  if (ret<0) {
    printf("La lecture du fichier 'covalence' ne s'est pas bien passé.\n");
    exit(2);
  }
}
