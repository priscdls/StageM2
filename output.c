#include "output.h"

#include <sys/stat.h>
#include <string.h>

/**
*/
char* createDir(char *input) {
  char* dirName = malloc (256 * sizeof(char));
  
  mkdir("Results", 0755);
  sprintf(dirName,"Results/%s",input);
  mkdir(dirName,0755);

  return dirName;
}

/**
*/
char* getBasename (char * in) {
  char* r = malloc (256 * sizeof (char) );  
  char *start;
  int i;

  start = strrchr(in, '/');

  for (i=0; start[i+1] != '.'; ++i) {
    r[i] = start[i+1];
  }
  r[i] = '\0';
  return r;  
}

/**
*/
void copytoDir(char* InputFile, char* dirName, char* name) {
  char cmd[512];

  sprintf(cmd, "cp %s %s/%s.xyz", InputFile, dirName, name);
  if (system(cmd) == -1) {
    printf("%s",cmd);
    exit(0);
  }
}

void LST_write(List_t* l) {

  int i;

  printf("Size : %d, NbElement : %d\n", size(l), LST_nbElements(l));

  for (i=0; i<size(l); i++)
    printf(" %d,", elts(l,i));
  printf("\n");
}

void MOL_writeAtom(Atom_t* a, unsigned cycle) {
  int i;

  printf("%s (%d, %d, %d) (%8.4f, %8.4f, %8.4f) :", symbol(a),
      ligands(a), lonePairs(a), cycle,
      atomX(a), atomY(a), atomZ(a));

  for (i=0; i<neighborhoodSize(a); i++)
    if (neighbor(a,i) != -1)
      printf(" %d,", neighbor(a,i));

  printf("\n");
}

void MOL_write(Molecule_t* m) {
  int i;

  printf("Size : %d, Edges : %d\n", size(m), MOL_nbEdges(m));
  for (i=0; i<size(m); i++) {
    printf("%d ", i);
    MOL_writeAtom(atom(m,i), cycle(m, i));
  }
  printf("\n");
}

void SHL_writeAtom(AtomShl_t* a) {
  int i;

  printf("(%d) (%d) (%8.4f, %8.4f, %8.4f) :", parentAtom(a), flag(a),
      atomX(a), atomY(a), atomZ(a));

  for (i=0; i<neighborhoodSize(a); i++)
    //if (neighbor(a,i) != -1)
      printf(" %d,", neighbor(a,i));

  printf("\n");
}

void SHL_write(Shell_t* s) {
  int i;

  printf("Size : %d, SizeMax : %d\n", SHL_nbAtom(s), size(s));
  for (i=0; i<size(s); i++) {
    //if (flag(atom(s,i)) != -1) {
      printf("%d ", i);
      SHL_writeAtom(atom(s,i));
    //}
  }
  printf("\n");
}

void GPH_writeVertex(Vertex_t* v) {

  int i;

  for (i=0; i<nbNeighbors(v); i++) {
    if (neighbor(v, i) != -1)
      printf(" %d,", neighbor(v,i));
  }

  printf("\n");
}

void GPH_write(Graph_t* g) {

  int i;
  printf("Size : %d, SizeMax : %d\n", GPH_nbVertex(g), size(g));
  for (i=0; i<size(g); i++) {

    if (id(vertex(g,i)) != -1) {

      printf("%d :", id(vertex(g,i)));
      GPH_writeVertex(vertex(g,i));
    }
  }
}

void MOL_writeMol2(char* output, Molecule_t* m) {
  FILE* filestream = NULL;
  int ret, i, j, l;

  filestream = fopen(output, "w");
  ret = fprintf(filestream, "@<TRIPOS>MOLECULE\n*****\n");
  ret = fprintf(filestream, " %d %d 0 0 0\n", size(m), MOL_nbEdges(m));
  ret = fprintf(filestream, "SMALL\nGASTEIGER\n\n");

  //Ecriture des sommets
  ret = fprintf(filestream, "@<TRIPOS>ATOM\n");
  for (i=0; i<size(m); i++) {
    ret = fprintf(filestream, " %3d %s", i+1, symbol(atom(m,i)));


    if (atomX(atom(m,i)) < 0)
      ret = fprintf(filestream, "   %3.4f", atomX(atom(m,i)));
    else
      ret = fprintf(filestream, "    %3.4f", atomX(atom(m,i)));
    if (atomY(atom(m,i)) < 0)
      ret = fprintf(filestream, "   %3.4f", atomY(atom(m,i)));
    else
      ret = fprintf(filestream, "    %3.4f", atomY(atom(m,i)));
    if (atomZ(atom(m,i)) < 0)
      ret = fprintf(filestream, "   %3.4f", atomZ(atom(m,i)));
    else
      ret = fprintf(filestream, "    %3.4f", atomZ(atom(m,i)));

    ret = fprintf(filestream, " %s\n", symbol(atom(m,i)));
  }

  //Ecriture des liens
  ret = fprintf(filestream, "@<TRIPOS>BOND\n");
  for (i=0, l=1; i<size(m); i++)
    for (j=0; j<ligands(atom(m,i)); j++)
      if (i < neighbor(atom(m,i), j)) {
        ret = fprintf(filestream, " %3d %3d %3d %3d\n", l, i+1, neighbor(atom(m,i),j)+1, 1);
      }

  if (ret<0) {
    printf("L'écriture du fichier %s ne s'est pas bien passé.\n", output);
    exit(2);
  }

}

void SHL_writeMol2(char* output, Shell_t* s) {
  FILE* filestream = NULL;
  int ret, i, j, l;
  int* indice = malloc(size(s)*sizeof(int));

  filestream = fopen(output, "w");
  ret = fprintf(filestream, "@<TRIPOS>MOLECULE\n*****\n");
  ret = fprintf(filestream, " %d %d 0 0 0\n", SHL_nbAtom(s), SHL_nbEdges(s));
  ret = fprintf(filestream, "SMALL\nGASTEIGER\n\n");

  //Ecriture des sommets
  ret = fprintf(filestream, "@<TRIPOS>ATOM\n");
  for (i=0, j=1; i<size(s); i++) {

    if (flag(atom(s,i)) != -1) {
      indice[i] = j;
      if (flag(atom(s,i)) == 2)
        ret = fprintf(filestream, " %3d S", j);
      else if (flag(atom(s,i)) == 3)
        ret = fprintf(filestream, " %3d N", j);
      else if (flag(atom(s,i)) == 1)
        ret = fprintf(filestream, " %3d O", j);
      else 
        ret = fprintf(filestream, " %3d C", j);
      ret = fprintf(filestream, "   %3.4f", atomX(atom(s,i)));
      ret = fprintf(filestream, "   %3.4f", atomY(atom(s,i)));
      ret = fprintf(filestream, "   %3.4f", atomZ(atom(s,i)));
      if (flag(atom(s,i)) == 2)
        ret = fprintf(filestream, "   S\n");
      else if (flag(atom(s,i)) == 3)
        ret = fprintf(filestream, "   N\n");
      else if (flag(atom(s,i)) == 1)
        ret = fprintf(filestream, "   O\n");
      else
        ret = fprintf(filestream, "   C\n");
      j++;
    }
    else
      indice[i] = -1;
  }

  //Ecriture des liens
  ret = fprintf(filestream, "\n@<TRIPOS>BOND\n");
  for (i=0, l=1; i<size(s); i++)
    for (j=0; j<neighborhoodSize(atom(s,i)); j++)
      if (i < neighbor(atom(s,i), j)) {
        ret = fprintf(filestream, " %3d %3d %3d %3d\n",
          l, indice[i], indice[neighbor(atom(s,i),j)], 1);
        l++;
      }

  if (ret<0) {
    printf("L'écriture du fichier %s ne s'est pas bien passé.\n", output);
    exit(2);
  }

  free(indice);
}

void outputShell(char* InputFile, Shell_t* s) {
	char outputname[512];
	char* name = getBasename (InputFile);
	char* dirName = createDir(name);
	static int i = 0;
	
	//Sortie avec enveloppe
	sprintf(outputname, "%s/%s_moc%d.mol2", dirName, name, i);
	SHL_writeMol2(outputname, s);
	
	//Sortie sans enveloppe
	Shell_t* s2 = SHL_copy(s);
	for (int j = 0; j < SHL_nbAtom(s2); j++)
	{
		if (flag(atom(s2,j)) == 0)
		{
			SHL_removeAtom(s2, j); // Enleve les atomes de l'enveloppe qui ne sont pas des motifs
		}
		
	}
	
	sprintf(outputname, "%s/%s_mot%d.mol2", dirName, name, i);
	SHL_writeMol2(outputname, s2);
	printf("Result : %d\n", i);
	SHL_delete(s2);
	free(name);
	free(dirName);
	
	i++;
}

void output(char* InputFile, Main_t* m) {
  char outputname[512];
  char* name = getBasename (InputFile);
  char* dirName = createDir(name);
  //int i;

  printf("taille de %s = %d\n", dirName, (int)strlen(dirName));


  //copytoDir(InputFile, dirName, name);

  sprintf(outputname, "%s/%s.mol2", dirName, name);
  MOL_writeMol2(outputname, substrat(m));

  sprintf(outputname, "%s/%s_shell.mol2", dirName, name);
  SHL_writeMol2(outputname, envelope(m));

  sprintf(outputname, "%s/%s_aro.mol2", dirName, name);
  SHL_writeMol2(outputname, envarom(m));

 //for (i=0; i</*mocSize(m)*/1; i++) { // Neutralisé pour avoir une sortie adaptée à une cage finie
  //printf("GDD de %d\n", i);
  //GPH_write(bond(moc(m,i)));
  //sprintf(outputname, "%s/%s_moc%d.mol2", dirName, name, i);
  //SHL_writeMol2(outputname, moc(m,i));
//}

  //Ecrire les générés
  /*Shell_t* out = SHL_avoir(moc);
  sprintf(outputname, "%s/%s_moc2.mol2", dirName, name);
  SHL_writeMol2(outputname, out);


  SHL_delete(out);*/

  free(name);
  free(dirName);
}
