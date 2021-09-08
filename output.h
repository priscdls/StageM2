#ifndef __OUTPUT_H
#define __OUTPUT_H

#include "structure.h"

char* createDir(char *);
char* getBasename (char *);
void LST_write(List_t*);
void MOL_write(Molecule_t*);
void SHL_write(Shell_t*);
void GPH_write(Graph_t*) ;
void MOL_writeMol2(char*, Molecule_t*);
void SHL_writeMol2(char*, Shell_t*);
void output(char*, Main_t*);
void outputShell(char* InputFile, Shell_t* s);

#endif
