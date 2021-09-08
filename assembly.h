#ifndef __ASSEMBLY_H
#define __ASSEMBLY_H

#include "structure.h"
#include "utile.h"

void assemblage(Main_t* m);
void assemblage2(char* InputFile, Main_t* m, double alpha, Ashape_t* as3d);
void testEnveloppe2(Main_t* m, double alpha);
void testEnveloppe3(Main_t* m, double alpha, Ashape_t* as3d);

#endif
