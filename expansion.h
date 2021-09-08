#ifndef __EXPANSION_H
#define __EXPANSION_H

#include "structure.h"

Shell_t* createShell(Molecule_t* m, double alpha);
Shell_t* createShell2(Molecule_t* m, double alpha, Ashape_t** as3d);
void expansion(Molecule_t* m, Shell_t* s);
Ashape_t* alphaShape2(Shell_t* s, double alpha);

#endif
