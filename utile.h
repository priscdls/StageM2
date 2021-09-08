#ifndef __UTILE_H
#define __UTILE_H

#include "structure.h"
#include <math.h>

float radianToDegre(float);
float degreToRadian(float);
Point_t initPoint(float);
Point_t addPoint(Point_t, Point_t);
Point_t subPoint(Point_t, Point_t);
Point_t mulPoint(Point_t, float);
Point_t divPoint(Point_t, float);
Point_t merPoint(Point_t, Point_t);

float dist(Point_t, Point_t);
Point_t normalization(Point_t, float);
float angle(Point_t, Point_t, Point_t);
Point_t vector(Point_t, Point_t);
Point_t autre(Point_t, Point_t, Point_t, float);
Point_t planNormal(Point_t, Point_t, Point_t);
Point_t rotation(Point_t, float, Point_t);

Point_t AX1E1(Point_t, Point_t, float);
Point_t AX2E1(Point_t, Point_t, Point_t, float);
Point_t AX1E2(Point_t, Point_t, Point_t, float);
Point_t AX3E1(Point_t, Point_t, Point_t, Point_t, float);
Point_t AX2E2(Point_t, Point_t, Point_t, float);
Point_t AX1E3(Point_t, Point_t, Point_t, float);

#endif
