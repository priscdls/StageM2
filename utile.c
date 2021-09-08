#include "utile.h"

float radianToDegre(float a) {
	return a * 180 / M_PI;
}

float degreToRadian(float a) {
	return a * M_PI / 180;
}

Point_t initPoint(float scal) {
	Point_t _new;

	_new.x = scal;
	_new.y = scal;
	_new.z = scal;

	return _new;
}

Point_t addPoint(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = A.x + B.x;
	_new.y = A.y + B.y;
	_new.z = A.z + B.z;

	return _new;
}

Point_t subPoint(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = A.x - B.x;
	_new.y = A.y - B.y;
	_new.z = A.z - B.z;

	return _new;
}

Point_t mulPoint(Point_t A, float scal) {
	Point_t _new;

	_new.x = scal * A.x;
	_new.y = scal * A.y;
	_new.z = scal * A.z;

	return _new;
}

Point_t divPoint(Point_t A, float scal) {
	Point_t _new;

	_new.x = A.x / scal;
	_new.y = A.y / scal;
	_new.z = A.z / scal;

	return _new;
}

Point_t merPoint(Point_t A, Point_t B) {

	Point_t _new;

	_new.x = (A.x + B.x)/2;
	_new.y = (A.y + B.y)/2;
	_new.z = (A.z + B.z)/2;

	return _new;
}

//Calcul de la distance entre deux points.
float dist(Point_t A, Point_t B) {
	return sqrt(pow((A.x - B.x), 2) + pow((A.y - B.y), 2)	+ pow((A.z - B.z), 2));
}

//Normaliser un vecteur à la longueur length.
Point_t normalization(Point_t normal, float length) {
	///A vérifier
	Point_t a;
	float z;

	z = sqrt(pow(length,2) / (pow(normal.x,2) + pow(normal.y,2) + pow(normal.z,2)));
	a.x = z * normal.x;
	a.y = z * normal.y;
	a.z = z * normal.z;
	return a;
}

float angle(Point_t A, Point_t B, Point_t C) {
	float AB = dist(A,B), AC = dist(A,C), BC = dist(B,C);

	return acos( (pow(AC,2)+pow(AB,2)-pow(BC,2)) / (2*AC*AB) ) * 180 / M_PI;
}

Point_t vector(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = B.x - A.x;
	_new.y = B.y - A.y;
	_new.z = B.z - A.z;

	return _new;
}



//Retourne la normal d'un plan à partir de 3 points du plan.
Point_t planNormal(Point_t A, Point_t B, Point_t C) {
	Point_t normal;

	normal.x = (B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y);
  normal.y = (B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z);
  normal.z = (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);

  return normalization(normal, 1);
}

//Rotation à partir d'un vecteur rotation, d'un angle et d'un point.
//Alpha doit être en degree
Point_t rotation(Point_t vec, float alpha, Point_t A) {
	Point_t rot;
	alpha = degreToRadian(alpha);
	vec = normalization(vec, 1);

	rot.x = (vec.x*vec.x + (1 - vec.x*vec.x)*cos(alpha)) * A.x +
			(vec.x*vec.y * (1 - cos(alpha)) - vec.z * sin(alpha)) * A.y +
			(vec.x*vec.z * (1 - cos(alpha)) + vec.y * sin(alpha)) * A.z;
	rot.y = (vec.y*vec.y + (1 - vec.y*vec.y)*cos(alpha)) * A.y +
			(vec.x*vec.y * (1 - cos(alpha)) + vec.z * sin(alpha)) * A.x +
			(vec.y*vec.z * (1 - cos(alpha)) - vec.x * sin(alpha)) * A.z;
	rot.z = (vec.z*vec.z + (1 - vec.z*vec.z)*cos(alpha)) * A.z +
			(vec.x*vec.z * (1 - cos(alpha)) - vec.y * sin(alpha)) * A.x +
			(vec.y*vec.z * (1 - cos(alpha)) + vec.x * sin(alpha)) * A.y;

	return rot;
}

Point_t autre(Point_t A, Point_t B, Point_t C, float scal) {
	Point_t normal;

	normal.x = 2 * A.x - B.x - C.x;
 	normal.y = 2 * A.y - B.y - C.y;
 	normal.z = 2 * A.z - B.z - C.z;
 	normal = normalization(normal, scal);
  			
  return addPoint(A, normal);
}

Point_t AX1E1(Point_t a, Point_t x1, float length) {

	return addPoint(a, normalization(vector(x1, a), length));
}

Point_t AX2E1(Point_t a, Point_t x1, Point_t x2, float length) {

	Point_t v1, v2;

	v1 = normalization(vector(x1, a), 1);
	v2 = normalization(vector(x2, a), 1);

	return addPoint(a, normalization(addPoint(v1, v2), length));
}

Point_t AX1E2(Point_t a, Point_t x1, Point_t normal, float length) {

	Point_t v1;

	v1 = normalization(vector(a, x1), 1);

	return addPoint(a, normalization(rotation(normal, 120, v1), length));
}

Point_t AX3E1(Point_t a, Point_t x1, Point_t x2, Point_t x3, float length) {

	Point_t v1, v2, v3;

	v1 = normalization(vector(x1, a), 1);
	v2 = normalization(vector(x2, a), 1);
	v3 = normalization(vector(x3, a), 1);

	return addPoint(a, normalization(addPoint(v1, addPoint(v2,v3)), length));
}

Point_t AX2E2(Point_t a, Point_t x1, Point_t x2, float length) {

	Point_t v1, v2, other, normal, zero = {0, 0, 0};
	float angle = 180 - (109.47/2);

	v1 = normalization(vector(a,x1), 1);
	v2 = normalization(vector(a,x2), 1);

	other = normalization(addPoint(v1,v2), 1);
	normal = normalization(planNormal(zero, planNormal(a, x1, x2), other), 1);


	return addPoint(a, normalization(rotation(normal, angle, other), length));
}

Point_t AX1E3(Point_t a, Point_t x1, Point_t normal, float length) {

	Point_t v1;

	v1 = normalization(vector(a, x1), 1);

	return addPoint(a, normalization(rotation(normal, 109.47, v1), length));
}
