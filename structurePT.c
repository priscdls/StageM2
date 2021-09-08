#include "structure.h"
#include <math.h>

Point_t PT_init() {

	Point_t _new;

	_new.x = 0;
	_new.y = 0;
	_new.z = 0;

	return _new;
}

Point_t PT_add(Point_t A, Point_t B) {

	Point_t _new;

	_new.x = A.x + B.x;
	_new.y = A.y + B.y;
	_new.z = A.z + B.z;

	return _new;
}

Point_t PT_sub(Point_t A, Point_t B) {

	Point_t _new;

	_new.x = A.x - B.x;
	_new.y = A.y - B.y;
	_new.z = A.z - B.z;

	return _new;
}

Point_t PT_mul(Point_t A, float scal) {

	Point_t _new;

	_new.x = scal * A.x;
	_new.y = scal * A.y;
	_new.z = scal * A.z;

	return _new;
}

Point_t PT_div(Point_t A, float scal) {

	Point_t _new;

	if (scal == 0)
		return PT_init();

	_new.x = A.x / scal;
	_new.y = A.y / scal;
	_new.z = A.z / scal;

	return _new;
}

float PT_distance(Point_t A, Point_t B) {

	return sqrt(pow((A.x - B.x), 2) + pow((A.y - B.y), 2)	+ pow((A.z - B.z), 2));
}
