#ifndef UTILS_H
#define UTILS_H

#include "shape.h"


constexpr double PI = 3.1415926535897932384626433832795028841971693993751058209;
constexpr double FOV = 150 * PI / 180.0;
constexpr double MIN_DOUBLE_THRESHOLD = .00001;

double degrees_to_radians(double degrees);

void shape_min_max_coords(Shape* shape, double* x_min_max, double* y_min_max, double* z_min_max);

void triangle_min_max_coords(double *side1, double *side2, double *side3, double *x_min_max, double *y_min_max, double *z_min_max, double *origin_coord, double scale_factor);

void cross(double* input, double px1, double py1, double pz1, double px2, double py2, double pz2);

double dot(double *parr1, double *parr2);

bool same_side(double px1, double py1, double pz1, double px2, double py2, double pz2, double ax, double ay, double az, double bx, double by, double bz);

bool plane_point_in_triangle(double *side1, double *side2, double *side3, double x, double y, double z, double *tri_origin);


#endif //UTILS_H
