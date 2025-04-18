#include <algorithm>
#include <cmath>
#include "utils.h"
#include "movement.h"

double degrees_to_radians(double degrees)
{
    return degrees * (PI / 180.0);
}

// finds the coords of a shape's bounding box
void shape_min_max_coords(Shape* shape, double* x_min_max, double* y_min_max, double* z_min_max)
{
    for (unsigned int i = 0; i < shape->vertex_count; i++) {
        double *current_vector = shape->vector_array[i];

        x_min_max[0] = min(x_min_max[0], current_vector[0] + shape->origin[0]);
        x_min_max[1] = max(x_min_max[1], current_vector[0] + shape->origin[0]);

        y_min_max[0] = min(y_min_max[0], current_vector[1] + shape->origin[1]);
        y_min_max[1] = max(y_min_max[1], current_vector[1] + shape->origin[1]);

        z_min_max[0] = min(z_min_max[0], current_vector[2] + shape->origin[2]);
        z_min_max[1] = max(z_min_max[1], current_vector[2] + shape->origin[2]);
    }
}

// finds the bounding box of a triangle
void triangle_min_max_coords(double *side1, double *side2, double *side3, double *x_min_max, double *y_min_max, double *z_min_max, double *origin_coord, double scale_factor)
{
    // set initial min/max values to starting coord
    for (unsigned int i = 0; i < 2; i++) {
        x_min_max[i] = origin_coord[0];
        y_min_max[i] = origin_coord[1];
        z_min_max[i] = origin_coord[2];
    }

    double current_coord[3] = { origin_coord[0], origin_coord[1], origin_coord[2] };
    double *tri_side_ptrs[3] = { side1, side2, side3 };

    // side3 not checked; endpoint is assumed to be equal to origin
    for (unsigned int side_index = 0; side_index < 2; side_index++) {
        // find next coord by moving to vector's endpoint
        for (unsigned int i = 0; i < 3; i++) {
            current_coord[i] += tri_side_ptrs[side_index][i] * scale_factor;
        }

        // compare coord to existing value
        x_min_max[0] = min(x_min_max[0], current_coord[0]);
        x_min_max[1] = max(x_min_max[1], current_coord[0]);
        y_min_max[0] = min(y_min_max[0], current_coord[1]);
        y_min_max[1] = max(y_min_max[1], current_coord[1]);
        z_min_max[0] = min(z_min_max[0], current_coord[2]);
        z_min_max[1] = max(z_min_max[1], current_coord[2]);
    }
}

// find cross product of two vectors
void cross(double* input, double px1, double py1, double pz1, double px2, double py2, double pz2) {
    input[0] = py1 * pz2 - pz1 * py2;
    input[1] = pz1 * px2 - px1 * pz2;
    input[2] = px1 * py2 - py1 * px2;
}

// find dot product of two vectors
double dot(double *parr1, double *parr2) {
    return parr1[0] * parr2[0] + parr1[1] * parr2[1] + parr1[2] * parr2[2];
}

// helper function that checks if a point is on the same side of two vectors using the right-hand rule
bool same_side(double px1, double py1, double pz1, double px2, double py2, double pz2, double ax, double ay, double az, double bx, double by, double bz) {
    double cross1[3];
    double cross2[3];
    cross(cross1, bx - ax, by - ay, bz - az, px1 - ax, py1 - ay, pz1 - az);
    cross(cross2, bx - ax, by - ay, bz - az, px2 - ax, py2 - ay, pz2 - az);

    if (dot(cross1, cross2) >= 0) {
        return true;
    }
    return false;
}

bool plane_point_in_triangle(double *side1, double *side2, double *side3, double x, double y, double z, double *tri_origin) {
    // calculate coords of triangle vertices
    double ax = tri_origin[0];
    double ay = tri_origin[1];
    double az = tri_origin[2];

    double bx = ax + side1[0];
    double by = ay + side1[1];
    double bz = az + side1[2];

    double cx = bx + side2[0];
    double cy = by + side2[1];
    double cz = bz + side2[2];

    // check if point is on the same side for each pair of triangle sides
    bool side1_pass = same_side(x, y, z, ax, ay, az, bx, by, bz, cx, cy, cz);
    bool side2_pass = same_side(x, y, z, bx, by, bz, ax, ay, az, cx, cy, cz);
    bool side3_pass = same_side(x, y, z, cx, cy, cz, ax, ay, az, bx, by, bz);

    return side1_pass && side2_pass && side3_pass;

}
