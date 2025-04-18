#ifndef DRAW_H
#define DRAW_H
#include "shape.h"

void draw_shape(Shape *shape, unsigned int triangle_index, Pixel *bitmap, double **depth_map);

void draw_texture(Shape *shape, unsigned int tri_index, Pixel *bitmap, double **depth_map);

#endif //DRAW_H
