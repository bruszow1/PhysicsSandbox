#ifndef MOVEMENT_H
#define MOVEMENT_H
#include "shape.h"

void model_movement(Shape **input_shapes, unsigned int shape_count, int x_rotate, int y_rotate, int z_rotate, double initial_timestep, double final_timestep);

void draw_shapes(Pixel *input_bitmap, Shape **input_shapes, unsigned int shape_count);

#endif //MOVEMENT_H
