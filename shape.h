#ifndef SHAPE_H
#define SHAPE_H

#include "bitmap.h"

// struct origin_vector
// {
//     double vector[3] = {0, 0, 0};
// };


struct triangle
{
    // double origin[3]; // coordinates to start drawing from
    // double sides[3][3]; // vectors drawn in order
    int side_colors[3][3] = {{255, 255, 255}, {255, 255, 255}, {255, 255, 255}};
    // double *origin_vector = nullptr; // dependent vector
    unsigned int origin_vector_index = 0;
    unsigned int side_vector_index = 0;
    Pixel *bitmap = nullptr;
};

struct texture
{
    double *origin; // coordinates to start drawing from
    double *ref_x_vector = nullptr;
    double *ref_y_vector = nullptr;
    Pixel *bitmap = nullptr;
    double x_size;
    double y_size;
    bool invert_x = false;
    bool invert_y = false;
    unsigned int origin_vector_index = 0;
    unsigned int ref_x_triangle = 0;
    unsigned int ref_x_side = 0;
    unsigned int ref_x_axis = 0;
    unsigned int ref_y_triangle = 0;
    unsigned int ref_y_side = 0;
    unsigned int ref_y_axis = 0;
    unsigned int ref_triangle_index = 0;
    triangle *ref_triangle = nullptr;

    double x_vector[3] = {0, 0, 0}; // x-axis
    double y_vector[3] = {0, 0, 0}; // y-axis
};

class Shape {
    public:
        unsigned int triangle_count;
        unsigned int texture_count;
        unsigned int vertex_count;

        triangle *triangles;
        texture *textures;
        // origin_vector *vertices;
        // double (*vertices)[3];
        double origin[3] = { 0.0, 0.0, 0.0 };
        double (*vector_array)[3];

        double mass;
        double velocity[3] = { 0.0, 0.0, 0.0 }; // unit?
        double angular_velocity[3] = { 0.0, 0.0, 0.0 };

        Shape(unsigned int triangle_count, unsigned int texture_count, unsigned int vertex_count);

        Shape(unsigned int triangle_count, unsigned int texture_count, unsigned int vertex_count, double mass);

        Shape(const Shape &shape);

        Shape& operator=(const Shape &shape);

        ~Shape();
};

void assemble_cube(Shape *cube, double *origin_coord, double side_len, double *velocity, double *angular_velocity);

void assemble_cube45(Shape *cube, double *origin_coord, double side_len, double *velocity, double *angular_velocity);

#endif //SHAPE_H
