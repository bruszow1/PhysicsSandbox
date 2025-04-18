#include "shape.h"

Shape::Shape(unsigned int triangle_count, unsigned int texture_count, unsigned int vertex_count, double mass)
{
    this->triangle_count = triangle_count;
    this->texture_count = texture_count;
    this->vertex_count = vertex_count;
    this->mass = mass;

    this->textures = new texture[texture_count];
    this->triangles = new triangle[triangle_count];
    this->vector_array = new double[3 * triangle_count + vertex_count][3];
}


Shape::Shape(unsigned int triangle_count, unsigned int texture_count, unsigned int vertex_count) : Shape(triangle_count, texture_count, vertex_count, 1.0) {}

Shape::Shape(const Shape& other) : Shape(other.triangle_count, other.texture_count, other.vertex_count, other.mass)
{
    this->origin[0] = other.origin[0];
    this->origin[1] = other.origin[1];
    this->origin[2] = other.origin[2];

    for (unsigned int i = 0; i < 3; i++) {
        this->velocity[i] = other.velocity[i];
        this->angular_velocity[i] = other.angular_velocity[i];
    }

    for (unsigned int i = 0; i < this->triangle_count; i++) {
        this->triangles[i].origin_vector_index = other.triangles[i].origin_vector_index;
        for (unsigned int j = 0; j < 3; j++) {
            this->triangles[i].side_vector_index = other.triangles[i].side_vector_index;
            this->triangles[i].bitmap = other.triangles[i].bitmap;

            this->triangles[i].side_colors[j][0] = other.triangles[i].side_colors[j][0];
            this->triangles[i].side_colors[j][1] = other.triangles[i].side_colors[j][1];
            this->triangles[i].side_colors[j][2] = other.triangles[i].side_colors[j][2];
        }
    }
    for (unsigned int i = 0; i < this->texture_count; i++) {
        this->textures[i].bitmap = other.textures[i].bitmap;
        this->textures[i].x_size = other.textures[i].x_size;
        this->textures[i].y_size = other.textures[i].y_size;
        this->textures[i].invert_x = other.textures[i].invert_x;
        this->textures[i].invert_y = other.textures[i].invert_y;
        this->textures[i].ref_x_triangle = other.textures[i].ref_x_triangle;
        this->textures[i].ref_x_side = other.textures[i].ref_x_side;
        this->textures[i].ref_x_axis = other.textures[i].ref_x_axis;
        this->textures[i].ref_y_triangle = other.textures[i].ref_y_triangle;
        this->textures[i].ref_y_side = other.textures[i].ref_y_side;
        this->textures[i].ref_y_axis = other.textures[i].ref_y_axis;
        this->textures[i].ref_triangle = other.textures[i].ref_triangle;
        this->textures[i].ref_triangle_index = other.textures[i].ref_triangle_index;

        this->textures[i].origin_vector_index = other.textures[i].origin_vector_index;
    }
    for (unsigned int i = 0; i < 3 * this->triangle_count + this->vertex_count; i++) {
        this->vector_array[i][0] = other.vector_array[i][0];
        this->vector_array[i][1] = other.vector_array[i][1];
        this->vector_array[i][2] = other.vector_array[i][2];
    }
}

Shape& Shape::operator= (const Shape &other)
{
    if (this != &other) {
        this->triangle_count = other.triangle_count;
        this->texture_count = other.texture_count;
        this->vertex_count = other.vertex_count;
        this->mass = other.mass;

        this->triangles = new triangle[triangle_count];
        this->textures = new texture[texture_count];
        // this->vertices = new double[vertex_count][3];
        this->vector_array = new double[3 * triangle_count + vertex_count][3];

        for (unsigned int i = 0; i < 3; i++) {
            this->velocity[i] = other.velocity[i];
            this->angular_velocity[i] = other.angular_velocity[i];
        }

        this->origin[0] = other.origin[0];
        this->origin[1] = other.origin[1];
        this->origin[2] = other.origin[2];

        for (unsigned int i = 0; i < this->triangle_count; i++) {
            this->triangles[i].origin_vector_index = other.triangles[i].origin_vector_index;
            for (unsigned int j = 0; j < 3; j++) {
                this->triangles[i].side_vector_index = other.triangles[i].side_vector_index;
                this->triangles[i].bitmap = other.triangles[i].bitmap;

                this->triangles[i].side_colors[j][0] = other.triangles[i].side_colors[j][0];
                this->triangles[i].side_colors[j][1] = other.triangles[i].side_colors[j][1];
                this->triangles[i].side_colors[j][2] = other.triangles[i].side_colors[j][2];
            }
        }
        for (unsigned int i = 0; i < this->texture_count; i++) {
            this->textures[i].bitmap = other.textures[i].bitmap;
            this->textures[i].x_size = other.textures[i].x_size;
            this->textures[i].y_size = other.textures[i].y_size;
            this->textures[i].invert_x = other.textures[i].invert_x;
            this->textures[i].invert_y = other.textures[i].invert_y;
            this->textures[i].ref_x_triangle = other.textures[i].ref_x_triangle;
            this->textures[i].ref_x_side = other.textures[i].ref_x_side;
            this->textures[i].ref_x_axis = other.textures[i].ref_x_axis;
            this->textures[i].ref_y_triangle = other.textures[i].ref_y_triangle;
            this->textures[i].ref_y_side = other.textures[i].ref_y_side;
            this->textures[i].ref_y_axis = other.textures[i].ref_y_axis;
            this->textures[i].ref_triangle = other.textures[i].ref_triangle;
            this->textures[i].ref_triangle_index = other.textures[i].ref_triangle_index;

            this->textures[i].origin_vector_index = other.textures[i].origin_vector_index;
        }

        for (unsigned int i = 0; i < 3 * this->triangle_count + this->vertex_count; i++) {
            this->vector_array[i][0] = other.vector_array[i][0];
            this->vector_array[i][1] = other.vector_array[i][1];
            this->vector_array[i][2] = other.vector_array[i][2];
        }
    }
    return *this;
}

Shape::~Shape()
{
    delete[] this->triangles;
    delete[] this->textures;
    delete[] this->vector_array;
}

// creates a square using two triangles; hypotenuse goes from upper left to lower right
void create_square(Shape *shape, unsigned int tri_index, unsigned int vector_index, double *bottom_left, double *top_left, double *top_right, double *bottom_right, unsigned int bottom_left_vertex_index, unsigned int top_right_vertex_index) {
    for (unsigned int i = 0; i < 3; i++) {
        shape->vector_array[vector_index][i] = top_left[i] - bottom_left[i];
        shape->vector_array[vector_index + 1][i] = bottom_right[i] - top_left[i];
        shape->vector_array[vector_index + 2][i] = bottom_left[i] - bottom_right[i];

        shape->vector_array[vector_index + 3][i] = top_left[i] - top_right[i];
        shape->vector_array[vector_index + 4][i] = bottom_right[i] - top_left[i];
        shape->vector_array[vector_index + 5][i] = top_right[i] - bottom_right[i];
    }
    shape->triangles[tri_index].origin_vector_index = bottom_left_vertex_index;
    shape->triangles[tri_index + 1].origin_vector_index = top_right_vertex_index;
    shape->triangles[tri_index].side_vector_index = vector_index;
    shape->triangles[tri_index + 1].side_vector_index = vector_index + 3;
}

void create_cube(Shape *shape, double side_len)
{
    // calculate coordinates of vertices
    double shape_origin[] = { 0, 0, 0 };
    double front_bottom_left[] = { shape_origin[0] - side_len / 2.0, shape_origin[1] - side_len / 2.0,
                                    shape_origin[2] - side_len / 2.0 };
    double front_top_left[] = { shape_origin[0] - side_len / 2.0, shape_origin[1] + side_len / 2.0,
                                    shape_origin[2] - side_len / 2.0 };
    double front_top_right[] = { shape_origin[0] + side_len / 2.0, shape_origin[1] + side_len / 2.0,
                                    shape_origin[2] - side_len / 2.0 };
    double front_bottom_right[] = { shape_origin[0] + side_len / 2.0, shape_origin[1] - side_len / 2.0,
                                    shape_origin[2] - side_len / 2.0 };

    double back_bottom_left[] = { shape_origin[0] - side_len / 2.0, shape_origin[1] - side_len / 2.0,
                                    shape_origin[2] + side_len / 2.0 };
    double back_top_left[] = { shape_origin[0] - side_len / 2.0, shape_origin[1] + side_len / 2.0,
                                    shape_origin[2] + side_len / 2.0 };
    double back_top_right[] = { shape_origin[0] + side_len / 2.0, shape_origin[1] + side_len / 2.0,
                                    shape_origin[2] + side_len / 2.0 };
    double back_bottom_right[] = { shape_origin[0] + side_len / 2.0, shape_origin[1] - side_len / 2.0,
                                    shape_origin[2] + side_len / 2.0 };

    // set vertex vectors
    for (unsigned int i = 0; i < 3; i++) {
        shape->vector_array[0][i] = 0;
        shape->vector_array[1][i] = front_bottom_left[i];
        shape->vector_array[2][i] = front_top_left[i];
        shape->vector_array[3][i] = front_top_right[i];
        shape->vector_array[4][i] = front_bottom_right[i];
        shape->vector_array[5][i] = back_bottom_left[i];
        shape->vector_array[6][i] = back_top_left[i];
        shape->vector_array[7][i] = back_top_right[i];
        shape->vector_array[8][i] = back_bottom_right[i];
    }

    int tri_index = 0;
    int vector_index = 9;
    // create front
    create_square(shape, tri_index, vector_index, front_bottom_left, front_top_left,
                    front_top_right, front_bottom_right, 1, 3);
    tri_index += 2;
    vector_index += 6;
    // create back
    create_square(shape, tri_index, vector_index, back_bottom_left, back_top_left,
                    back_top_right, back_bottom_right, 5, 7);
    tri_index += 2;
    vector_index += 6;
    // create left
    create_square(shape, tri_index, vector_index, back_bottom_left, back_top_left,
                    front_top_left, front_bottom_left, 5, 2);
    tri_index += 2;
    vector_index += 6;
    // create right
    create_square(shape, tri_index, vector_index, front_bottom_right, front_top_right,
                    back_top_right, back_bottom_right,  4, 7);
    tri_index += 2;
    vector_index += 6;
    // create bottom
    create_square(shape, tri_index, vector_index, back_bottom_left, front_bottom_left,
                    front_bottom_right, back_bottom_right,  5, 4);
    tri_index += 2;
    vector_index += 6;
    // create top
    create_square(shape, tri_index, vector_index, front_top_left, back_top_left,
            back_top_right, front_top_right,  2, 7);

    // don't draw diagonals
    for (unsigned int i = 0; i < shape->triangle_count; i++) {
        shape->triangles[i].side_colors[1][0] = 0;
        shape->triangles[i].side_colors[1][1] = 0;
        shape->triangles[i].side_colors[1][2] = 0;
    }
}

void create_cube_45(Shape *shape, double side_len)
{
    // calculate coordinates of vertices
    double side_slanted = sqrt(2) / 2.0 * side_len;
    double box_side = side_len - side_slanted;

    double shape_origin[] = { 0, 0, 0 };
    double front_bottom_left[] = { shape_origin[0] - side_slanted - box_side / 2.0,
                                    shape_origin[1] - side_slanted / 2.0, shape_origin[2] - side_len / 2.0  };
    double front_top_left[] = { shape_origin[0] - box_side / 2.0, shape_origin[1] + side_slanted / 2.0,
                                shape_origin[2] - side_len / 2.0 };
    double front_top_right[] = { shape_origin[0] + side_slanted + box_side / 2.0, shape_origin[1] + side_slanted / 2.0,
                                    shape_origin[2] - side_len / 2.0 };
    double front_bottom_right[] = { shape_origin[0] + box_side / 2.0, shape_origin[1] - side_slanted / 2.0,
                                    shape_origin[2] - side_len / 2.0 };

    double back_bottom_left[] = { shape_origin[0] - side_slanted - box_side / 2.0, shape_origin[1] - side_slanted / 2.0,
                                    shape_origin[2] + side_len / 2.0 };
    double back_top_left[] = { shape_origin[0] - box_side / 2.0, shape_origin[1] + side_slanted / 2.0,
                                shape_origin[2] + side_len / 2.0 };
    double back_top_right[] = { shape_origin[0] + side_slanted + box_side / 2.0, shape_origin[1] + side_slanted / 2.0,
                                shape_origin[2] + side_len / 2.0 };
    double back_bottom_right[] = { shape_origin[0] + box_side / 2.0, shape_origin[1] - side_slanted / 2.0,
                                    shape_origin[2] + side_len / 2.0 };

    // set vertex vectors
    for (unsigned int i = 0; i < 3; i++) {
        shape->vector_array[0][i] = 0;
        shape->vector_array[1][i] = front_bottom_left[i];
        shape->vector_array[2][i] = front_top_left[i];
        shape->vector_array[3][i] = front_top_right[i];
        shape->vector_array[4][i] = front_bottom_right[i];
        shape->vector_array[5][i] = back_bottom_left[i];
        shape->vector_array[6][i] = back_top_left[i];
        shape->vector_array[7][i] = back_top_right[i];
        shape->vector_array[8][i] = back_bottom_right[i];
    }

    int tri_index = 0;
    int vector_index = 9;
    // create front
    create_square(shape, tri_index, vector_index, front_bottom_left, front_top_left,
            front_top_right, front_bottom_right, 1, 3);
    tri_index += 2;
    vector_index += 6;
    // create back
    create_square(shape, tri_index, vector_index, back_bottom_left, back_top_left,
            back_top_right, back_bottom_right, 5, 7);
    tri_index += 2;
    vector_index += 6;
    // create left
    create_square(shape, tri_index, vector_index, back_bottom_left, back_top_left,
            front_top_left, front_bottom_left, 5, 2);
    tri_index += 2;
    vector_index += 6;
    // create right
    create_square(shape, tri_index, vector_index, front_bottom_right, front_top_right,
            back_top_right, back_bottom_right,  4, 7);
    tri_index += 2;
    vector_index += 6;
    // create bottom
    create_square(shape, tri_index, vector_index, back_bottom_left, front_bottom_left,
            front_bottom_right, back_bottom_right,  5, 4);
    tri_index += 2;
    vector_index += 6;
    // create top
    create_square(shape, tri_index, vector_index, front_top_left, back_top_left,
            back_top_right, front_top_right,  2, 7);

    // don't draw diagonals
    for (unsigned int i = 0; i < shape->triangle_count; i++) {
        shape->triangles[i].side_colors[1][0] = 0;
        shape->triangles[i].side_colors[1][1] = 0;
        shape->triangles[i].side_colors[1][2] = 0;
    }
}

void assemble_cube(Shape *cube, double *origin_coord, double side_len, double *velocity, double *angular_velocity) {
    create_cube(cube, side_len);

    for (unsigned int i = 0; i < 3; i++) {
        cube->origin[i] = origin_coord[i];
        cube->velocity[i] = velocity[i];
        cube->angular_velocity[i] = angular_velocity[i];
    }

}

void assemble_cube45(Shape * cube, double *origin_coord, double side_len, double *velocity, double *angular_velocity) {
    create_cube_45(cube, side_len);

    for (unsigned int i = 0; i < 3; i++) {
        cube->origin[i] = origin_coord[i];
        cube->velocity[i] = velocity[i];
        cube->angular_velocity[i] = angular_velocity[i];
    }
}
