#include "draw.h"
#include "BitMap.h"
#include "movement.h"
#include "utils.h"

//constexpr double PI = 3.1415926535897932384626433832795028841971693993751058209;
//constexpr double FOV = 150 * PI / 180.0;
//constexpr double MIN_DOUBLE_THRESHOLD = .00001;


void draw_shape(Shape *shape, unsigned int triangle_index, Pixel *bitmap, double **depth_map)
{
    double current_x = shape->origin[0];
    double current_y = shape->origin[1];
    double current_z = shape->origin[2];

    unsigned int triangle_origin_vertex = shape->triangles[triangle_index].origin_vector_index;

    double shape_z = shape->origin[2];
    double z_scaling = IMAGE_WIDTH / (IMAGE_WIDTH + 2 * shape_z * tan(PI / 2.0 - FOV / 2.0));
    current_x = shape->origin[0] + shape->vector_array[triangle_origin_vertex][0] * z_scaling;
    current_y = shape->origin[1] + shape->vector_array[triangle_origin_vertex][1] * z_scaling;
    current_z = shape->origin[2] + shape->vector_array[triangle_origin_vertex][2] * z_scaling;
    double perception_x_side_change = shape_z / tan(FOV / 2.0);
    double perception_x_range = IMAGE_WIDTH + 2 * perception_x_side_change;
    double corrected_x_origin = (shape->origin[0] + perception_x_side_change) / perception_x_range * IMAGE_WIDTH;
    current_x = corrected_x_origin + shape->vector_array[triangle_origin_vertex][0] * z_scaling;

    int y_image_size = IMAGE_HEIGHT;
    double perception_y_side_change = shape_z / tan(FOV / 2.0) * 9.0 / 16.0;
    double perception_y_range = y_image_size + 2 * perception_y_side_change;
    double corrected_y_origin = (shape->origin[1] + perception_y_side_change) / perception_y_range * IMAGE_HEIGHT;
    current_y = corrected_y_origin + shape->vector_array[triangle_origin_vertex][1] * z_scaling;

    // loop over each side vector
    for (unsigned int side_index = 0; side_index < 3; side_index++) {
        unsigned int vector_index = shape->triangles[triangle_index].side_vector_index + side_index;
        double vector_x = shape->vector_array[vector_index][0];
        double vector_y = shape->vector_array[vector_index][1];
        double vector_z = shape->vector_array[vector_index][2];

        // normalize vector increments to max size .5; prevents skipping pixels in bitmap
        double increment_sum = abs(vector_x) + abs(vector_y) + abs(vector_z);
        double x_increment = (increment_sum == 0) ? 0 : vector_x / increment_sum / 2.0;
        double y_increment = (increment_sum == 0) ? 0 : vector_y / increment_sum / 2.0;
        double z_increment = (increment_sum == 0) ? 0 : vector_z / increment_sum / 2.0;

        double dest_x = current_x + vector_x * z_scaling;
        double dest_y = current_y + vector_y * z_scaling;
        double dest_z = current_z + vector_z * z_scaling;
        bool x_reached = false, y_reached = false, z_reached = false;
        bool no_color = shape->triangles[triangle_index].side_colors[side_index][0] == 0 &&
                        shape->triangles[triangle_index].side_colors[side_index][1] == 0 &&
                        shape->triangles[triangle_index].side_colors[side_index][2] == 0;
        while (!x_reached || !y_reached || !z_reached) {
            const int x_pixel = (int) current_x;
            const int y_pixel = (int) current_y;

            // do not color if index out of bounds or color is black (same as background; can cover drawn lines)
            if (!no_color && current_x >= 0 && current_x <= IMAGE_WIDTH && current_y >= 0 &&
                current_y <= IMAGE_HEIGHT && current_z > -1) {

                const unsigned int pixel_index = y_pixel * IMAGE_WIDTH + x_pixel;
                if (depth_map[y_pixel][x_pixel] < 0 || current_z < depth_map[y_pixel][x_pixel]) {
                    unsigned int bitmap_index = y_pixel * IMAGE_WIDTH + x_pixel;
                    bitmap[bitmap_index].blue = shape->triangles[triangle_index].side_colors[side_index][0];
                    bitmap[bitmap_index].green = shape->triangles[triangle_index].side_colors[side_index][1];
                    bitmap[bitmap_index].red = shape->triangles[triangle_index].side_colors[side_index][2];

                    depth_map[y_pixel][x_pixel] = current_z;
                }
            }

            // update progress booleans
            x_reached = abs(x_increment) < 0.01 ||
                        (x_increment > 0 ? trunc(current_x) >= trunc(dest_x) : trunc(current_x) <= trunc(dest_x));
            y_reached = abs(y_increment) < 0.01 ||
                        (y_increment > 0 ? trunc(current_y) >= trunc(dest_y) : trunc(current_y) <= trunc(dest_y));
            z_reached = abs(z_increment) < .01 ||
                        (z_increment > 0 ? trunc(current_z) >= trunc(dest_z) : trunc(current_z) <= trunc(dest_z));

            // increment coords until new pixel is reached or z-axis destination is reached
            do {
                current_x += x_increment;
                current_y += y_increment;
                current_z += z_increment;
                z_reached = abs(z_increment) < .01 ||
                            (z_increment > 0 ? trunc(current_z) >= trunc(dest_z) : trunc(current_z) <= trunc(dest_z));
            } while (trunc(current_x) == x_pixel && trunc(current_y) == y_pixel && !z_reached);
        }
    }
}

void draw_texture(Shape *shape, unsigned int tri_index, Pixel *bitmap, double **depth_map)
{
    // find triangle's origin
    double shape_z = shape->origin[2];
    double z_scaling = IMAGE_WIDTH / (IMAGE_WIDTH + 2 * shape_z * tan(PI / 2.0 - FOV / 2.0));

    double calc_shape_origin_coords[3];
    unsigned int origin_vertex_index = shape->triangles[tri_index].origin_vector_index;
    calc_shape_origin_coords[0] = shape->origin[0] + shape->vector_array[origin_vertex_index][0] * z_scaling;
    calc_shape_origin_coords[1] = shape->origin[1] + shape->vector_array[origin_vertex_index][1] * z_scaling;
    calc_shape_origin_coords[2] = shape->origin[2] + shape->vector_array[origin_vertex_index][2] * z_scaling;

    // map shape coord origin to bitmap coord origin
    double perception_x_side_change = shape_z / tan(FOV / 2.0);
    double perception_x_range = IMAGE_WIDTH + 2 * perception_x_side_change;
    double corrected_x_origin = (shape->origin[0] + perception_x_side_change) / perception_x_range * IMAGE_WIDTH;
    calc_shape_origin_coords[0] = corrected_x_origin + shape->vector_array[origin_vertex_index][0] * z_scaling;

    int y_image_size = IMAGE_HEIGHT;
    double perception_y_side_change = shape_z / tan(FOV / 2.0) * 9.0 / 16.0;
    double perception_y_range = y_image_size + 2 * perception_y_side_change;
    double corrected_y_origin = (shape->origin[1] + perception_y_side_change) / perception_y_range * IMAGE_HEIGHT;
    calc_shape_origin_coords[1] = corrected_y_origin + shape->vector_array[origin_vertex_index][1] * z_scaling;

    // save scaled triangle sides in local variable
    double side_1[3];
    double side_2[3];
    double side_3[3];
    unsigned int side_vertex_index = shape->triangles[tri_index].side_vector_index;
    for (unsigned int i = 0; i < 3; i++) {
        side_1[i] = shape->vector_array[side_vertex_index][i] * z_scaling;
        side_2[i] = shape->vector_array[side_vertex_index + 1][i] * z_scaling;
        side_3[i] = shape->vector_array[side_vertex_index + 2][i] * z_scaling;

    }

    // return if triangle's x/y area is 0
    if ((abs(side_1[0]) <= MIN_DOUBLE_THRESHOLD && abs(side_1[1]) <= MIN_DOUBLE_THRESHOLD) ||
        (abs(side_2[0]) <= MIN_DOUBLE_THRESHOLD && abs(side_2[1]) <= MIN_DOUBLE_THRESHOLD) ||
        (abs(side_3[0]) <= MIN_DOUBLE_THRESHOLD && abs(side_3[1]) <= MIN_DOUBLE_THRESHOLD)) {
        return;
    }

    // find bounding box for triangle
    double coord_x_limit[2];
    double coord_y_limit[2];
    double coord_z_limit[2];
    triangle_min_max_coords(side_1, side_2, side_3, coord_x_limit, coord_y_limit, coord_z_limit,
                            calc_shape_origin_coords, 1);

    // z-coord must be ignored for triangle bounds check but needed to determine closer shape if overlap exists
    double true_x_z = side_1[2];
    double true_y_z = side_2[2];
    double true_z_z = side_3[2];
    side_1[2] = 0;
    side_2[2] = 0;
    side_3[2] = 0;

    // find normal vector of triangle's plane
    double normal_vector[3];
    cross(normal_vector, side_1[0], side_1[1], true_x_z, side_2[0], side_2[1], true_y_z);

    // loop over bounding box
    for (int bitmap_y_index = coord_y_limit[0]; bitmap_y_index < coord_y_limit[1]; bitmap_y_index++) {
        for (int bitmap_x_index = coord_x_limit[0]; bitmap_x_index < coord_x_limit[1]; bitmap_x_index++) {
            if (plane_point_in_triangle(side_1, side_2, side_3, bitmap_x_index, bitmap_y_index,
                calc_shape_origin_coords[2], calc_shape_origin_coords)) {

                // find z coordinate
                double z_coord_bitmap = (-1 * normal_vector[0] * (bitmap_x_index - calc_shape_origin_coords[0])
                                        - normal_vector[1] * (bitmap_y_index - calc_shape_origin_coords[1]))
                                        / (normal_vector[2]);
                z_coord_bitmap += calc_shape_origin_coords[2];

                int texture_bitmap_index = 0; // placeholder pending mapping adjustments for bitmaps

                bool valid_coord = bitmap_x_index >= 0 && bitmap_y_index >= 0 &&
                                    bitmap_x_index < IMAGE_WIDTH && bitmap_y_index < IMAGE_HEIGHT;
                if (valid_coord && (depth_map[bitmap_y_index][bitmap_x_index] < 0 ||
                                    z_coord_bitmap < depth_map[bitmap_y_index][bitmap_x_index])) {

                    // update bitmap pixel
                    unsigned int bitmap_index = bitmap_y_index * IMAGE_WIDTH + bitmap_x_index;
                    bitmap[bitmap_index].blue = shape->triangles[tri_index].bitmap[texture_bitmap_index].blue;
                    bitmap[bitmap_index].green = shape->triangles[tri_index].bitmap[texture_bitmap_index].green;
                    bitmap[bitmap_index].red = shape->triangles[tri_index].bitmap[texture_bitmap_index].red;

                    // update depth map
                    depth_map[bitmap_y_index][bitmap_x_index] = z_coord_bitmap;
                }
            }
        }
    }
}
