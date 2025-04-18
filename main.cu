#include "bitmap.h"
#include "movement.h"
#include "draw.h"


int main(int argc, char *argv[]) {
    int iteration_count = 50;

    Pixel *bitmap = new Pixel[IMAGE_WIDTH * IMAGE_HEIGHT];

    // create shapes
    Shape **shape_arr;
    unsigned int shape_count = 2;
    shape_arr = new Shape*[shape_count];

    shape_arr[0] = new Shape(12, 12, 9);
    double origin1[] = {200, 175, 425};
    double velocity1[] = {5, 0, 0};
    double angular_velocity1[] = {0, 0, 0};
    assemble_cube45(shape_arr[0], origin1, 300, velocity1, angular_velocity1);

    shape_arr[1] = new Shape(12, 12, 9);
    double origin2[] = {650, 300, 500};
    double velocity2[] = {0, 0, 0};
    double angular_velocity2[] = {0, 0, 0};
    assemble_cube(shape_arr[1], origin2, 300, velocity2, angular_velocity2);

    draw_shapes(bitmap, shape_arr, shape_count);
    const std::string file_name_start = "initial_state.bmp";
    write_bitmap_to_file(bitmap, file_name_start);

    for (unsigned int i = 0; i < iteration_count; i++) {
        model_movement(shape_arr, shape_count, 0, 0, 0, 0, 1);
    }

    draw_shapes(bitmap, shape_arr, shape_count);
    const std::string file_name_end = "final_state.bmp";
    write_bitmap_to_file(bitmap, file_name_end);

    delete[] shape_arr;
    return 0;
}