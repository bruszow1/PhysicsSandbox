#include <iostream>
#include "BitMap.h"

int main() {
    std::cout << "Hello, World!" << std::endl;

    constexpr unsigned int pixel_count = IMAGE_HEIGHT * IMAGE_WIDTH;
    Pixel *test_bitmap = new Pixel[pixel_count];
    for (int i = 0; i < pixel_count; i++) {
        Pixel pixel;
        test_bitmap[i] = pixel;
        test_bitmap[i].red = 110;
        test_bitmap[i].green = 110;
        test_bitmap[i].blue = 110;
    }


    std::string file_name = "test.png";
    BITMAP_H::write_bitmap_to_file();
    delete[] test_bitmap;
    return 0;
}