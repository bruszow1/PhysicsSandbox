#ifndef BITMAP_H
#define BITMAP_H

#include <cstdint> // for specific size integers
#include <string>

constexpr unsigned int IMAGE_WIDTH = 1280;
constexpr unsigned int IMAGE_HEIGHT = 720;

struct Pixel {
    uint8_t blue = 0;
    uint8_t green = 0;
    uint8_t red = 0;
};

int write_bitmap_to_file(Pixel *test_bitmap, const std::string &file_name);

#endif //BITMAP_H
