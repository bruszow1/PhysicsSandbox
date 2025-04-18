#include <cstdint> // for specific size integers
#include <string>
#include <fstream> // for file handling
#include "bitmap.h"

#include <iostream>

using std::string;

struct BmpHeader { // indicates file type
    char bitmapSignatureBytes[2] = {'B', 'M'}; // signature to indicate file is a bitmap
    uint32_t sizeOfBitmapFile = 54 + (IMAGE_WIDTH * IMAGE_HEIGHT) * 3; // metadata + (pixel_width * pixel_height) * color_bytes
    uint32_t reservedBytes = 0;
    uint32_t pixelDataOffset = 54;
} bmp_header;

struct BmpInfoHeader { // additional image info
    uint32_t sizeOfThisHeader = 40;
    int32_t  width = IMAGE_WIDTH; // in pixels
    int32_t  height = IMAGE_HEIGHT; // in pixels; if negative, first pixel is top left : otherwise, starts bottom left
    uint16_t numberOfColorPlanes = 1; // must be 1
    uint16_t colorDepth = 24; // bits per pixel
    uint32_t compressionMethod = 0;
    uint32_t rawBitmapDataSize = 0; // generally ignored
    int32_t  horizontalResolution = 3780; // in pixel per meter
    int32_t  verticalResolution = 3780; // in pixel per meter
    uint32_t colorTableEntries = 0;
    uint32_t importantColors = 0;
} bmp_info_header;



int write_bitmap_to_file(Pixel *test_bitmap, const string &file_name) {
    // std::string file_name = "test.bmp";
    // constexpr unsigned int pixel_count = IMAGE_HEIGHT * IMAGE_WIDTH;
    // Pixel *test_bitmap = new Pixel[pixel_count];
    // for (int i = 0; i < pixel_count; i++) {
    //     Pixel pixel;
    //     test_bitmap[i] = pixel;
    //     test_bitmap[i].red = 110;
    //     test_bitmap[i].green = 110;
    //     test_bitmap[i].blue = 110;
    // }

    std::ofstream file_out(file_name, std::ofstream::binary);

    // file_out.write((char *) (&bmp_header), 14);

    file_out.write((char *) (&bmp_header.bitmapSignatureBytes), 2);
    file_out.write((char *) (&bmp_header.sizeOfBitmapFile), sizeof(uint32_t));
    file_out.write((char *) (&bmp_header.reservedBytes), sizeof(uint32_t));
    file_out.write((char *) (&bmp_header.pixelDataOffset), sizeof(uint32_t));

    // file_out.write((char *) (&bmp_info_header), 40);

    file_out.write((char *) (&bmp_info_header.sizeOfThisHeader), sizeof(uint32_t));
    file_out.write((char *) (&bmp_info_header.width), sizeof(int32_t));
    file_out.write((char *) (&bmp_info_header.height), sizeof(int32_t));
    file_out.write((char *) (&bmp_info_header.numberOfColorPlanes), sizeof(uint16_t));
    file_out.write((char *) (&bmp_info_header.colorDepth), sizeof(uint16_t));
    file_out.write((char *) (&bmp_info_header.compressionMethod), sizeof(uint32_t));
    file_out.write((char *) (&bmp_info_header.rawBitmapDataSize), sizeof(uint32_t));
    file_out.write((char *) (&bmp_info_header.horizontalResolution), sizeof(int32_t));
    file_out.write((char *) (&bmp_info_header.verticalResolution), sizeof(int32_t));
    file_out.write((char *) (&bmp_info_header.colorTableEntries), sizeof(uint32_t));
    file_out.write((char *) (&bmp_info_header.importantColors), sizeof(uint32_t));


    const size_t pixel_count = IMAGE_WIDTH * IMAGE_HEIGHT;
    for (int i = 0; i < pixel_count; i++) {
        // file_out.write((char *) (&test_bitmap[i]), 3);
        file_out.write(((char *) &test_bitmap[i].blue), sizeof(uint8_t));
        file_out.write(((char *) &test_bitmap[i].green), sizeof(uint8_t));
        file_out.write(((char *) &test_bitmap[i].red), sizeof(uint8_t));
    }
    file_out.close();

    return 0;
}
