#pragma once

#include <vector>
#include <cstdint>
#include <sstream>
#include <algorithm>
#include <functional>
#include <iostream>

#include <fstream>
#include <nlohmann/json.hpp>

class Image
{
public:

    struct pixel_t
    {
        std::uint8_t r, g, b, a;

        bool operator==(const pixel_t& other);
    };

    struct BmpHeader {
        char bitmapSignatureBytes[2] = { 'B', 'M' };
        uint32_t sizeOfBitmapFile = 54; // width * height * channels + 14 + 40
        uint32_t reservedBytes = 0;
        uint32_t pixelDataOffset = 54;
    } bmpHeader;

    struct BmpInfoHeader {
        uint32_t sizeOfThisHeader = 40;
        int32_t width = 512;
        int32_t height = 512;
        uint16_t numberOfColorPlanes = 1;
        uint16_t colorDepth = 24;
        uint32_t compressionMethod = 0;
        uint32_t rawBitmapDataSize = 512 * 512 * 3;
        int32_t horizontalResolution = 0;
        int32_t verticalResolution = 0;
        uint32_t colorTableEntries = 0;
        uint32_t importantColors = 0;
    } bmpInfoHeader;

    Image();

    Image(std::int32_t width, std::int32_t height);

    Image(std::string filename);

    ~Image();

    std::vector<pixel_t> data;

    std::int32_t width;

    std::int32_t height;

    bool operator==(const Image* other);

    friend std::ostream& operator<<(std::ostream& output, const Image* other);

    friend std::istream& operator>>(std::istream& input, Image* other);

};