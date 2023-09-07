#include "Image.h"

Image::Image()
    : width(0), height(0)
{
}

Image::Image(std::int32_t width, std::int32_t height)
    : width(width), height(height)
{
    data.resize(width * height);
    bmpHeader.sizeOfBitmapFile = width * height * 3 + 54;
    bmpInfoHeader.width = width;
    bmpInfoHeader.height = height;
    bmpInfoHeader.rawBitmapDataSize = 0;
    bmpInfoHeader.colorDepth = 32;
    bmpInfoHeader.numberOfColorPlanes = 1;
    bmpInfoHeader.importantColors = 0;
    bmpInfoHeader.colorTableEntries = 0;
}

Image::Image(std::string filename) 
    : width(0), height(0)
{
    std::ifstream ifstream(filename);
    nlohmann::json x = nlohmann::json::parse(ifstream);

    for (int i = 0; i < x["Palette"].size(); i++)
    {
        Image::pixel_t p;
        p.b = x["Palette"][i][0];
        p.g = x["Palette"][i][1];
        p.r = x["Palette"][i][2];
        p.a = x["Palette"][i][3];

        this->data.push_back(p);
    }
}

Image::~Image()
{
    data.swap(data);
}

bool Image::operator==(const Image* other)
{
    return std::equal(data.begin(), data.end(), other->data.begin());
}

std::ostream& operator<<(std::ostream& output, const Image* image)
{
    output.write((char*)&image->bmpHeader.bitmapSignatureBytes, sizeof(char) * 2);
    output.write((char*)&image->bmpHeader.sizeOfBitmapFile, sizeof(uint32_t));
    output.write((char*)&image->bmpHeader.reservedBytes, sizeof(uint32_t));
    output.write((char*)&image->bmpHeader.pixelDataOffset, sizeof(uint32_t));

    output.write((char*)&image->bmpInfoHeader.sizeOfThisHeader, sizeof(uint32_t));
    output.write((char*)&image->bmpInfoHeader.width, sizeof(int32_t));
    output.write((char*)&image->bmpInfoHeader.height, sizeof(int32_t));
    output.write((char*)&image->bmpInfoHeader.numberOfColorPlanes, sizeof(uint16_t));
    output.write((char*)&image->bmpInfoHeader.colorDepth, sizeof(uint16_t));
    output.write((char*)&image->bmpInfoHeader.compressionMethod, sizeof(uint32_t));
    output.write((char*)&image->bmpInfoHeader.rawBitmapDataSize, sizeof(uint32_t));
    output.write((char*)&image->bmpInfoHeader.horizontalResolution, sizeof(int32_t));
    output.write((char*)&image->bmpInfoHeader.verticalResolution, sizeof(int32_t));
    output.write((char*)&image->bmpInfoHeader.colorTableEntries, sizeof(uint32_t));
    output.write((char*)&image->bmpInfoHeader.importantColors, sizeof(uint32_t));

    std::size_t row_size = image->width * sizeof(Image::pixel_t);

    for (int i = image->height - 1; i >= 0; i--)
    {
        output.write((char*)&image->data[i * image->width], row_size);
    }

    return output;
}

std::istream& operator>>(std::istream& input, Image* image)
{
    std::string format;
    input >> format;

    std::string width, height;
    input >> width >> height;

    std::string tmp;
    input >> tmp;

    image->width = std::stoi(width);
    image->height = std::stoi(height);

    image->data.resize(image->width * image->height);

    input.read((char*)&image->data[0], image->data.size() * sizeof(Image::pixel_t));

    return input;
}

bool Image::pixel_t::operator==(const pixel_t& other)
{
    return r == other.r && g == other.g && b == other.b;
}