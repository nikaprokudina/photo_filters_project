// #define _CRT_SECURE_NO_WARNINGS
// #define _USE_MATH_DEFINES

#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "util.h"

// typedef unsigned long DWORD;
// typedef unsigned short WORD;
// typedef long LONG;

using WORD = uint16_t;
using DWORD = uint16_t;
using LONG = uint16_t;

#pragma pack(push)
#pragma pack(1)
struct BitmapFileHeader {
    using WORD = uint16_t;
    using DWORD = unsigned int;

    WORD bfType;
    DWORD bfSize;
    WORD bfReserved1;
    WORD bfReserved2;
    DWORD bfOffBits;
};

struct BitmapInfoHeader {
    using DWORD = unsigned int;
    using LONG = int;
    using WORD = uint16_t;

    DWORD biSize;
    LONG biWidth;
    LONG biHeight;
    WORD biPlanes;
    WORD biBitCount;
    DWORD biCompression;
    DWORD biSizeImage;
    LONG biXPelsPerMeter;
    LONG biYPelsPerMeter;
    DWORD biClrUsed;
    DWORD biClrImportant;
};

#pragma pack(pop)

class BMPFile {
    unsigned char *iz_rgb_;
    BitmapFileHeader file_header_;
    BitmapInfoHeader info_header_;
    uint32_t bytes_per_pixel_;
    uint32_t padding_bytes_;
    uint32_t stride_;
    const int count_bytes_ = 8;
    const unsigned int dvoichn_ = 0xFFFFFFFC;

public:
    uint32_t Height() const;

    uint32_t Width() const;

    uint32_t BitCount() const;

    void SetNewIzRGB(unsigned char *data);

    void SetDimensions(uint32_t new_w, uint32_t new_h);

    BMPFile(BMPFile &other);

    explicit BMPFile(char *input_file_name);

    ~BMPFile();

    void Export(char *out_file_name);

    BMPFile *Grayscale();

    BMPFile *Negative();

    std::unique_ptr<BMPFile> Sharpening();

    BMPFile *Crop(uint32_t width, uint32_t height);

    // BMPFile *GaussianBlur(double sigma);

    // BMPFile *Sharpening();

    BMPFile *EdgeDetection(const double threshold);

    BMPFile *Space(double x0, double y0, double speed, double extent);
};