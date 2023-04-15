#include "image.h"
// #include <string.h>

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "неправильный ввод"
                  << "\n";
        return -1;
    }
    BMPFile *image = new BMPFile(argv[1]);
    BMPFile *processed_image = nullptr;
    uint32_t number = 3;
    double x0 = 0;
    double y0 = 0;
    double speed = 0;
    double extent = 0;

    while (number < argc) {
        if (StrCmp(argv[number], "-gs") == 0) {
            processed_image = image->Grayscale();
            delete image;
            image = processed_image;
        } else if (StrCmp(argv[number], "-neg") == 0) {
            processed_image = image->Negative();
            delete image;
            image = processed_image;
        } else if (StrCmp(argv[number], "-sharp") == 0) {
            processed_image = image->Sharpening();
            delete image;
            image = processed_image;
        } else if (StrCmp(argv[number], "-space") == 0) {
            if (number + 4 >= argc) {
                std::cerr << "Error: missing arguments for option -space\n";
                break;
            }
            x0 = std::strtod(argv[number + 1], nullptr);
            y0 = std::strtod(argv[number + 2], nullptr);
            speed = std::strtod(argv[number + 3], nullptr);
            extent = std::strtod(argv[number + 4], nullptr);
            processed_image = image->Space(x0, y0, speed, extent);
            delete image;
            image = processed_image;
            number += 4;
        } else if (StrCmp(argv[number], "-edge") == 0) {
            double t = 0;
            if (number + 1 >= argc) {
                std::cerr << "Error: missing argument for option -edge\n";
                break;
            }
            t = std::strtod(argv[number + 1], nullptr);
            processed_image = image->Grayscale()->EdgeDetection(t);
            delete image;
            image = processed_image;
            ++number;
        } else {
            std::cerr << "Error: unknown option " << argv[number] << "\n";
            break;
        }
        ++number;
    }

    if (image != nullptr) {
        image->Export(argv[2]);
        uint32_t w = image->Width();
        uint32_t h = image->Height();
        std::cout << "Output image dimensions: " << w << " x " << h << " pixels, " << image->BitCount()
                  << " bits per pixel\n";
        delete image;
    }

    return 0;
}
