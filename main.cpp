//#define _CRT_SECURE_NO_WARNINGS
//#define _USE_MATH_DEFINES

#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>

//typedef unsigned long DWORD;
//typedef unsigned short WORD;
//typedef long LONG;

typedef uint16_t WORD;
typedef uint32_t DWORD;
typedef int32_t LONG;

#pragma pack(push)
#pragma pack(1)
typedef struct tagBITMAPFILEHEADER {
    WORD bfType;
    DWORD bfSize;
    WORD bfReserved1;
    WORD bfReserved2;
    DWORD bfOffBits;
} BITMAPFILEHEADER, *LPBITMAPFILEHEADER, *PBITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
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
} BITMAPINFOHEADER, *LPBITMAPINFOHEADER, *PBITMAPINFOHEADER;
#pragma pack(pop)

class Point {
public:
    double x;
    double y;
    double ro;
    double phi;

    void MakePolar() {
        ro = sqrt(x * x + y * y);
        if (x != 0) {
            phi = atan(y / x);
        } else {
            phi = 0;
        }
    }

    void MakeGrizzly() {
        x = ro * cos(phi);
        y = ro * sin(phi);
    }

    Point(double a, double b, bool ispolar = false) {
        if (ispolar) {
            ro = a;
            phi = b;
            MakeGrizzly();
        } else {
            x = a;
            y = b;
            MakePolar();
        }
    }
};


class BMPFile {
    unsigned char *Iz_RGB;
    BITMAPFILEHEADER FileHeader;
    BITMAPINFOHEADER InfoHeader;
    uint32_t bytes_per_pixel;
    uint32_t padding_bytes;
    uint32_t stride;
public:
    uint32_t Height() {
        return InfoHeader.biHeight;
    }

    uint32_t Width() {
        return InfoHeader.biWidth;
    }

    uint32_t BitCount() {
        return InfoHeader.biBitCount;
    }

    void SetNewIzRGB(unsigned char *data) {
        Iz_RGB = data;
    }

    void SetDimensions(uint32_t new_w, uint32_t new_h) {
        InfoHeader.biWidth = new_w;
        InfoHeader.biHeight = new_h;
        bytes_per_pixel = BitCount() / 8;
        stride = new_w * bytes_per_pixel;
        if ((stride & 3) != 0) stride = (stride & 0xFFFFFFFC) + 4;
        uint32_t new_biSizeImage = stride * new_h;
        FileHeader.bfSize += (new_biSizeImage - InfoHeader.biSizeImage);
        InfoHeader.biSizeImage = new_biSizeImage;


    }

    BMPFile(BMPFile &other) {
        memcpy(&FileHeader, &(other.FileHeader), sizeof(BITMAPFILEHEADER));
        memcpy(&InfoHeader, &(other.InfoHeader), sizeof(BITMAPINFOHEADER));
        Iz_RGB = NULL;
    }

    BMPFile(char *input_file_name) {
        FILE *stream;
        stream = fopen(input_file_name, "rb");
        fread(&FileHeader, sizeof(FileHeader), 1, stream);
        fread(&InfoHeader, sizeof(InfoHeader), 1, stream);
        stride = InfoHeader.biWidth * (InfoHeader.biBitCount / 8);
        if ((stride & 3) != 0) stride = (stride & 0xFFFFFFFC) + 4;
        if (InfoHeader.biSizeImage == 0) {//### 7e315
            InfoHeader.biSizeImage = InfoHeader.biHeight * stride;
        }
        fseek(stream, FileHeader.bfOffBits, SEEK_SET); //указатель позиции в файле (сейчас стоит в начале)
        Iz_RGB = (unsigned char *) malloc(InfoHeader.biSizeImage); //выделяет блок памяти под размер файла
        fread(Iz_RGB, InfoHeader.biSizeImage, 1, stream);
        fclose(stream);
        bytes_per_pixel = BitCount() / 8;
        padding_bytes = (4 - (Width() * bytes_per_pixel) % 4) % 4;
    }

    ~BMPFile() { //деструктор
        free(Iz_RGB);
    }

    void Export(char *out_file_name) {
        FILE *outfile = fopen(out_file_name, "wb");
        fwrite(&FileHeader, sizeof(char), sizeof(BITMAPFILEHEADER), outfile);
        fwrite(&InfoHeader, sizeof(char), sizeof(BITMAPINFOHEADER), outfile);
        fwrite(Iz_RGB, sizeof(unsigned char), InfoHeader.biSizeImage, outfile);
        fclose(outfile);
    }

    BMPFile *Grayscale() {
        auto *new_Iz_RGB = (unsigned char *) malloc(InfoHeader.biSizeImage);
        uint32_t i;
        for (int y = 0; y < Height(); ++y) {
            for (int x = 0; x < Width(); ++x) {
                i = y * stride + x * bytes_per_pixel;
                unsigned char b = Iz_RGB[i];
                unsigned char g = Iz_RGB[i + 1];
                unsigned char r = Iz_RGB[i + 2];
                auto gray = (unsigned char) (0.299 * r + 0.587 * g + 0.114 * b);
                new_Iz_RGB[i] = gray; //b
                new_Iz_RGB[i + 1] = gray; //g
                new_Iz_RGB[i + 2] = gray; //r
            }
        }
        auto *f = new BMPFile(*this);
        f->SetNewIzRGB(new_Iz_RGB);
        return f;
    }

    BMPFile *Negative() {
        auto *new_Iz_RGB = (unsigned char *) malloc(InfoHeader.biSizeImage);
        uint32_t i;
        for (int y = 0; y < Height(); ++y) {
            for (int x = 0; x < Width(); ++x) {
                i = y * stride + x * bytes_per_pixel;
                unsigned char b = Iz_RGB[i];
                unsigned char g = Iz_RGB[i + 1];
                unsigned char r = Iz_RGB[i + 2];
                auto gray = (unsigned char) (0.299 * r + 0.587 * g + 0.114 * b);
                new_Iz_RGB[i] = (unsigned char) (255 - b); //b
                new_Iz_RGB[i + 1] = (unsigned char) (255 - g); //g
                new_Iz_RGB[i + 2] = (unsigned char) (255 - r); //r
            }
        }
        auto *f = new BMPFile(*this);
        f->SetNewIzRGB(new_Iz_RGB);
        return f;
    }

    BMPFile *Crop(uint32_t width, uint32_t height) {
        if (width > InfoHeader.biWidth) {
            width = InfoHeader.biWidth;
        }
        if (height > InfoHeader.biHeight) {
            height = InfoHeader.biHeight;
        }
        uint32_t newStride = width * bytes_per_pixel;
        if ((newStride & 3) != 0) newStride = (newStride & 0xFFFFFFFC) + 4;
        uint32_t new_biSizeImage = newStride * height;
        auto *newIz_RGB = (unsigned char *) malloc(new_biSizeImage);
        uint32_t i;
        uint32_t j;
        uint32_t y2 = 0;
        for (uint32_t y = InfoHeader.biHeight - height; y < InfoHeader.biHeight; ++y) {
            for (uint32_t x = 0; x < width; ++x) {
                i = y * stride + x * bytes_per_pixel;
                j = y2 * newStride + x * bytes_per_pixel;
                for (uint32_t k = 0; k < bytes_per_pixel; ++k) {
                    newIz_RGB[j + k] = Iz_RGB[i + k];
                }

            }
            ++y2;
        }
        auto *f = new BMPFile(*this);
        f->SetNewIzRGB(newIz_RGB);
        f->SetDimensions(width, height);
        return f;
    }

    BMPFile *Gaussian_Blur(double sigma) {
        auto *new_Iz_RGB = (unsigned char *) malloc(InfoHeader.biSizeImage);
        auto *new_Iz_RGB_2 = (double *) malloc(InfoHeader.biSizeImage * sizeof(double));
        uint32_t i;
        double xx;
        double yy;
        double t = 1 / (2 * sigma * sigma);
        double k = t / M_PI;
        double mq;
        double porog = 6 - log(k);
        double p;
        int d = 3;
        int x_min;
        int x_max;
        int y_min;
        int y_max;
        for (int y = 0; y < Height(); ++y) {
            for (int x = 0; x < Width(); ++x) {
                i = y * stride + x * bytes_per_pixel;
                new_Iz_RGB_2[i] = 0;
                new_Iz_RGB_2[i + 1] = 0;
                new_Iz_RGB_2[i + 2] = 0;
                x_min = x - 3;
                if (x_min < 0) {
                    x_min = 0;
                }
                y_min = y - 3;
                if (y_min < 0) {
                    y_min = 0;
                }
                x_max = x + 3;
                if (x_max >= Width()) {
                    x_max = Width() - 1;
                }
                y_max = y + 3;
                if (y_max >= Height()) {
                    y_max = Height() - 1;
                }
                for (int y1 = y_min; y1 <= y_max; ++y1) {
                    yy = y1 - y;
                    yy *= yy;
                    for (int x1 = x_min; x1 <= x_max; ++x1) {
                        xx = x1 - x;
                        xx *= xx;
                        p = (xx + yy) * t;
                        if (p > porog) {
                            continue;
                        }
                        mq = k * exp(-p);
                        uint32_t ii = (y1 * (Width() * bytes_per_pixel + padding_bytes)) + (x1 * bytes_per_pixel);
                        unsigned char b = Iz_RGB[ii];
                        unsigned char g = Iz_RGB[ii + 1];
                        unsigned char r = Iz_RGB[ii + 2];
                        new_Iz_RGB_2[i] += (b * mq); //b
                        new_Iz_RGB_2[i + 1] += (g * mq); //g
                        new_Iz_RGB_2[i + 2] += (r * mq); //r
                    }
                }
                new_Iz_RGB[i] = (unsigned char) std::round(new_Iz_RGB_2[i]);
                new_Iz_RGB[i + 1] = (unsigned char) std::round(new_Iz_RGB_2[i + 1]);
                new_Iz_RGB[i + 2] = (unsigned char) std::round(new_Iz_RGB_2[i + 2]);
            }
        };
        free(new_Iz_RGB_2);
        auto *f = new BMPFile(*this);
        f->SetNewIzRGB(new_Iz_RGB);
        return f;
    }

    void
    Matrix_Multiplication(const unsigned char *Iz, int32_t x, int32_t y, const double *matrix, int32_t w, double *Iz2,
                          int32_t size,
                          int32_t stride) {
        for (int32_t i = 0; i < w; ++i) {
            for (int32_t j = 0; j < w; ++j) {
                for (int32_t k = 0; k < w; ++k) {
                    if (((x + i) * 3 >= stride) || ((x + k) * 3 >= stride)) continue;
                    int32_t p1 = (y + j) * stride + (x + k) * 3;
                    int32_t p2 = (y + k) * stride + (x + i) * 3;
                    if ((p1 + 2 > size) || (p2 + 2 > size)) continue;
                    Iz2[p1 + 0] += Iz[p2 + 0] * matrix[k * w + i];
                    Iz2[p1 + 1] += Iz[p2 + 1] * matrix[k * w + i];
                    Iz2[p1 + 2] += Iz[p2 + 2] * matrix[k * w + i];
                }
            }
        }
    }

    BMPFile *Sharpening() {
        auto *new_Iz_RGB = (unsigned char *) malloc(InfoHeader.biSizeImage);
        auto *new_Iz_RGB_2 = (double *) malloc(InfoHeader.biSizeImage * sizeof(double));
        for (size_t p = 0; p < InfoHeader.biSizeImage; ++p) {
            new_Iz_RGB_2[p] = 0;
        }
        uint32_t i;
        double sharp_matr[3][3] = {{0,  -1, 0},
                                   {-1, 5,  -1},
                                   {0,  -1, 0}};
        for (int y = 0; y < Height(); ++y) {
            for (int x = 0; x < Width(); ++x) {
                Matrix_Multiplication(Iz_RGB, x, y, &(sharp_matr[0][0]), 3, new_Iz_RGB_2, InfoHeader.biSizeImage,
                                      stride);
                i = y * stride + x * bytes_per_pixel;
            }
        }
        for (auto j = 0; j < InfoHeader.biSizeImage; ++j) {
            if (new_Iz_RGB_2[j] > 255) new_Iz_RGB_2[j] = 255;
            if (new_Iz_RGB_2[j] < 0) new_Iz_RGB_2[j] = 0;
            new_Iz_RGB[j] = (unsigned char) (new_Iz_RGB_2[j]);
        }
        free(new_Iz_RGB_2);
        auto *f = new BMPFile(*this);
        f->SetNewIzRGB(new_Iz_RGB);
        return f;
    }


    BMPFile *Edge_Detection(const double threshold) {
        auto *new_Iz_RGB = (unsigned char *) malloc(InfoHeader.biSizeImage);
        auto *new_Iz_RGB_2 = (double *) malloc(InfoHeader.biSizeImage * sizeof(double));
        for (size_t p = 0; p < InfoHeader.biSizeImage; ++p) {
            new_Iz_RGB_2[p] = 0;
            new_Iz_RGB[p] = 0;
        }
        if ((stride & 3) != 0) stride = (stride & 0xFFFFFFFC) + 4;
        uint32_t i;
        double edge_matr[3][3] = {{0,  -1, 0},
                                  {-1, 4,  -1},
                                  {0,  -1, 0}};
        for (int y = 0; y < Height(); ++y) {
            for (int x = 0; x < Width(); ++x) {
                Matrix_Multiplication(Iz_RGB, x, y, &(edge_matr[0][0]), 3, new_Iz_RGB_2, InfoHeader.biSizeImage,
                                      stride);
                i = y * stride + x * bytes_per_pixel;
                if (new_Iz_RGB_2[i] > threshold) {
                    new_Iz_RGB[i] = 255;   // белок
                } else {
                    new_Iz_RGB[i] = 0;     // чернок
                }
            }
        }
        for (auto j = 0; j < InfoHeader.biSizeImage; ++j) {
            if (new_Iz_RGB_2[j] > 255) new_Iz_RGB_2[j] = 255;
            if (new_Iz_RGB_2[j] < 0) new_Iz_RGB_2[j] = 0;
            new_Iz_RGB[j] = (unsigned char) (new_Iz_RGB_2[j]);
        }
        free(new_Iz_RGB_2);
        auto *f = new BMPFile(*this);
        f->SetNewIzRGB(new_Iz_RGB);
        return f;
    }

    static double fade(double color, double numbofpoint, double kolvo, double extent) {
        double t = kolvo - numbofpoint + 1;
        if (t <= 0) {
            t = 1;
        }
        return extent * (color / t);
    }

    BMPFile *Light_Tunnel(double x0, double y0, double speed, double extent) {
        auto *new_Iz_RGB = (unsigned char *) malloc(InfoHeader.biSizeImage);
        auto *Iz_RGB_2 = (double *) malloc(InfoHeader.biSizeImage * sizeof(double));
        uint32_t i;
        for (int y = 0; y < Height(); ++y) {
            for (int x = 0; x < Width(); ++x) {
                i = y * stride + x * bytes_per_pixel;
                //Iz_RGB_2[i] = Iz_RGB[i];// ###7e315
                //Iz_RGB_2[i + 1] = Iz_RGB[i + 1];// ###7e315
                //Iz_RGB_2[i + 2] = Iz_RGB[i + 2];// ###7e315
                Iz_RGB_2[i] = 0;
                Iz_RGB_2[i + 1] = 0;
                Iz_RGB_2[i + 2] = 0;
            }
        }

        int flag = 0;
        auto *p = new Point(0, 0, true);
        auto *q = new Point(0.2 * speed, 0.2, true);
        auto *p2 = new Point(0, 0, true); // ###7e315
        auto *q2 = new Point(0.2 * speed, 0.2, true);// ###7e315
        auto *u = new Point(0, 0, true);
        double dphi = 0.1;
        int n = 20;
        double r;
        double b;
        double g;
        int index1;
        int index2;
        int xx;//### 7e315
        int yy;//### 7e315
        double weight; //### 7e315
        do {
            flag += 1;
            for (double rr = -speed / 2; rr < speed / 2; ++rr) {// ###7e315
                p2->phi = p->phi;
                q2->phi = q->phi;
                p2->ro = p2->phi * (speed + rr);// ###7e315
                q2->ro = q2->phi * (speed + rr);// ###7e315
                p2->MakeGrizzly();// ###7e315
                q2->MakeGrizzly();// ###7e315
                xx = (int) q2->x + x0;//### 7e315
                yy = (int) q2->y + y0;//### 7e315
                if ((xx < 0) || (xx >= InfoHeader.biWidth)) continue;//### 7e315
                index1 = yy * stride + xx * bytes_per_pixel;//### 7e315
                if ((index1 >= 0) && (index1 < InfoHeader.biSizeImage)) {
                    b = Iz_RGB[index1];
                    g = Iz_RGB[index1 + 1];
                    r = Iz_RGB[index1 + 2];
                    u->phi = p2->phi;
                    u->ro = u->phi * (speed + rr);
                    u->MakeGrizzly();//### 7e315
                    weight = 1; //### 7e315
                    for (size_t i = 0; i < n; ++i) {
                        xx = (int) u->x + x0;//### 7e315
                        yy = (int) u->y + y0;//### 7e315
                        if ((xx >= 0) || (xx < InfoHeader.biWidth)) {//### 7e315
                            index2 = yy * stride + xx * bytes_per_pixel;//### 7e315
                            if ((index2 >= 0) && (index2 < InfoHeader.biSizeImage)) {
                                b += fade(Iz_RGB[index2], i, n, extent);
                                g += fade(Iz_RGB[index2 + 1], i, n, extent);
                                r += fade(Iz_RGB[index2 + 2], i, n, extent);
                                weight += fade(1, i, n, extent);//### 7e315
                            }
                        }
                        u->phi += dphi;
                        u->ro = u->phi * (speed + rr);
                        u->MakeGrizzly();//### 7e315
                    }
                    Iz_RGB_2[index1] = b / weight;//### 7e315
                    Iz_RGB_2[index1 + 1] = g / weight;//### 7e315
                    Iz_RGB_2[index1 + 2] = r / weight;//### 7e315
                    //Iz_RGB_2[index1] = 0;
                    //Iz_RGB_2[index1 + 1] = 0;
                    //Iz_RGB_2[index1 + 2] = 192;
                    flag = 0;
                }
            }
            p->phi += dphi;
            p->ro = p->phi * speed;// ###7e315
            p->MakeGrizzly();
            q->phi += dphi;
            q->ro = q->phi * speed;// ###7e315
            q->MakeGrizzly();
        } while (flag < 1000);

        //new_Iz_RGB[i] = (unsigned char) std::round(new_Iz_RGB_2[i]);
        //new_Iz_RGB[i + 1] = (unsigned char) std::round(new_Iz_RGB_2[i + 1]);
        //new_Iz_RGB[i + 2] = (unsigned char) std::round(new_Iz_RGB_2[i + 2]);

        for (auto j = 0; j < InfoHeader.biSizeImage; ++j) {
            if (Iz_RGB_2[j] > 255) Iz_RGB_2[j] = 255;
            new_Iz_RGB[j] = (unsigned char) (Iz_RGB_2[j]);
        }
        free(Iz_RGB_2);
        auto *f = new BMPFile(*this);
        f->SetNewIzRGB(new_Iz_RGB);
        return f;
    }
};



//void ReadBMP(char *input_file_name) {}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "неправильный ввод" << "\n";
        return -1;
    }
    uint32_t w;
    uint32_t h;
    auto *r = new BMPFile(argv[1]);
    BMPFile *r0;
    uint32_t number = 3;
    double t;
    double x0;
    double y0;
    double speed;
    double extent;
    double lag;
    double delag;
    char *p; // отслеживаем ошибки
    //считывает фильтры и действия
    while (number < argc) {
        if (std::strcmp(argv[number], "-gs") == 0) {
            r0 = r;
            r = r0->Grayscale();
            delete r0;
            ++number;
        } else if (std::strcmp(argv[number], "-crop") == 0) {
            ++number;
            if (number == argc) {
                break;
            }
            w = std::strtoul(argv[number], &p, 10);
            ++number;
            if (number == argc) {
                break;
            }
            h = std::strtoul(argv[number], &p, 10);
            ++number;
            r0 = r;
            r = r0->Crop(w, h);
            delete r0;
        } else if (std::strcmp(argv[number], "-neg") == 0) {
            r0 = r;
            r = r0->Negative();
            delete r0;
            ++number;
        } else if (std::strcmp(argv[number], "-blur") == 0) {
            r0 = r;
            ++number;
            t = strtod(argv[number], &p);
            r = r0->Gaussian_Blur(t);
            delete r0;
            ++number;
        } else if (std::strcmp(argv[number], "-sharp") == 0) {
            r0 = r;
            r = r0->Sharpening();
            delete r0;
            ++number;
        } else if (std::strcmp(argv[number], "-lt") == 0) {
            r0 = r;
            ++number;
            x0 = strtod(argv[number], &p);
            ++number;
            y0 = strtod(argv[number], &p);
            ++number;
            speed = strtod(argv[number], &p);
            ++number;
            extent = strtod(argv[number], &p);
            r = r0->Light_Tunnel(x0, y0, speed, extent);
            delete r0;
            ++number;
        } else if (std::strcmp(argv[number], "-edge") == 0) {
            r0 = r->Grayscale();
            ++number;
            t = strtod(argv[number], &p);
            r = r0->Edge_Detection(t);
            delete r0;
            ++number;
        }

        r->Export(argv[2]);
        w = r->Width();
        h = r->Height();
        std::cout << '\n';
        std::cout << w << ' ' << h << ' ' << r->BitCount();
        delete r;
        return 0;
    }
}