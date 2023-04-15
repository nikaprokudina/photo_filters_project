#include <cmath>
#include <math.h>
#include <cstdint>

class Point {
public:
    double x;
    double y;
    double ro;
    double phi;
    uint32_t image_width;
    uint32_t image_height;

    void MakePolar();

    void MakeGrizzly();

    Point(double a, double b, bool ispolar);
};

double FadeF(double color, double numbofpoint, double kolvo, double extent);

void MatrixMultiplication(const unsigned char *iz, int32_t x, int32_t y, double *matrix, int32_t w, double *iz2,
                          int32_t size, int32_t stride, uint32_t image_width, uint32_t image_height);
// bool StrCmp(const char *a, const char *b);
int StrCmp(const char *str1, const char *str2);