#include "util.h"

void Point::MakePolar() {
    ro = sqrt(x * x + y * y);
    if (x != 0) {
        phi = atan(y / x);
    } else {
        phi = 0;
    }
}

void Point::MakeGrizzly() {
    x = ro * cos(phi);
    y = ro * sin(phi);
}

Point::Point(double a, double b, bool ispolar) {
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

double FadeF(double color, double numbofpoint, double kolvo, double extent) {
    double t = kolvo - numbofpoint + 1;
    if (t <= 0) {
        t = 1;
    }
    return extent * (color / t);
}

void MatrixMultiplication(const unsigned char *iz, int32_t x, int32_t y, double *matrix, int32_t w, double *iz2,
                          int32_t size, int32_t stride, uint32_t image_width, uint32_t image_height) {
    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < w; ++j) {
            int32_t new_x =
                static_cast<int>(std::fmax(0.0, std::fmin(static_cast<double>(image_width - 1), x + (j - w / 2))));
            int32_t new_y =
                static_cast<int>(std::fmax(0.0, std::fmin(static_cast<double>(image_height - 1), y + (i - w / 2))));
            for (int color = 0; color < 3; ++color) {
                int32_t p1 = y * stride + x * 3 + color;
                int32_t p2 = new_y * stride + new_x * 3 + color;
                iz2[p1] += matrix[i * w + j] * iz[p2];
            }
        }
    }
}

int StrCmp(const char *str1, const char *str2) {
    while (*str1 && (*str1 == *str2)) {
        str1++;
        str2++;
    }
    return static_cast<const int>(*reinterpret_cast<const unsigned char *>(str1)) -
           static_cast<const int>(*reinterpret_cast<const unsigned char *>(str2));
}
// bool StrCmp(const char *a, const char *b) {
// if (!a || !b) {
// return false;
//}
// while (true) {
// if (*a != *b) {
// return false;
//}
// if ((*a == 0) || (*b == 0)) {
// return true;
//}
//++a;
//++b;
//}
//}