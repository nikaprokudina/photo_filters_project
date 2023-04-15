#include "image.h"

uint32_t BMPFile::Height() const {
    return info_header_.biHeight;
}

uint32_t BMPFile::Width() const {
    return info_header_.biWidth;
}

uint32_t BMPFile::BitCount() const {
    return info_header_.biBitCount;
}

void BMPFile::SetNewIzRGB(unsigned char *new_iz_rgb) {
    delete[] iz_rgb_;
    iz_rgb_ = new unsigned char[static_cast<size_t>(info_header_.biSizeImage)];
    std::copy(new_iz_rgb, new_iz_rgb + info_header_.biSizeImage, iz_rgb_);
}

void BMPFile::SetDimensions(uint32_t new_w, uint32_t new_h) {
    info_header_.biWidth = static_cast<LONG>(new_w);
    info_header_.biHeight = static_cast<LONG>(new_h);
    bytes_per_pixel_ = BitCount() / count_bytes_;
    stride_ = new_w * bytes_per_pixel_;
    if ((stride_ & 3) != 0) {
        stride_ = (stride_ & dvoichn_) + 4;
    }
    uint32_t new_bi_size_image = stride_ * new_h;
    file_header_.bfSize += (new_bi_size_image - info_header_.biSizeImage);
    info_header_.biSizeImage = new_bi_size_image;
}

BMPFile::BMPFile(BMPFile &other) {
    char *t = reinterpret_cast<char *>(&file_header_);
    char *u = reinterpret_cast<char *>(&(other.file_header_));
    for (int i = 0; i < sizeof(file_header_); ++i) {
        t[i] = u[i];
    }
    t = reinterpret_cast<char *>(&info_header_);
    u = reinterpret_cast<char *>(&(other.info_header_));
    for (int i = 0; i < sizeof(info_header_); ++i) {
        t[i] = u[i];
    }
    iz_rgb_ = NULL;
    stride_ = other.stride_;
    bytes_per_pixel_ = other.bytes_per_pixel_;
}

BMPFile::BMPFile(char *input_file_name) {
    FILE *stream = nullptr;
    stream = fopen(input_file_name, "rb");
    fread(&file_header_, sizeof(file_header_), 1, stream);
    fread(&info_header_, sizeof(info_header_), 1, stream);
    stride_ = info_header_.biWidth * (info_header_.biBitCount / count_bytes_);
    if ((stride_ & 3) != 0) {
        stride_ = (stride_ & dvoichn_) + 4;
    }
    if (info_header_.biSizeImage == 0) {  // ### 7e315
        info_header_.biSizeImage = info_header_.biHeight * stride_;
    }
    fseek(stream, file_header_.bfOffBits, SEEK_SET);  // указатель позиции в файле (сейчас стоит в начале)
    iz_rgb_ =
        new unsigned char[static_cast<size_t>(info_header_.biSizeImage)];  // выделяет блок памяти под размер файла
    fread(iz_rgb_, info_header_.biSizeImage, 1, stream);
    fclose(stream);
    bytes_per_pixel_ = BitCount() / count_bytes_;
    padding_bytes_ = (4 - (Width() * bytes_per_pixel_) % 4) % 4;
}

BMPFile::~BMPFile() {  // destructor
    if (iz_rgb_ != nullptr) {
        delete[] iz_rgb_;
    }
}

void BMPFile::Export(char *out_file_name) {
    FILE *outfile = fopen(out_file_name, "wb");
    fwrite(&file_header_, sizeof(char), sizeof(BitmapFileHeader), outfile);
    fwrite(&info_header_, sizeof(char), sizeof(BitmapInfoHeader), outfile);
    fwrite(iz_rgb_, sizeof(unsigned char), info_header_.biSizeImage, outfile);
    fclose(outfile);
}

BMPFile *BMPFile::Grayscale() {
    auto *new_iz_rgb = new unsigned char[static_cast<size_t>(info_header_.biSizeImage)];
    // std::copy(iz_rgb_, iz_rgb_ + info_header_.biSizeImage, new_iz_rgb);
    uint64_t i = 0;
    const double amount_of_red = 0.299;
    const double amount_of_green = 0.587;
    const double amount_of_blue = 0.114;
    const double max_color = 255;
    for (int y = 0; y < Height(); ++y) {
        for (int x = 0; x < Width(); ++x) {
            i = y * stride_ + x * bytes_per_pixel_;
            unsigned char b = iz_rgb_[i];
            unsigned char g = iz_rgb_[i + 1];
            unsigned char r = iz_rgb_[i + 2];
            auto gray = static_cast<unsigned char>(
                std::clamp(amount_of_red * r + amount_of_green * g + amount_of_blue * b, 0.0, max_color));
            new_iz_rgb[i] = gray;
            new_iz_rgb[i + 1] = gray;
            new_iz_rgb[i + 2] = gray;
            if (bytes_per_pixel_ == 4) {
                new_iz_rgb[i + 3] = iz_rgb_[i + 3];
            }
        }
    }
    // std::copy(new_iz_rgb, new_iz_rgb + info_header_.biSizeImage, iz_rgb_);
    auto *f = new BMPFile(*this);
    f->SetNewIzRGB(new_iz_rgb);
    delete[] new_iz_rgb;
    return f;
}

BMPFile *BMPFile::Negative() {
    auto *new_iz_rgb = new unsigned char[static_cast<size_t>(info_header_.biSizeImage)];
    // std::copy(iz_rgb_, iz_rgb_ + info_header_.biSizeImage, new_iz_rgb);
    uint32_t i = 0;
    const int max_color = 255;
    for (int y = 0; y < Height(); ++y) {
        for (int x = 0; x < Width(); ++x) {
            i = y * stride_ + x * bytes_per_pixel_;
            unsigned char b = iz_rgb_[i];
            unsigned char g = iz_rgb_[i + 1];
            unsigned char r = iz_rgb_[i + 2];
            new_iz_rgb[i] = static_cast<unsigned char>(max_color - b);      // b
            new_iz_rgb[i + 1] = static_cast<unsigned char>(max_color - g);  // g
            new_iz_rgb[i + 2] = static_cast<unsigned char>(max_color - r);  // r
        }
    }
    // std::copy(new_iz_rgb, new_iz_rgb + info_header_.biSizeImage, iz_rgb_);
    auto *f = new BMPFile(*this);
    f->SetNewIzRGB(new_iz_rgb);
    delete[] new_iz_rgb;
    return f;
}

// BMPFile *BMPFile::Crop(uint32_t width, uint32_t height) {
//     if (width > info_header_.biWidth) {
//         width = info_header_.biWidth;
//     }
//     if (height > info_header_.biHeight) {
//         height = info_header_.biHeight;
//     }
//     uint32_t new_stride = width * bytes_per_pixel_;
//     if ((new_stride & 3) != 0) {
//         new_stride = (new_stride & dvoichn_) + 4;
//     }
//     uint32_t new_bi_size_image = new_stride * height;
//     auto *new_iz_rgb = new unsigned char[new_bi_size_image];
//     uint32_t i = 0;
//     uint32_t j = 0;
//     uint32_t y2 = 0;
//     for (uint32_t y = info_header_.biHeight - height; y < info_header_.biHeight; ++y) {
//         for (uint32_t x = 0; x < width; ++x) {
//             i = y * stride_ + x * bytes_per_pixel_;
//             j = y2 * new_stride + x * bytes_per_pixel_;
//             for (uint32_t k = 0; k < bytes_per_pixel_; ++k) {
//                 new_iz_rgb[j + k] = iz_rgb_[i + k];
//             }
//         }
//         ++y2;
//     }
//     auto *f = new BMPFile(*this);
//     f->SetNewIzRGB(new_iz_rgb);
//     f->SetDimensions(width, height);
//     return f;
// }

// BMPFile *BMPFile::GaussianBlur(double sigma) {
//     auto *new_iz_rgb = new unsigned char[info_header_.biSizeImage];
//     auto *new_iz_rgb_2 = new double[info_header_.biSizeImage];
//     uint32_t i = 0;
//     double xx = 0;
//     double yy = 0;
//     double t = 1 / (2 * sigma * sigma);
//     double k = t / M_PI;
//     double mq = 0;
//     const int max_porog = 6;
//     double porog = max_porog - log(k);
//     double p = 0;
//     uint32_t x_min = 0;
//     uint32_t x_max = 0;
//     uint32_t y_min = 0;
//     uint32_t y_max = 0;
//     for (int y = 0; y < Height(); ++y) {
//         for (int x = 0; x < Width(); ++x) {
//             i = y * stride_ + x * bytes_per_pixel_;
//             if (new_iz_rgb_2 != nullptr) {
//                 new_iz_rgb_2[i] = 0;
//                 new_iz_rgb_2[i + 1] = 0;
//                 new_iz_rgb_2[i + 2] = 0;
//             }
//             x_min = x - 3;
//             if (x - 3 < 0) {
//                 x_min = 0;
//             }
//             y_min = y - 3;
//             if (y - 3 < 0) {
//                 y_min = 0;
//             }
//             x_max = x + 3;
//             if (x_max >= Width()) {
//                 x_max = Width() - 1;
//             }
//             y_max = y + 3;
//             if (y_max >= Height()) {
//                 y_max = Height() - 1;
//             }
//             for (auto y1 = y_min; y1 <= y_max; ++y1) {
//                 yy = y1 - y;
//                 yy *= yy;
//                 for (auto x1 = x_min; x1 <= x_max; ++x1) {
//                     xx = x1 - x;
//                     xx *= xx;
//                     p = (xx + yy) * t;
//                     if (p > porog) {
//                         continue;
//                     }
//                     mq = k * exp(-p);
//                     uint32_t ii = (y1 * (Width() * bytes_per_pixel_ + padding_bytes_)) + (x1 * bytes_per_pixel_);
//                     unsigned char b = iz_rgb_[ii];
//                     unsigned char g = iz_rgb_[ii + 1];
//                     unsigned char r = iz_rgb_[ii + 2];
//                     new_iz_rgb_2[i] += (b * mq);      // b
//                     new_iz_rgb_2[i + 1] += (g * mq);  // g
//                     new_iz_rgb_2[i + 2] += (r * mq);  // r
//                 }
//             }
//             new_iz_rgb[i] = static_cast<unsigned char>(std::round(new_iz_rgb_2[i]));
//             new_iz_rgb[i + 1] = static_cast<unsigned char>(std::round(new_iz_rgb_2[i + 1]));
//             new_iz_rgb[i + 2] = static_cast<unsigned char>(std::round(new_iz_rgb_2[i + 2]));
//         }
//     }
//     delete[]new_iz_rgb_2;
//     auto *f = new BMPFile(*this);
//     f->SetNewIzRGB(new_iz_rgb);
//     return f;
// }

std::unique_ptr<BMPFile> BMPFile::Sharpening() {
    auto *new_iz_rgb = new unsigned char[static_cast<size_t>(info_header_.biSizeImage)];
    auto *new_iz_rgb_2 = new double[static_cast<size_t>(info_header_.biSizeImage)];
    for (size_t p = 0; p < info_header_.biSizeImage; ++p) {
        new_iz_rgb_2[p] = 0;
    }
    constexpr int central_matrix = 5;
    constexpr double sharp_matr[3][3] = {{0, -1, 0}, {-1, central_matrix, -1}, {0, -1, 0}};

    for (int y = 1; y < Height() - 1; ++y) {
        for (int x = 1; x < Width() - 1; ++x) {
            MatrixMultiplication(iz_rgb_, x, y, &sharp_matr[0][0], 3, &new_iz_rgb_2[0],
                                 static_cast<int32_t>(info_header_.biSizeImage), static_cast<int32_t>(stride_), Width(),
                                 Height());
        }
    }

    const int max_color = 255;
    std::transform(new_iz_rgb_2.begin(), new_iz_rgb_2.end(), new_iz_rgb.begin(), [max_color](double v) {
        return static_cast<unsigned char>(std::clamp(v, 0.0, static_cast<double>(max_color)));
    });

    auto f = std::make_unique<BMPFile>(*this);
    f->SetNewIzRGB(new_iz_rgb.data());
    return f;
}

BMPFile *BMPFile::EdgeDetection(const double threshold) {
    auto *new_iz_rgb = new unsigned char[static_cast<size_t>(info_header_.biSizeImage)];
    auto *new_iz_rgb_2 = new double[static_cast<size_t>(info_header_.biSizeImage)];
    for (size_t p = 0; p < info_header_.biSizeImage; ++p) {
        new_iz_rgb_2[p] = 0;
        new_iz_rgb[p] = 0;
    }
    if ((stride_ & 3) != 0) {
        stride_ = (stride_ & dvoichn_) + 4;
    }
    uint32_t i = 0;
    double edge_matr[3][3] = {{0, -1, 0}, {-1, 4, -1}, {0, -1, 0}};
    const int max_color = 255;
    for (int y = 0; y < Height(); ++y) {
        for (int x = 0; x < Width(); ++x) {
            MatrixMultiplication(iz_rgb_, x, y, &(edge_matr[0][0]), 3, new_iz_rgb_2,
                                 static_cast<int32_t>(info_header_.biSizeImage), static_cast<int32_t>(stride_), Width(),
                                 Height());
            i = y * stride_ + x * bytes_per_pixel_;
            if (new_iz_rgb_2[i] > threshold) {
                new_iz_rgb[i] = max_color;  // белок
            } else {
                new_iz_rgb[i] = 0;  // чернок
            }
        }
    }
    for (auto j = 0; j < info_header_.biSizeImage; ++j) {
        if (new_iz_rgb_2[j] > max_color) {
            new_iz_rgb_2[j] = max_color;
        }
        if (new_iz_rgb_2[j] < 0) {
            new_iz_rgb_2[j] = 0;
        }
        new_iz_rgb[j] = static_cast<unsigned char>(new_iz_rgb_2[j]);
    }
    std::copy(new_iz_rgb, new_iz_rgb + info_header_.biSizeImage, iz_rgb_);
    delete[] new_iz_rgb_2;
    auto *f = new BMPFile(*this);
    f->SetNewIzRGB(new_iz_rgb);
    delete[] new_iz_rgb;
    return f;
}

BMPFile *BMPFile::Space(double x0, double y0, double speed, double extent) {
    auto *new_iz_rgb = new unsigned char[info_header_.biSizeImage];
    auto *new_iz_rgb_2 = new double[info_header_.biSizeImage];
    uint32_t i = 0;
    for (int y = 0; y < Height(); ++y) {
        for (int x = 0; x < Width(); ++x) {
            i = y * stride_ + x * bytes_per_pixel_;
            // Iz_RGB_2[i] = Iz_RGB[i];// ###7e315
            // Iz_RGB_2[i + 1] = Iz_RGB[i + 1];// ###7e315
            // Iz_RGB_2[i + 2] = Iz_RGB[i + 2];// ###7e315
            new_iz_rgb_2[i] = 0;
            new_iz_rgb_2[i + 1] = 0;
            new_iz_rgb_2[i + 2] = 0;
        }
    }

    int flag = 0;
    const int max_flag = 1000;
    const int max_color = 255;
    const double znach_q = 0.2;
    const double znach_dphi = 0.1;
    const int znach_n = 20;
    auto *p = new Point(0, 0, true);
    auto *q = new Point(znach_q * speed, znach_q, true);
    auto *p2 = new Point(0, 0, true);
    auto *q2 = new Point(znach_q * speed, znach_q, true);
    auto *u = new Point(0, 0, true);
    double dphi = znach_dphi;
    int n = znach_n;
    double r = 0;
    double b = 0;
    double g = 0;
    uint32_t index1 = 0;
    uint32_t index2 = 0;
    int xx = 0;
    int yy = 0;
    double weight = 0;
    do {
        flag += 1;
        for (double rr = -speed / 2; rr < speed / 2; ++rr) {
            p2->phi = p->phi;
            q2->phi = q->phi;
            p2->ro = p2->phi * (speed + rr);
            q2->ro = q2->phi * (speed + rr);
            p2->MakeGrizzly();
            q2->MakeGrizzly();
            xx = static_cast<int>(q2->x + x0);
            yy = static_cast<int>(q2->y + y0);
            if ((xx < 0) || (xx >= info_header_.biWidth)) {
                continue;
            }
            index1 = yy * stride_ + xx * bytes_per_pixel_;
            if ((yy * stride_ + xx * bytes_per_pixel_ >= 0) && (index1 < info_header_.biSizeImage)) {
                b = iz_rgb_[index1];
                g = iz_rgb_[index1 + 1];
                r = iz_rgb_[index1 + 2];
                u->phi = p2->phi;
                u->ro = u->phi * (speed + rr);
                u->MakeGrizzly();
                weight = 1;
                for (auto h = 0; h < n; ++h) {
                    xx = static_cast<int>(u->x + x0);
                    yy = static_cast<int>(u->y + y0);
                    if ((xx >= 0) || (xx < info_header_.biWidth)) {
                        index2 = yy * stride_ + xx * bytes_per_pixel_;
                        if ((yy * stride_ + xx * bytes_per_pixel_ >= 0) && (index2 < info_header_.biSizeImage)) {
                            b += FadeF(iz_rgb_[index2], h, n, extent);
                            g += FadeF(iz_rgb_[index2 + 1], h, n, extent);
                            r += FadeF(iz_rgb_[index2 + 2], h, n, extent);
                            weight += FadeF(1, h, n, extent);  // ### 7e315
                        }
                    }
                    u->phi += dphi;
                    u->ro = u->phi * (speed + rr);
                    u->MakeGrizzly();
                }
                new_iz_rgb_2[index1] = b / weight;
                new_iz_rgb_2[index1 + 1] = g / weight;
                new_iz_rgb_2[index1 + 2] = r / weight;
                // Iz_RGB_2[index1] = 0;
                // Iz_RGB_2[index1 + 1] = 0;
                // Iz_RGB_2[index1 + 2] = 192;
                flag = 0;
            }
        }
        p->phi += dphi;
        p->ro = p->phi * speed;  // ###7e315
        p->MakeGrizzly();
        q->phi += dphi;
        q->ro = q->phi * speed;  // ###7e315
        q->MakeGrizzly();
    } while (flag < max_flag);

    // new_iz_rgb[i] = (unsigned char) std::round(new_iz_rgb_2[i]);
    // new_iz_rgb[i + 1] = (unsigned char) std::round(new_iz_rgb_2[i + 1]);
    // new_iz_rgb[i + 2] = (unsigned char) std::round(new_iz_rgb_2[i + 2]);

    for (auto j = 0; j < info_header_.biSizeImage; ++j) {
        if (new_iz_rgb_2[j] > max_color) {
            new_iz_rgb_2[j] = max_color;
        }
        new_iz_rgb[j] = static_cast<unsigned char>(new_iz_rgb_2[j]);
    }
    std::copy(new_iz_rgb, new_iz_rgb + info_header_.biSizeImage, iz_rgb_);
    delete[] new_iz_rgb_2;
    auto *f = new BMPFile(*this);
    f->SetNewIzRGB(new_iz_rgb);
    delete[] new_iz_rgb;
    return f;
}