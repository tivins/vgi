#ifndef _vgi_tuple_h_
#define _vgi_tuple_h_

#include <array>
#include <string>
#include <ostream>
#include <iostream>
#include <fstream>
#include <cmath>
// #include <stdarg.h>
#include <cxxabi.h> // demangled type names

// ---------------------------------------

namespace vgi {

inline double real(double r) { return r; }
inline double to_radians(const double& a) { return a * M_PI / 180; }

// ---------------------------------------
enum Tuple_Type {
    kTuple, kPoint, kDirection, kSize, kColor, kLast
    };
static const char * Tuple_Type_Names[kLast] = {
    "Tuple", "Point", "Direction", "Size", "Color"
    };

template <typename T>
std::string type_name()
{
    int status;
    std::string tname = typeid(T).name();
    char *demangled_name = abi::__cxa_demangle(tname.c_str(), NULL, NULL, &status);
    if (status == 0)
    {
        tname = demangled_name;
        std::free(demangled_name);
    }
    return tname;
}

// ---------------------------------------
// Data Structures
// ---------------------------------------

/**
 * Linear vector of TDim dimension of a TType values.
 * TClass is used to enclose the tuple's types.
 */
template <int TDim, typename TType, Tuple_Type TClass>
struct Tuple {
    std::array<TType, TDim> data;
    Tuple() { data.fill(static_cast<TType>(0)); }
    explicit Tuple(const TType &a) { data.fill(a); }
    ~Tuple() { std::cout<<"~"<<*this<<std::endl; }
    Tuple(const TType &a, const TType &b) { data[0] = a; data[1] = b; }
    Tuple(const TType &a, const TType &b, const TType &c) { data[0] = a; data[1] = b; data[2] = c; }
    Tuple(const TType &a, const TType &b, const TType &c, const TType &d) { data[0] = a; data[1] = b; data[2] = c; data[3] = d; }
    inline const TType& at(int idx) const { return data.at(idx); }
    inline const TType& operator[](int idx) const { return data.at(idx); }
};

/**
 * Square matrix
 */
template <int TDim, typename TType>
struct Matrix {
    std::array<Tuple<TDim, TType, kTuple>, TDim> data;
    Matrix() {
        for (int i=0;i<TDim;i++)
            data[i].data.fill(static_cast<TType>(0));
    }
    explicit Matrix(const TType& val) {
        for (int i=0;i<TDim;i++)
            data[i].data.fill(val);
    }
    inline const Tuple<TDim, TType, kTuple>& row(int row) const {
        return data[row] ;
    }
    inline Tuple<TDim, TType, kTuple> col(int col) const {
        Tuple<TDim, TType, kTuple> out;
        for (int row = 0; row < TDim; row++)
            out.data[row] = data[row].data[col];
        return out;
    }
    inline void identity(TType val0, TType val1) {
        for (int row = 0; row < TDim; row++) {
            data[row].data.fill(val0);
            data[row].data[row] = val1;
        }
    }
};

template <int TDim, typename TType>
struct Ray {
    Tuple<TDim, TType, kPoint> origin;
    Tuple<TDim, TTYpe, kDirection> direction;
    Ray() : origin(0), direction(0,0,1) { }
    Ray(Tuple<TDim, TType, kPoint> origin, Tuple<TDim, TTYpe, kDirection> direction) : origin(origin), direction(direction) { }
};

template <int TDim, typename TType>
struct Triangle {
    Tuple<TDim, TType, kPoint> vertices[3];
};

const int Alpha = 3;

template <int TDim, int TDim2, typename TType>
Tuple<TDim, TType, kColor> mix(const Tuple<TDim, TType, kColor>& under, const Tuple<TDim2, TType, kColor>& upper) {
    if (TDim2 == 3) return Tuple<TDim, TType, kColor>(
            upper[0],
            upper[1],
            upper[2]
        );
    if (TDim2 != 4) return Tuple<TDim, TType, kColor>(
            upper[0],
            upper[1],
            upper[2]
        );
    if (upper[Alpha] == static_cast<TType>(1))
        return Tuple<TDim, TType, kColor>(
            upper[0],
            upper[1],
            upper[2]
        );
    return Tuple<TDim, TType, kColor>(
        under[0] * (static_cast<TType>(1) - upper[Alpha]) + upper[0] * (upper[Alpha]),
        under[1] * (static_cast<TType>(1) - upper[Alpha]) + upper[1] * (upper[Alpha]),
        under[2] * (static_cast<TType>(1) - upper[Alpha]) + upper[2] * (upper[Alpha])
    );
}

/**
 * Image
 */
template<int TDim, typename TType>
struct Image {
    Tuple<2, int, kSize> size;
    Tuple<TDim, TType, kColor> * data = nullptr;
    void allocate(const Tuple<2, int, kSize>& size) {
        this->size = size;
        data = new Tuple<TDim, TType, kColor>[area()];
    }
    Image() { }
    ~Image() {
        if (data != nullptr) {
            delete [] data;
        }
    }
    inline size_t area() const { return size[0] * size[1]; }

    template <int extdim>
    void set(const Tuple<2, int, kPoint>& at, const Tuple<extdim, TType, kColor>& color) {
        size_t index = at[1] * size[0] + at[0];
        data[index] = mix(data[index], color);
    }
};

// ---------------------------------------
// Methods
// ---------------------------------------

template<int TDim, typename TType, typename TTypeOut>
Tuple<TDim, TTypeOut, kColor> convert(const Tuple<TDim, TType, kColor>& source) {
    Tuple<TDim, uint8_t, kColor> out;
    for (int it=0;it<TDim;it++)
        out.data[it] = static_cast<TTypeOut>(source[it]);
    return out;
}

template<int TDim, typename TType>
void convert(const Image<TDim, TType>& image, Image<TDim, uint8_t>& out) {
    out.allocate(image.size);
    for (size_t it = 0; it < image.area(); it++) {
        out.data[it] = convert<3, double, uint8_t>(image.data[it] * 255.0);
    }
}

void save(const Image<3, uint8_t>& image, const char * filename) {
    std::ofstream file;
    file.open(filename);
    file << "P" << 3 << "\n\n";
    file << image.size[0] << " " << image.size[1] << "\n255\n";
    for (size_t it = 0; it < image.area(); it++) {
        file << (int)image.data[it][0] << " " << (int)image.data[it][1] << " " << (int)image.data[it][2] << "\n";
    }
    file.close();
}

/**
 * Dot product between 2 directions.
 */
template <int TDim, typename TType>
inline TType dot(const Tuple<TDim, TType, kDirection>& a, const Tuple<TDim, TType, kDirection>& b) {
    TType out = static_cast<TType>(0);
    for (int i = 0; i < TDim; i++)
        out += a[i] * b[i];
    return out;
}

/**
 * Cross product between 2 directions with 3 dimensions.
 */
template <typename TType>
inline Tuple<3, TType, kDirection> cross(const Tuple<3, TType, kDirection>& a, const Tuple<3, TType, kDirection>& b) {
    return Tuple<3, TType, kDirection>(a[1]*b[2]-a[2]*b[1], a[0]*b[2]+a[2]*b[0], a[0]*b[1]-a[1]*b[0]);
}

/**
 * Length of a direction.
 */
template <int TDim, typename TType>
inline TType length(const Tuple<TDim, TType, kDirection>& v) {
    return sqrt(dot(v, v));
}

/**
 * Normalize a direction.
 */
template <int TDim, typename TType>
inline Tuple<TDim, TType, kDirection> normalize(const Tuple<TDim, TType, kDirection>& v) {
    double length = sqrt(dot(v, v));
    if (length == 0) return v;
    return v / length;
}

/**
 * Gets the sum of all components.
 */
template <int TDim, typename TType, Tuple_Type TClass>
inline TType sum(const Tuple<TDim, TType, TClass>& tpl) {
    TType out = static_cast<TType>(0);
    for (int i = 0; i < TDim; i++)
        out += tpl.data[i];
    return out;
}

/**
 * Gets the average value of all components.
 */
template <int TDim, typename TType, Tuple_Type TClass>
inline TType average(const Tuple<TDim, TType, TClass>& value) {
    return sum(value) / static_cast<double>(TDim); // double? or TType. If double, change type of return ?
}

/**
 * Negate a Tuple
 */
template <int TDim, typename TType, Tuple_Type TClass>
inline Tuple<TDim, TType, TClass> operator -(const Tuple<TDim, TType, TClass>& a)
{
    Tuple<TDim, TType, TClass> out;
    for (int it = 0; it < TDim; it++) {
        out.data[it] = -a.at(it);
    }
    return out;
};

/**
 * Tuple x Tuple
 */
template <int TDim, typename TType, Tuple_Type TClass>
inline Tuple<TDim, TType, TClass> operator *(const Tuple<TDim, TType, TClass>& a, const Tuple<TDim, TType, TClass>& b)
{
    Tuple<TDim, TType, TClass> out;
    for (int it = 0; it < TDim; it++) {
        out.data[it] = a.at(it) * b.at(it);
    }
    return out;
};

/**
 * Tuple x Type
 */
template <int TDim, typename TType, Tuple_Type TClass>
inline Tuple<TDim, TType, TClass> operator *(const Tuple<TDim, TType, TClass>& a, const TType& b)
{
    Tuple<TDim, TType, TClass> out;
    for (int it = 0; it < TDim; it++) {
        out.data[it] = a.at(it) * b;
    }
    return out;
};

/**
 * Tuple / Type
 */
template <int TDim, typename TType, Tuple_Type TClass>
inline Tuple<TDim, TType, TClass> operator /(const Tuple<TDim, TType, TClass>& a, const TType& b)
{
    Tuple<TDim, TType, TClass> out;
    for (int it = 0; it < TDim; it++) {
        out.data[it] = a.at(it) / b;
    }
    return out;
};

/**
 * Direction = Point - Point;
 */
template <int TDim, typename TType>
inline Tuple<TDim, TType, kDirection> operator - (const Tuple<TDim, TType, kPoint>& a, const Tuple<TDim, TType, kPoint>& b) {
    Tuple<TDim, TType, kDirection> out;
    for (int td = 0; td < TDim; td++) {
        out.data[td] = a.at(td) - b.at(td);
    }
    return out;
}

/**
 * Matrix x Matrix
 */
template <int TDimMat, typename TType>
inline Matrix<TDimMat, TType> operator *(const Matrix<TDimMat, TType>& a, const Matrix<TDimMat, TType>& b)
{
    Matrix<TDimMat, TType> out;
    for (int it1 = 0; it1 < TDimMat; it1++) {
        out.data[it1] = a.row(it1) * b.col(it1);
    }
    return out;
}

/**
 * Matrix x Tuple
 */
template <int TDimMat, int TDimTuple, typename TType, Tuple_Type TClass>
inline Tuple<TDimTuple, TType, TClass> operator *(const Matrix<TDimMat, TType>& mat, const Tuple<TDimTuple, TType, TClass>& tuple)
{
    Tuple<TDimTuple, TType, TClass> out;
    for (int td = 0; td < TDimTuple; td++) {
        out.data[td] = static_cast<TType>(0);
        for (int td2 = 0; td2 < TDimTuple; td2++) {
            out.data[td] += mat.data[td2].data[td] * tuple.data[td2];
        }
    }
    return out;
}

/**
 * Returns a perspective matrix.
 */
template <typename TType>
inline Matrix<4, TType> make_perspective(TType fov, TType aspect, TType znear, TType zfar)
{
    // | w, 0,  0,  0 |
    // | 0, h,  0,  0 |
    // | 0, 0,  q, -1 |
    // | 0, 0, qn,  0 |

    TType yScale = 1.0 / tan(to_radians(fov) / static_cast<TType>(2));
    TType xScale = yScale / aspect;
    TType nearmfar = znear - zfar;

    Matrix<4, TType> mat;
    mat.data[0] = Tuple<4, TType, kTuple>(xScale, 0, 0, 0);
    mat.data[1] = Tuple<4, TType, kTuple>(0, yScale, 0, 0);
    mat.data[2] = Tuple<4, TType, kTuple>(0, 0, (zfar + znear) / nearmfar, -1);
    mat.data[3] = Tuple<4, TType, kTuple>(0, 0, 2 * zfar * znear / nearmfar, 0);
    return mat;
}

/**
 * Returns a view matrix.
 */
template <typename TType>
Matrix<4, TType> make_look_at(const Tuple<3, TType, kPoint> &eye,
                                const Tuple<3, TType, kPoint> &at,
                                const Tuple<3, TType, kDirection> &up)
{
    Tuple<3, TType, kDirection> zaxis = normalize(eye - at);
    Tuple<3, TType, kDirection> xaxis = cross(normalize(up), zaxis);
    Tuple<3, TType, kDirection> yaxis = cross(zaxis, xaxis);
    Tuple<3, TType, kPoint> O(0);

    Matrix<4, TType> mat;
    mat.data[0] = Tuple<4, TType, kTuple>(xaxis[0], yaxis[0], zaxis[0], 0);
    mat.data[1] = Tuple<4, TType, kTuple>(xaxis[1], yaxis[1], zaxis[1], 0);
    mat.data[2] = Tuple<4, TType, kTuple>(xaxis[2], yaxis[2], zaxis[2], 0);
    mat.data[3] = Tuple<4, TType, kTuple>(dot(xaxis, O - eye), dot(yaxis, O - eye), dot(zaxis, O - eye), 1);
    return mat;
}

// ---------------------------------------
// TO STRING
// ---------------------------------------
namespace term {
enum Code
{
    reset = 0,
    black = 30,
    red = 31,
    green = 32,
    blue = 34,
    yellow = 33,
    magenta = 35,
    cyan = 36,
    light_gray = 37,
    dark_gray = 90,
    light_red = 91,
    light_green = 92,
    light_yellow = 93,
    light_blue = 94,
    light_magenta = 95,
    light_cyan = 96,
    white = 97,
};
}; // namespace term

static bool term_color_activated = true;
std::ostream &operator<<(std::ostream &os, const term::Code &code)
{
    if (! term_color_activated) return os;
    return os << "\033[" << (int)code << "m";
}

/**
 * Tuple to string
 */
template <int TDim, typename TType, Tuple_Type TClass>
std::ostream& operator <<(std::ostream& stream, const Tuple<TDim, TType, TClass>& value) {
    stream << term::cyan << Tuple_Type_Names[TClass] << term::reset << "<" << term::dark_gray << TDim << "," << type_name<TType>() << term::reset << ">(";
    stream << term::red;
    for (int i = 0; i < TDim; i++) {
        //stream.width(8);
        stream << +value.at(i);
        if (i < TDim -1) stream << term::reset << ", " << term::red;
    }
    stream << term::reset << ")";
    return stream;
}

/**
 * Matrix to string
 */
template <int TDim, typename TType>
std::ostream& operator <<(std::ostream& stream, const Matrix<TDim, TType>& value) {
    stream << term::cyan << "Matrix" << term::reset <<  "<" << term::dark_gray << TDim << "," << type_name<TType>() << term::reset << ">(\n";
    for (int i = 0; i < TDim; i++) {
        stream << "\t" << value.data[i] <<  (i < TDim -1 ?  ", " : "") << std::endl;
    }
    stream << ")";
    return stream;
}

// ---------------------------------------
// Predefined types
// ---------------------------------------
typedef Tuple<3, double, kPoint>        Point3d;
typedef Tuple<3, double, kDirection>    Direction3d;
typedef Tuple<3, double, kColor>        Color3d;
typedef Tuple<4, double, kColor>        Color4d;
typedef Tuple<2, int, kSize>            Size2i;
typedef Tuple<2, double, kPoint>        Point2d;
typedef Tuple<2, int, kPoint>           Point2i;
typedef Matrix<4, double>             Matrix4d;
typedef Image<3, double>                Image3d;
typedef Image<3, uint8_t>               Image3c;
typedef Image<4, double>                Image4d;
typedef Image<4, uint8_t>               Image4c;
typedef Ray<3, double>                Ray3d;
typedef Triangle<3, double>             Triangle3d;

// ---------------------------------------
} // end of namespace
// ---------------------------------------

#endif // _vgi_tuple_h_

#ifdef TESTING_TUPLE

using namespace vgi;

void test_1()
{
    term_color_activated = false;

    Point3d pt1(4,8,9);
    Point3d pt2(43,7.5,12.7);
    Point3d pt3(4.2);
    Vector3d vec1(10,45,32);
    Direction3d dir1(0,1,0);
    Direction3d dir2(1,0,0);
    Matrix4d mat1;
    mat1.identity(0.0, 2.0);
    Matrix4d mat2;
    mat2.identity(0.0, 3.0);
    Matrix4d mat3 = make_perspective<double>(60.0, 4/3.0, 0.1, 100.0);
    Matrix4d mat4 = make_look_at<double>(Point3d(5, 5, 5), Point3d(0, 0, 0), Direction3d(0, 1,0));

    // std::cout << "dot product of " << pt1 << " and " << pt2 << " is " << dot(pt1, pt2) << std::endl; // no matching function
    // std::cout << "pt1 * vec1 is " << (pt1 * vec1) << std::endl; // no matching function

    std::cout << "dot product of " << dir1 << " and " << dir2 << " is " << term::light_yellow << dot(dir1, dir2) << term::reset << std::endl;
    std::cout << "sum of " << pt1 << " is " << term::light_yellow << sum(pt1) << term::reset << std::endl;
    std::cout << "average of " << pt1 << " is " << term::light_yellow << average(pt1) << term::reset << std::endl;
    std::cout << "pt3 is " << pt3 << std::endl;
    std::cout << "mat1 is " << mat1 << std::endl;
    std::cout << "mat2 is " << mat2 << std::endl;
    std::cout << "mat3 is " << mat3 << std::endl;
    std::cout << "mat4 is " << mat4 << std::endl;
    std::cout << "pt1 * pt2 is " << (pt1 * pt2) << std::endl;
    std::cout << "mat1 * mat2 is " << (mat1 * mat2) << std::endl;
    std::cout << "vec1 is " << vec1 << std::endl;
    std::cout << "mat2 * vec1 is " << (mat2 * vec1) << std::endl;
    std::cout << "mat2 * pt1 is " << (mat2 * pt1) << std::endl;
    std::cout << "vec1 * 3 is " << (vec1 * real(3)) << std::endl;

    Image3d img;
    img.allocate(Size2i(2, 2));
    img.set(Point2i(0, 0), Color3d(1, 0, 0));
    img.set(Point2i(1, 0), Color3d(0, 1, 0));
    img.set(Point2i(0, 1), Color3d(0, 0, 1));
    img.set(Point2i(1, 1), Color3d(1, 0, 1));
    img.set(Point2i(1, 1), Color4d(0, 0, 0, .5));

    Image3c imgbyte;
    convert(img, imgbyte);
    save(imgbyte, "out.ppm");
}

int main()
{
    std::cout.precision(8);
    // std::cout.setf(std::ios::fixed, std:: ios::floatfield);
    test_1();
    return 0;
}

#endif // TESTING_TUPLE