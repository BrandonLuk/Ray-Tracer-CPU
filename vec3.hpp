#pragma once

#include <cmath>
#include <ostream>

#ifndef TOLERANCE
#define DOUBLE_TOLERANCE 10e-7
#endif

template <typename T>
struct Vec3{

public:
    T x, y, z;

    Vec3() = default;
    Vec3(T x, T y, T z) : x(x), y(y), z(z) {}
    Vec3(const Vec3<T>& other){
        x = other.x;
        y = other.y;
        z = other.z;
    }
    Vec3(const Vec3<T>& left, const Vec3<T>& right){
        x = right.x - left.x;
        y = right.y - left.y;
        z = right.z - left.z;
    }
    ~Vec3() = default;

    Vec3& operator=(const Vec3& v);
    Vec3 operator+(const Vec3& v) const;
    Vec3 operator-(const Vec3& v) const;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const Vec3<U>& v);

    double Magnitude();
    Vec3<T> Normalize();
    Vec3<T> CrossProduct(const Vec3<T>& v);
    double DotProduct(const Vec3<T>& v) const;
    double Distance(const Vec3<T>& v);

    Vec3<T> Rotate_x(double theta);
    Vec3<T> Rotate_y(double theta);
};

template <typename T>
Vec3<T> operator*(double d, const Vec3<T>& v);

template <typename T>
Vec3<T> operator*(const Vec3<T>& v, double d);












template <typename T>
Vec3<T>& Vec3<T>::operator=(const Vec3<T>& v)
{
    x = v.x;
    y = v.y;
    z = v.z;
    return *this;
}

template <typename T>
Vec3<T> Vec3<T>::operator+(const Vec3& v) const
{
    Vec3<T> result;

    result.x = x + v.x;
    result.y = y + v.y;
    result.z = z + v.z;

    return result;
}

template <typename T>
Vec3<T> Vec3<T>::operator-(const Vec3& v) const
{
    Vec3<T> result;

    result.x = x - v.x;
    result.y = y - v.y;
    result.z = z - v.z;

    return result;
}

template <typename T>
Vec3<T> operator*(double d, const Vec3<T>& v)
{
    Vec3<T> result;

    result.x = d * v.x;
    result.y = d * v.y;
    result.z = d * v.z;

    return result;
}

template <typename T>
Vec3<T> operator*(const Vec3<T>& v, double d)
{
    return d * v;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vec3<T>& v)
{
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";

    return os;
}

template <typename T>
double Vec3<T>::Magnitude()
{
    return std::sqrt(x * x + y * y + z * z);
}

template <typename T>
Vec3<T> Vec3<T>::Normalize()
{
    Vec3<T> v;
    double magnitude = Magnitude();

    v.x = x / magnitude;
    v.y = y / magnitude;
    v.z = z / magnitude;

    return v;
}

template <typename T>
Vec3<T> Vec3<T>::CrossProduct(const Vec3<T>& v)
{
    Vec3<T> result;

    result.x = (y * v.z) - (z * v.y);
    result.y = (z * v.x) - (x * v.z);
    result.z = (x * v.y) - (y * v.x);

    return result;
}

template <typename T>
double Vec3<T>::DotProduct(const Vec3<T>& v) const
{
    return (x * v.x) + (y * v.y) + (z * v.z);
}

template <typename T>
double Vec3<T>::Distance(const Vec3<T>& v)
{
    return std::sqrt(std::pow(v.x - x, 2) + std::pow(v.y - y, 2) + std::pow(v.z - z, 2));
}

template <typename T>
Vec3<T> Vec3<T>::Rotate_x(double theta)
{
    Vec3<double> rotated;

    rotated.x = x;
    rotated.y = (std::cos(theta) * y) + (-1 * std::sin(theta) * z);
    rotated.z = (std::sin(theta) * y) + (std::cos(theta) * z);

    return rotated;
}

template <typename T>
Vec3<T> Vec3<T>::Rotate_y(double theta)
{
    Vec3<double> rotated;

    rotated.x = (std::cos(theta) * x) + (std::sin(theta) * z);
    rotated.y = y;
    rotated.z = (-1 * std::sin(theta) * x) + (std::cos(theta) * z);

    return rotated;
}