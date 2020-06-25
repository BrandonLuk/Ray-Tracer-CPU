#pragma once

#include "vec3.hpp"

#include <cmath>

#ifndef TOLERANCE
#define DOUBLE_TOLERANCE 10e-7
#endif

struct Ray{
    Vec3<double> origin;
    Vec3<double> direction;

    Ray() = default;
    Ray(Vec3<double> o, Vec3<double> d) : origin(o), direction(d) {}
    Ray(const Ray& other){
        origin = other.origin;
        direction = other.direction;
    }

    Ray& operator=(const Ray& l);
};

struct Sphere{
    Vec3<double> origin;
    double radius;

    Ray NormalAtPoint(const Vec3<double>& point) const;
    Ray ReflectedRay(const Ray& incoming_ray, const Vec3<double>& intersection) const;
};

struct Plane{
    Vec3<double> origin;
    Vec3<double> normal;

    Ray NormalAtPoint(const Vec3<double>& point) const;
    Ray ReflectedRay(const Ray& incoming_ray, const Vec3<double>& intersection) const;
};

bool RaySphereIntersection(const Ray& l, const Sphere& s, Vec3<double>& intersection);
bool RayPlaneIntersection(const Ray& l, const Plane& p, Vec3<double>& intersection);

