/*
    Defines a couple shapes that are used as objects in a scene.
*/

#pragma once

#include "color.hpp"
#include "ray_math.hpp"
#include "vec3.hpp"

enum Surface_Type {OPAQUE, REFLECTIVE};

struct Entity{
    Color color;
    Surface_Type surface;
    double reflectivity;

    Entity(Color c, Surface_Type st, double r = 0.0) : color{c}, surface{st}, reflectivity{r} {}

    virtual bool RayIntersect(const Ray& ray, Vec3<double>& intersection) = 0;
    virtual Ray NormalAtPoint(const Vec3<double>& point) const = 0;
    virtual Ray ReflectedRay(const Ray& incoming_ray, const Vec3<double>& intersection) const = 0;
};

struct Sphere_Entity : Entity{
    Sphere sphere;

    Sphere_Entity(Sphere s, Color c, Surface_Type st, double r) : Entity{c, st, r}, sphere{s} {}

    bool RayIntersect(const Ray& ray, Vec3<double>& intersection) override{
        return RaySphereIntersection(ray, sphere, intersection);
    }

    Ray NormalAtPoint(const Vec3<double>& point) const override{
        return sphere.NormalAtPoint(point);
    }

    Ray ReflectedRay(const Ray& incoming_ray, const Vec3<double>& intersection) const override{
        return sphere.ReflectedRay(incoming_ray, intersection);
    }
};

struct Rectangle_Entity : Entity{
    Plane plane;  // The plane on which the rectangle exists
    Vec3<double> normal; // Normal vector of the plane
    Vec3<double> A; // "Main" point which meets at the right angle made by B and C
    Vec3<double> B; // Connected to A
    Vec3<double> C; // Connected to A

    Rectangle_Entity(Vec3<double> p1, Vec3<double> p2, Vec3<double> p3, Color c, Surface_Type st, double r) : Entity{c, st, r}, A(p1), B(p2), C(p3) {
        normal = ((B - A).CrossProduct(C - A)).Normalize();
        plane = Plane{A, normal};
    }

    bool RayIntersect(const Ray& ray, Vec3<double>& intersection) override{
        if(RayPlaneIntersection(ray, plane, intersection))
        {
           Vec3<double> AB{A, B};
           Vec3<double> AC{A, C};
           Vec3<double> AI{A, intersection};

           double ABAB_dot = AB.DotProduct(AB);
           double ACAC_dot = AC.DotProduct(AC);
           double AIAB_dot = AI.DotProduct(AB);
           double AIAC_dot = AI.DotProduct(AC);

           return(  AIAB_dot >= 0.0 &&
                    AIAC_dot >= 0.0 &&
                    AIAB_dot <= ABAB_dot &&
                    AIAC_dot <= ACAC_dot);
        }
        return false;
    }

    Ray NormalAtPoint(const Vec3<double>& point) const override{
        return plane.NormalAtPoint(point);
    }

    Ray ReflectedRay(const Ray& incoming_ray, const Vec3<double>& intersection) const override{
        return plane.ReflectedRay(incoming_ray, intersection);
    }

};
