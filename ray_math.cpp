#include "ray_math.hpp"

#include <algorithm>

#ifndef TOLERANCE
#define DOUBLE_TOLERANCE 10e-7
#endif

Ray& Ray::operator=(const Ray& l)
{
    origin = l.origin;
    direction = l.direction;

    return *this;
}

Ray Sphere::NormalAtPoint(const Vec3<double>& point) const
{
    Ray normal;
    normal.origin = origin;
    normal.direction = (point - origin).Normalize();

    return normal;
}

Ray Sphere::ReflectedRay(const Ray& incoming_ray, const Vec3<double>& intersection) const
{
    Ray reflected_ray;
    reflected_ray.origin = intersection;

    Ray normal = NormalAtPoint(intersection);

    double theta = -1 * normal.direction.DotProduct(incoming_ray.direction);
    reflected_ray.direction = incoming_ray.direction + (2 * theta * normal.direction);
    reflected_ray.direction = reflected_ray.direction.Normalize();

    return reflected_ray;
}

Ray Plane::NormalAtPoint(const Vec3<double>& point) const
{
    return Ray{point, normal};
}

Ray Plane::ReflectedRay(const Ray& incoming_ray, const Vec3<double>& intersection) const
{
    Vec3<double> reflected_ray_direction = incoming_ray.direction - 2 * (incoming_ray.direction.DotProduct(normal)) * normal;
    reflected_ray_direction = reflected_ray_direction.Normalize();
    return Ray{intersection, reflected_ray_direction};
}


bool RaySphereIntersection(const Ray& r, const Sphere& s, Vec3<double>& intersection)
{
    double v = r.direction.DotProduct(s.origin - r.origin);
    double discriminant = std::pow(s.radius, 2.0) - ((s.origin - r.origin).DotProduct(s.origin - r.origin) - std::pow(v, 2.0));

    // There are no intersections
    if(discriminant < 0.0)
        return false;
    else
    {
        double distance = std::sqrt(discriminant);
        double t_param = v - distance;

        // If the t_param is negative, the intersection is "behind" the ray. For our purpose we throw this out and return false
        if(t_param < 0.0)
            return false;

        // An intersection exists and it is "in front" of the ray
        intersection = r.origin + ((t_param) * r.direction);
        return true;
    }
}

bool RayPlaneIntersection(const Ray& r, const Plane& p, Vec3<double>& intersection)
{
    double dividend = (p.origin - r.origin).DotProduct(p.normal);
    double divisor = r.direction.DotProduct(p.normal);

    // The line and the plane are parallel
    if(fabs(divisor) < DOUBLE_TOLERANCE)
    {
        return false;
    }
    else
    {
        double distance = dividend / divisor;
        if(distance < 0.0)
            return false;
        intersection = r.origin + (distance * r.direction);
        return true;
    }
}