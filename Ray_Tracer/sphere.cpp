#include "sphere.h"
#include "ray.h"

// Determine if the ray intersects with the sphere
Hit Sphere::Intersection(const Ray& ray, int part) const
{
    TODO;
    double discriminant = pow(dot(ray.direction, (ray.endpoint - center)), 2) - ((ray.endpoint - center).magnitude_squared() - pow(radius, 2));

    if (discriminant > 0)
    {
        double d = -(dot(ray.direction, (ray.endpoint - center))) - sqrt(discriminant);

        if (d >= small_t)
            return { this, d, part };
    }
    
    return { 0, 0, 0 };
}

vec3 Sphere::Normal(const vec3& point, int part) const
{
    vec3 normal;
    TODO; // compute the normal direction
    normal = point - center;

    return normal.normalized();
}

Box Sphere::Bounding_Box(int part) const
{
    Box box;
    TODO; // calculate bounding box
    
    
    return box;
}
