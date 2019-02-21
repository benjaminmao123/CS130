#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"

vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
    color = color_ambient * world.ambient_color * world.ambient_intensity;

    TODO; //determine the color
    for (auto i : world.lights)
    {
        Ray l_ray;
        l_ray.direction = i->position - intersection_point;

        if (world.enable_shadows) 
        {
			Ray s_ray;
            s_ray.endpoint = intersection_point;
            s_ray.direction = i->position - intersection_point;
            double s_mag = s_ray.direction.magnitude();
            s_ray.direction = s_ray.direction.normalized();

            Hit hit = world.Closest_Intersection(s_ray);
            vec3 col = s_ray.direction * hit.dist;

			if (!(hit.object != 0 && s_mag > col.magnitude()))
            {
                vec3 d_intensity = i->Emitted_Light(l_ray.direction);
                double d_clamp = std::max(0.0, dot(l_ray.direction.normalized(), normal));
                color += color_diffuse * d_intensity * d_clamp;

                vec3 s_intensity = (2.0 * (dot(l_ray.direction.normalized(), normal) * normal) - l_ray.direction.normalized());
                double s_clamp = pow(std::max(0.0, dot(s_intensity, -ray.direction.normalized())), specular_power);
                color += color_specular * d_intensity * s_clamp;
			}
        }
        else
        {
            vec3 d_intensity = i->Emitted_Light(l_ray.direction);
            double d_clamp = std::max(0.0, dot(l_ray.direction.normalized(), normal));
            color += color_diffuse * d_intensity * d_clamp;

            vec3 s_intensity = (2.0 * (dot(l_ray.direction.normalized(), normal) * normal) - l_ray.direction.normalized());
            double s_clamp = pow(std::max(0.0, dot(s_intensity, -ray.direction.normalized())), specular_power);
            color += color_specular * d_intensity * s_clamp;
        }
    }

    return color;
}