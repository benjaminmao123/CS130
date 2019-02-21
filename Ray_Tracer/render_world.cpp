#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"

#include <climits>
#include <iostream>

extern bool disable_hierarchy;

Render_World::Render_World()
    :background_shader(0),ambient_intensity(0),enable_shadows(true),
    recursion_depth_limit(3)
{}

Render_World::~Render_World()
{
    delete background_shader;
    for(size_t i=0;i<objects.size();i++) delete objects[i];
    for(size_t i=0;i<lights.size();i++) delete lights[i];
}

// Find and return the Hit structure for the closest intersection.  Be careful
// to ensure that hit.dist>=small_t.
Hit Render_World::Closest_Intersection(const Ray& ray)
{
    TODO;
    double min_t = std::numeric_limits<double>::max();
    Hit h_closest = { 0, 0, 0 };

    for (unsigned int i = 0; i < objects.size(); ++i)
    {
        for (int j = 0; j < objects[i]->number_parts; ++j)
        {
            Hit hit = objects[i]->Intersection(ray, j);

            if (hit.dist >= small_t && hit.dist < min_t)
            {
                h_closest = hit;
                min_t = h_closest.dist;
            }
        }
    }
    
    return h_closest;
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
    TODO;
    Ray ray;
    ray.endpoint = camera.position;
    ray.direction = (camera.World_Position(pixel_index) - camera.position).normalized();    
    vec3 color = Cast_Ray(ray, 1);
    camera.Set_Pixel(pixel_index, Pixel_Color(color));
}

void Render_World::Render()
{
    if(!disable_hierarchy)
        Initialize_Hierarchy();

    for(int j=0;j<camera.number_pixels[1];j++)
        for(int i=0;i<camera.number_pixels[0];i++)
            Render_Pixel(ivec2(i,j));
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray, int recursion_depth)
{
    vec3 color;
    TODO; // determine the color here
    Hit h = Closest_Intersection(ray);
    
    if (recursion_depth > recursion_depth_limit)
		color = background_shader->Shade_Surface(ray, ray.endpoint, ray.endpoint, recursion_depth);
	else if (h.object != nullptr)
	{
	    vec3 intersection = ray.Point(h.dist);
	    vec3 normal = h.object->Normal(intersection, h.part);
	    color = h.object->material_shader->Shade_Surface(ray, intersection, normal, recursion_depth);
	}
	else
	    color = background_shader->Shade_Surface(ray, ray.endpoint, ray.endpoint, recursion_depth);
    
    return color;
}

void Render_World::Initialize_Hierarchy()
{
    //TODO Fill in hierarchy.entries; there should be one entry for
    // each part of each object.

    hierarchy.Reorder_Entries();
    hierarchy.Build_Tree();
}
