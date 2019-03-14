#ifndef APPLICATION_H
#define APPLICATION_H

#include <vector>
#include <map>

#include "gl_texture.h"
#include "gl_viewer.h"
#include "obj.h"
#include "timer.h"

using namespace std;

typedef std::map<std::string, gl_image_texture*>
    gl_image_texture_map;

struct Particle
{
	Particle(const vec3 &Position, const vec3 &Velocity, const vec3 &Color, const float &Mass)
		: Position(Position), Velocity(Velocity), Color(Color), Mass(Mass), AppliedForces(vec3(0, -Mass * 9.8, 0)), Duration(0)
	{
		
	}
	
	void Euler_Step(const float &h)
	{
		Position = Position + h * Velocity;
		Velocity = Velocity + (h / Mass) * AppliedForces;
	}
	
	void Reset_Forces()
	{
		AppliedForces = vec3(0, 0, 0);
	}
	
	void Handle_Collision(const float &damping, const float &coeff_restitution)
	{
		if (Position[1] < 0)
		{
			Position[1] = 0;
			Velocity[1] = -coeff_restitution * Velocity[1];
			Velocity[0] = damping * Velocity[0];
			Velocity[2] = damping * Velocity[2];
		}
	}
	
	vec3 Position;
	vec3 Velocity;
	vec3 Color;
	float Mass;
	vec3 AppliedForces;
	float Duration;
};

extern vector<Particle> ParticleList;

class application : public gl_viewer
{
public:
    application();
    ~application();
    void init_event();
    void draw_event();
    void mouse_click_event(int button, int button_state, int x, int y);
    void mouse_move_event(int x, int y);
    void keyboard_event(unsigned char key, int x, int y);

private:
    bool raytrace;
    int rendmode;
    int npart;
    timer t;
    obj o;
    gl_image_texture_map texs;
    bool paused;
    float sim_t;
    bool draw_volcano;
    float h;
};

#endif
