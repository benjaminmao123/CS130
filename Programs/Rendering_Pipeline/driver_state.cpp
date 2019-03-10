#include "driver_state.h"
#include <cstring>
#include <limits>
#include <vector>
#include <cmath>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete[] image_color;
    delete[] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state &state, int width, int height)
{
    state.image_width = width;
    state.image_height = height;
    state.image_color = 0;
    state.image_depth = 0;
    std::cout << "TODO: allocate and initialize state.image_color and state.image_depth." << std::endl;

    state.image_color = new pixel[width * height];

    for (int i = 0; i < width * height; ++i)
        state.image_color[i] = make_pixel(0, 0, 0);

    state.image_depth = new float[width * height];

    for (int i = 0; i < width * height; ++i)
        state.image_depth[i] = 1;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state &state, render_type type)
{
    // std::cout << "TODO: implement rendering." << std::endl;
    switch (type)
    {
    case render_type::triangle:
    {
        const data_geometry *d_GeometryArr[3];
        data_geometry g[3];
        data_vertex v[3];
        int num_triangles = state.num_vertices / 3;

        int k = 0;

        for (int i = 0; i < num_triangles; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                v[j].data = &state.vertex_data[k];
                g[j].data = v[j].data;
                state.vertex_shader(v[j], g[j], state.uniform_data);
                d_GeometryArr[j] = &g[j];
                k += state.floats_per_vertex;
            }

            clip_triangle(state, d_GeometryArr, 0);
        }
        break;
    }
    case render_type::indexed:
    {
        const data_geometry *d_GeometryArr[3];
        data_geometry g[3];
        data_vertex v[3];

        for (int i = 0; i < state.num_triangles * 3; i += 3)
        {
            for (int j = 0; j < 3; ++j)
            {
                v[j].data = &state.vertex_data[state.index_data[i + j] * state.floats_per_vertex];
                g[j].data = v[j].data;
                state.vertex_shader(v[j], g[j], state.uniform_data);
                d_GeometryArr[j] = &g[j];
            }

            clip_triangle(state, d_GeometryArr, 0);
        }    
        break;
    }
    case render_type::fan:
    {
        const data_geometry *d_GeometryArr[3];
        data_geometry g[3];
        data_vertex v[3];

        for (int i = 0; i < state.num_vertices; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {   
                if (j) v[j].data = state.vertex_data + ((i + j) * state.floats_per_vertex);
                else v[j].data = state.vertex_data + ((0 + j) * state.floats_per_vertex);
                g[j].data = v[j].data;
                state.vertex_shader(v[j], g[j], state.uniform_data);
                d_GeometryArr[j] = &g[j];
            }

            clip_triangle(state, d_GeometryArr, 0); 
        }
        break; 
    }
    case render_type::strip:
    {
        const data_geometry *d_GeometryArr[3];
        data_geometry g[3];
        data_vertex v[3];

        for (int i = 0; i < state.num_vertices - 2; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                v[j].data = &state.vertex_data[(i + j) * state.floats_per_vertex];
                g[j].data = v[j].data;
                state.vertex_shader(v[j], g[j], state.uniform_data);
                d_GeometryArr[j] = &g[j];
            }

            clip_triangle(state, d_GeometryArr, 0);
        }
        break;               
    }
    default:
        break;
    }
}

void Interpolate(driver_state &state, const Edge & edge, data_geometry & new_data, const int & Index, const data_geometry *in[3])
{
    for (int i = 0; i < state.floats_per_vertex; ++i)
    {
        switch (state.interp_rules[i])
        {
        case interp_type::flat:
            new_data.data[i] = in[Index]->data[i];
            break;
        case interp_type::smooth:
            new_data.data[i] = edge.beta * edge.start_data.data[i] + (1 - edge.beta) * edge.end_data.data[i];
            break;
        case interp_type::noperspective:
        {
            float alpha = edge.beta * edge.start_data.gl_Position[3] / (edge.beta * edge.start_data.gl_Position[3] + (1 - edge.beta) * edge.end_data.gl_Position[3]);
            new_data.data[i] = alpha * edge.start_data.data[i] + (1 - alpha) * edge.end_data.data[i];
            break;
        }
        default:
            break;
        }
    }
}

vector<vector<data_geometry>> Generate_Triangles(driver_state &state, vector<Edge> & edges, const data_geometry *in[3], const int & Index, const int & Case)
{
    vec4 new_p0 = edges[0].beta * edges[0].start_vertex + (1 - edges[0].beta) * edges[0].end_vertex;
    vec4 new_p1 = edges[1].beta * edges[1].start_vertex + (1 - edges[1].beta) * edges[1].end_vertex;

    data_geometry new_v0[3];
    data_geometry new_v1[3];
    vector<vector<data_geometry>> new_triangles;

    switch (Case)
    {
    case 0:
    {
        new_v0[1].data = new float[state.floats_per_vertex];
        new_v0[2].data = new float[state.floats_per_vertex];
        Interpolate(state, edges[0], new_v0[1], 1, in);
        Interpolate(state, edges[1], new_v0[2], 2, in);
        new_v0[1].gl_Position = new_p0;
        new_v0[2].gl_Position = new_p1;
        vector<data_geometry> tri0 = {*in[0], new_v0[1], new_v0[2]};
        new_triangles.emplace_back(tri0);          
        break;
    }
    case 1:
    {
        new_v0[0].data = new float[state.floats_per_vertex];
        new_v0[2].data = new float[state.floats_per_vertex];
        Interpolate(state, edges[0], new_v0[0], 0, in);
        Interpolate(state, edges[1], new_v0[2], 2, in);
        new_v0[0].gl_Position = new_p0;
        new_v0[2].gl_Position = new_p1;
        vector<data_geometry> tri0 = {new_v0[0], *in[1], new_v0[2]};
        new_triangles.emplace_back(tri0);   
        break;
    }
    case 2:
    {
        new_v0[0].data = new float[state.floats_per_vertex];
        new_v0[1].data = new float[state.floats_per_vertex];
        Interpolate(state, edges[1], new_v0[0], 0, in);
        Interpolate(state, edges[0], new_v0[1], 1, in);
        new_v0[0].gl_Position = new_p1;
        new_v0[1].gl_Position = new_p0;
        vector<data_geometry> tri0 = {new_v0[0], new_v0[1], *in[2]};
        new_triangles.emplace_back(tri0);              
        break;
    }
    case 3:
    {
        new_v0[2].data = new float[state.floats_per_vertex];
        Interpolate(state, edges[1], new_v0[2], 2, in);
        new_v0[2].gl_Position = new_p1;
        vector<data_geometry> tri0 = {*in[0], *in[1], new_v0[2]};
        new_triangles.emplace_back(tri0);     

        new_v1[2].data = new float[state.floats_per_vertex];
        Interpolate(state, edges[0], new_v1[2], 2, in);
        new_v1[2].gl_Position = new_p0;
        vector<data_geometry> tri1 = {new_v0[2], *in[1], new_v1[2]};
        new_triangles.emplace_back(tri1); 
        break;
    }
    case 4:
    {
        new_v0[1].data = new float[state.floats_per_vertex];
        Interpolate(state, edges[0], new_v0[1], 1, in);
        new_v0[1].gl_Position = new_p0;
        vector<data_geometry> tri0 = {*in[0], new_v0[1], *in[2]};
        new_triangles.emplace_back(tri0);     

        new_v1[1].data = new float[state.floats_per_vertex];
        Interpolate(state, edges[1], new_v1[1], 1, in);
        new_v1[1].gl_Position = new_p1;
        vector<data_geometry> tri1 = {new_v0[1], new_v1[1], *in[2]};
        new_triangles.emplace_back(tri1);          
        break;
    }
    case 5:
    {
        new_v0[0].data = new float[state.floats_per_vertex];
        Interpolate(state, edges[1], new_v0[0], 0, in);
        new_v0[0].gl_Position = new_p1;
        vector<data_geometry> tri0 = {new_v0[0], *in[1], *in[2]};
        new_triangles.emplace_back(tri0);     

        new_v1[0].data = new float[state.floats_per_vertex];
        Interpolate(state, edges[0], new_v1[0], 0, in);
        new_v1[0].gl_Position = new_p0;
        vector<data_geometry> tri1 = {new_v1[0], *in[1], new_v0[0]};
        new_triangles.emplace_back(tri1);   
        break; 
    }   
    default:
        break;   
    }

    return new_triangles;
}

vector<vector<data_geometry>> Cut_Triangle(driver_state &state, const int & Index, const data_geometry *in[3])
{
    vec4 a = in[0]->gl_Position;
    vec4 b = in[1]->gl_Position;
    vec4 c = in[2]->gl_Position;
    vector<vector<data_geometry>> new_triangles = {{*in[0], *in[1], *in[2]}};

    if (a[Index] >= -a[3] && b[Index] < -b[3] && c[Index] < -c[3])
    {
        Edge e0(a, b, *in[0], *in[1], Index);
        Edge e1(c, a, *in[2], *in[0], Index);
        vector<Edge> edges = {e0, e1};
        new_triangles = Generate_Triangles(state, edges, in, Index, 0);        
    }
    else if (a[Index] < -a[3] && b[Index] >= -b[3] && c[Index] < -c[3])
    {
        Edge e0(a, b, *in[0], *in[1], Index);
        Edge e1(b, c, *in[1], *in[2], Index);
        vector<Edge> edges = {e0, e1};
        new_triangles = Generate_Triangles(state, edges, in, Index, 1);          
    }
    else if (a[Index] < -a[3] && b[Index] < -b[3] && c[Index] >= -c[3])
    {
        Edge e0(b, c, *in[1], *in[2], Index);
        Edge e1(c, a, *in[2], *in[0], Index);
        vector<Edge> edges = {e0, e1};
        new_triangles = Generate_Triangles(state, edges, in, Index, 2);          
    }
    else if (a[Index] >= -a[3] && b[Index] >= -b[3] && c[Index] < -c[3])
    {
        Edge e0(b, c, *in[1], *in[2], Index);
        Edge e1(c, a, *in[2], *in[0], Index);
        vector<Edge> edges = {e0, e1};
        new_triangles = Generate_Triangles(state, edges, in, Index, 3);         
    } 
    else if (a[Index] >= -a[3] && b[Index] < -b[3] && c[Index] >= -c[3])
    {
        Edge e0(a, b, *in[0], *in[1], Index);
        Edge e1(b, c, *in[1], *in[2], Index);
        vector<Edge> edges = {e0, e1};
        new_triangles = Generate_Triangles(state, edges, in, Index, 4);         
    }           
    else if (a[Index] < -a[3] && b[Index] >= -b[3] && c[Index] >= -c[3])
    {
        Edge e0(a, b, *in[0], *in[1], Index);
        Edge e1(c, a, *in[2], *in[0], Index);
        vector<Edge> edges = {e0, e1};
        new_triangles = Generate_Triangles(state, edges, in, Index, 5);
    }

    return new_triangles;
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state &state, const data_geometry *in[3], int face)
{
    if (face == 3)
    {
        rasterize_triangle(state, in);
        return;
    }
    // std::cout << "TODO: implement clipping. (The current code passes the triangle through without clipping them.)" << std::endl;
    vec4 a = in[0]->gl_Position;
    vec4 b = in[1]->gl_Position;
    vec4 c = in[2]->gl_Position;    
    const data_geometry *input[3];
    vector<vector<data_geometry>> new_vertices;

    switch (face)
	{
    case 0:
        if (a[0] >= -a[3] && b[0] >= -b[3] && c[0] >= -c[3])
            clip_triangle(state, in, face + 1);
        else if (a[0] < -a[3] && b[0] < -b[3] && c[0] < -c[3])
            return;
        else
        {
            new_vertices = Cut_Triangle(state, 0, in);
            for (unsigned int i = 0; i < new_vertices.size(); ++i)
            {
                for (unsigned int j = 0; j < new_vertices[i].size(); ++j)
                    input[j] = &new_vertices[i][j];

                clip_triangle(state, input, face + 1);
            }            
        }
        break;
    case 1:
        if (a[1] >= -a[3] && b[1] >= -b[3] && c[1] >= -c[3])
            clip_triangle(state, in, face + 1);
        else if (a[1] < -a[3] && b[1] < -b[3] && c[1] < -c[3])
            return;
        else
        {
            new_vertices = Cut_Triangle(state, 1, in);
            for (unsigned int i = 0; i < new_vertices.size(); ++i)
            {
                for (unsigned int j = 0; j < new_vertices[i].size(); ++j)
                    input[j] = &new_vertices[i][j];

                clip_triangle(state, input, face + 1);
            }            
        }   
        break;
    case 2:
    {
        if (a[2] >= -a[3] && b[2] >= -b[3] && c[2] >= -c[3])
            clip_triangle(state, in, face + 1);
        else if (a[2] < -a[3] && b[2] < -b[3] && c[2] < -c[3])
            return;
        else
        {
            new_vertices = Cut_Triangle(state, 2, in);
            for (unsigned int i = 0; i < new_vertices.size(); ++i)
            {
                for (unsigned int j = 0; j < new_vertices[i].size(); ++j)
                    input[j] = &new_vertices[i][j];

                clip_triangle(state, input, face + 1);
            }
        }
        break;
    }
	default:
		break;
	}
}

void Interpolate(driver_state &state, const vec3 &bary_prime, const vec3 &w, data_fragment &d_frag, const data_geometry *in[3])
{
    for (int i = 0; i < state.floats_per_vertex; ++i)
        switch (state.interp_rules[i])
        {
        case interp_type::flat:
            d_frag.data[i] = in[0]->data[i];
            break;
        case interp_type::smooth:
        {
            float c = (bary_prime[0] / w[0]) + (bary_prime[1] / w[1]) + (bary_prime[2] / w[2]);
            vec3 bary = {(bary_prime[0] / w[0]) / c, (bary_prime[1] / w[1]) / c, (bary_prime[2] / w[2]) / c};
            d_frag.data[i] = (bary[0] * in[0]->data[i]) + (bary[1] * in[1]->data[i]) + (bary[2] * in[2]->data[i]);
            break;
        }
        case interp_type::noperspective:
            d_frag.data[i] = bary_prime[0] * in[0]->data[i] + bary_prime[1] * in[1]->data[i] + bary_prime[2] * in[2]->data[i];
            break;
        default:
            break;
        }
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state &state, const data_geometry *in[3])
{
    //std::cout << "TODO: implement rasterization" << std::endl;
    vec3 w = {in[0]->gl_Position[3], in[1]->gl_Position[3], in[2]->gl_Position[3]};
    vec3 x = {in[0]->gl_Position[0] / w[0], in[1]->gl_Position[0] / w[1], in[2]->gl_Position[0] / w[2]};
    vec3 y = {in[0]->gl_Position[1] / w[0], in[1]->gl_Position[1] / w[1], in[2]->gl_Position[1] / w[2]};
    vec3 z = {in[0]->gl_Position[2] / w[0], in[1]->gl_Position[2] / w[1], in[2]->gl_Position[2] / w[2]};

    vec2 v0 = {(float)((state.image_width / 2) * x[0]) + (float)((state.image_width / 2) - 0.5), 
               (float)((state.image_height / 2) * y[0]) + (float)((state.image_height / 2) - 0.5)};
    vec2 v1 = {(float)((state.image_width / 2) * x[1]) + (float)((state.image_width / 2) - 0.5), 
               (float)((state.image_height / 2) * y[1]) + (float)((state.image_height / 2) - 0.5)};
    vec2 v2 = {(float)((state.image_width / 2) * x[2]) + (float)((state.image_width / 2) - 0.5), 
               (float)((state.image_height / 2) * y[2]) + (float)((state.image_height / 2) - 0.5)};

    vec3 area = {(v2[0] - v0[0]) * (v1[1] - v0[1]) - (v2[1] - v0[1]) * (v1[0] - v0[0]), 0, 0};

    vec2 min_index = {min(v0[0], min(v1[0], v2[0])), min(v0[1], min(v1[1], v2[1]))};
    vec2 max_index = {max(v0[0], max(v1[0], v2[0])), max(v0[1], max(v1[1], v2[1]))};

    if (min_index[0] < 0) min_index[0] = 0;
    if (min_index[1] < 0) min_index[1] = 0;
    if (max_index[0] > state.image_width) max_index[0] = state.image_width;
    if (max_index[1] > state.image_height) max_index[1] = state.image_height;

    for (int i = min_index[0]; i < max_index[0]; ++i)
        for (int j = min_index[1]; j < max_index[1]; ++j)
        {
            vec2 p = {(float)i, (float)j};
            area[1] = (p[0] - v1[0]) * (v2[1] - v1[1]) - (p[1] - v1[1]) * (v2[0] - v1[0]);
            area[2] = (p[0] - v2[0]) * (v0[1] - v2[1]) - (p[1] - v2[1]) * (v0[0] - v2[0]);
			vec3 bary_prime = {area[1] / area[0], area[2] / area[0], 0};
            bary_prime[2] = 1 - bary_prime[0] - bary_prime[1];

            if (bary_prime[0] >= 0 && bary_prime[1] >= 0 && bary_prime[2] >= 0)
            {
                float zDepth = (bary_prime[0] * z[0]) + (bary_prime[1] * z[1]) + (bary_prime[2] * z[2]);

                if (state.image_depth[i + j * state.image_width] > zDepth)
                {
                    data_fragment d_frag;
                    d_frag.data = new float[state.floats_per_vertex];
                    data_output d_output;

                    Interpolate(state, bary_prime, w, d_frag, in);

                    state.fragment_shader(d_frag, d_output, state.uniform_data);

                    float r = d_output.output_color[0] * 255;
                    float g = d_output.output_color[1] * 255;
                    float b = d_output.output_color[2] * 255;

                    state.image_color[i + j * state.image_width] = make_pixel(r, g, b);
                    state.image_depth[i + j * state.image_width] = zDepth;
                }
            }
        }
}
