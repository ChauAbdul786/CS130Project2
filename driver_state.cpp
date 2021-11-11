#include "driver_state.h"
#include <cstring>
#include <limits>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];

    for (int i = 0; i < width * height; i++) {
        state.image_depth[i] = 1.0; //For z-buffering
        state.image_color[i] = make_pixel(0, 0, 0); //Make pixel black 
    }
    

}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    data_geometry* geo = new data_geometry[MAX_FLOATS_PER_VERTEX];

    for (int i = 0; i < state.num_vertices; i++) {
        geo[i].data = new float[MAX_FLOATS_PER_VERTEX];
        data_vertex data;
        data.data = (&state.vertex_data[state.floats_per_vertex * i]);
        state.vertex_shader(data, geo[i], state.uniform_data);
    }

    switch (type) {
    case render_type::triangle:
        for (int i = 0; i < state.num_vertices - 2; i += 3) {
            rasterize_triangle(state, geo[i], geo[i + 1], geo[i + 2]);
        }
        break;
    case render_type::indexed:
        break; //TODO
    case render_type::fan:
        break; //TODO
    case render_type::strip:
        break; //TODO
    default:
        break;
    }

    for (int i = 0; i < state.num_vertices; i++) {
        delete[] geo[i].data;
    }
    delete[] geo;

    std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,v0,v1,v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    vec3 a = vec3(v0.gl_Position[0] / v0.gl_Position[3], v0.gl_Position[1] / v0.gl_Position[3], v0.gl_Position[2] / v0.gl_Position[3]);
    a[0] = (state.image_width - 1) * ((a[0] + 1) / 2);
    a[1] = (state.image_height - 1) * ((a[1] + 1) / 2);

    vec3 b = vec3(v1.gl_Position[0] / v1.gl_Position[3], v1.gl_Position[1] / v1.gl_Position[3], v1.gl_Position[2] / v1.gl_Position[3]);
    b[0] = (state.image_width - 1) * ((b[0] + 1) / 2);
    b[1] = (state.image_height - 1) * ((b[1] + 1) / 2);

    vec3 c = vec3(v2.gl_Position[0] / v2.gl_Position[3], v2.gl_Position[1] / v2.gl_Position[3], v2.gl_Position[2] / v2.gl_Position[3]);
    c[0] = (state.image_width - 1) * ((c[0] + 1) / 2);
    c[1] = (state.image_height - 1) * ((c[1] + 1) / 2);

    float area = (.5 * ((b[0] * c[1]) - (c[0] * b[1]) + (c[0] * a[1]) - (a[0] * c[1]) + (a[0] * b[1]) - (b[0] * a[1])));

    float alpha, beta, gamma;
    for (int i = 0; i < state.image_width; i++) {
        for (int j = 0; j < state.image_height; j++) {
            vec3 p = vec3(i, j, 0);
            alpha = (0.5 * ((b[0] * c[1] - c[0] * b[1]) + (c[0] * p[1] - p[0] * c[1]) + (p[0] * b[1] - b[0] * p[1])) / area);
            beta = (0.5 * ((p[0] * c[1] - c[0] * p[1]) + (c[0] * a[1] - a[0] * c[1]) + (a[0] * p[1] - p[0] * a[1])) / area);
            gamma = (0.5 * ((b[0] * p[1] - p[0] * b[1]) + (p[0] * a[1] - a[0] * p[1]) + (a[0] * b[1] - b[0] * a[1])) / area);
        
            if (((alpha >= 0) && (beta >= 0)) && (gamma >= 0)) {
                state.image_color[i + j * state.image_width] = make_pixel(255, 255, 255); //Simply make white for now 
            }
        }
    }

    std::cout<<"TODO: implement rasterization"<<std::endl;
}
