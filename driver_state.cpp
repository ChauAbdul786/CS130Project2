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
    if (face == 1) {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }

    data_geometry one[3];
    data_geometry two[3];

    //Temp Floats
    float Ta = 0;
    float Tb1 = 0;
    float Tb2 = 0;

    vec4 p1;
    vec4 p2;

    vec4 a = (*v0).gl_Position;
    vec4 b = (*v1).gl_Position;
    vec4 c = (*v2).gl_Position;
    
    //Nothing to clip
    if ((a[2] < -a[3]) && (b[2] < -b[3]) && (c[2] < -c[2])) {
        return;
    }

    if ((a[2] < a[3]) && (b[2] >= b[3]) && (c[2] >= c[3])) {
        one[0].data = new float[state.float_per_vertex];
        one[1].data = *in[1];
        one[2].data = *in[2];

        Tb1 = (-b[3] - b[2]) / (a[2] + a[3] - b[3] - b[2]);
        Tb2 = (-a[3] - a[2]) / (c[2] + c[3] - a[3] - a[2]);

        p1 = Tb1 * a + (1 - Tb1) * b;
        p2 = Tb2 * c + (1 - Tb2) * a;

        for (int i = 0; i < state.floats_per_vertex; i++) {
            switch (state.interp_rules[i]) {
            case interp_type::flat:
                one[0].data[i] = v0.data[i];
                break;
            
            case interp_type::smooth:
                one[0].data[i] = Tb2 
            }
        }
    }

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
    vec4 vertices[3];
    vertices[0] = v0.gl_Position / v0.gl_Position[3];
    vertices[1] = v1.gl_Position / v1.gl_Position[3];
    vertices[2] = v2.gl_Position / v2.gl_Position[3];

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
                //state.image_color[i + j * state.image_width] = make_pixel(255, 255, 255); //Simply make white for now 
                data_fragment df;
                df.data = new float[MAX_FLOATS_PER_VERTEX];

                for (int i = 0; i < state.floats_per_vertex; i++) {
                    switch (state.interp_rules[i]) {
                    case(interp_type::invalid):
                        break;
                    case(interp_type::flat):
                        df.data[i] = v0.data[i];
                        break;
                    case(interp_type::smooth):
                        float pointAlph, pointBet, pointGam;
                        pointAlph = ((alpha / v0.gl_Position[3]) / ((alpha / v0.gl_Position[3]) + (beta / v1.gl_Position[3]) + (gamma / v2.gl_Position[3])));
                        pointBet = ((beta / v1.gl_Position[3]) / ((alpha / v0.gl_Position[3]) + (beta / v1.gl_Position[3]) + (gamma / v2.gl_Position[3])));
                        pointGam = ((gamma / v2.gl_Position[3]) / ((alpha / v0.gl_Position[3]) + (beta / v1.gl_Position[3]) + (gamma / v2.gl_Position[3])));

                        df.data[i] = (v0.data[i] * pointAlph) + (v1.data[i] * pointBet) + (v2.data[i] * pointGam);
                        break;
                    case(interp_type::noperspective):
                        df.data[i] = ((v0.data[i] * alpha) + (v1.data[i] * beta) + (v2.data[i] * gamma));
                        break;
                    }
                }

                data_output output;
                state.fragment_shader(df, output, state.uniform_data);

                float z;
                z = (alpha * vertices[0][2]) + (beta * vertices[1][2]) + (gamma * vertices[2][2]);

                if (z < state.image_depth[(state.image_width * j) + (i)]) {
                    state.image_color[(state.image_width * j) + (i)] = make_pixel(255 * output.output_color[0], 255 * output.output_color[1], 255 * output.output_color[2]);
                    state.image_depth[(state.image_width) * j + (i)] = z;
                }

                delete df.data;
            }
        }
    }

    std::cout << "TODO: implement rasterization" << std::endl;
}

