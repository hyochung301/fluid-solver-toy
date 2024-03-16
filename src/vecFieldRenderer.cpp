#include "vecFieldRenderer.h"
#include <flgl/tools.h>
#include <random>
#include <flgl/logger.h>
#include <cmath>
LOG_MODULE(vec_renderer)

/*
C Specification
void glLineWidth(	GLfloat width);
 
Parameters
width
Specifies the width of rasterized lines. The initial value is 1.

Description
glLineWidth specifies the rasterized width of both aliased and antialiased lines. 
Using a line width other than 1 has different effects, depending on whether line antialiasing is enabled. 
To enable and disable line antialiasing, call glEnable and glDisable with argument GL_LINE_SMOOTH. 
Line antialiasing is initially disabled.

If line antialiasing is disabled, the actual width is determined by rounding the supplied width to the nearest integer. 
(If the rounding results in the value 0, it is as if the line width were 1.) If ∣Δx∣>=∣Δy∣, 
i pixels are filled in each column that is rasterized, where i is the rounded value of width. 
Otherwise, i pixels are filled in each row that is rasterized.

If antialiasing is enabled, line rasterization produces a fragment for each pixel square that intersects the region 
lying within the rectangle having width equal to the current line width, length equal to the actual length of the line, 
and centered on the mathematical line segment. The coverage value for each fragment is the window coordinate area 
of the intersection of the rectangular region with the corresponding pixel square. This value is saved and used
in the final rasterization step.

Not all widths can be supported when line antialiasing is enabled.
If an unsupported width is requested, the nearest supported width is used. 
Only width 1 is guaranteed to be supported; others depend on the implementation. 
Likewise, there is a range for aliased line widths as well. To query the range of supported widths 
and the size difference between supported widths within the range, call glGet with arguments 
GL_ALIASED_LINE_WIDTH_RANGE, GL_SMOOTH_LINE_WIDTH_RANGE, and GL_SMOOTH_LINE_WIDTH_GRANULARITY.
*/

template <>
void VertexBuffer<vecfield_vert>::attach_to_vao(VertexArray const& vao) const {
	vao.bind();
	this->bind();
	vao.attrib(0,							// layout index
			   2, GL_FLOAT,					// dimension and repr
			   sizeof(vecfield_vert)/2, 				// size of vertex
			   0);		// offset of member

	// vao.attrib(0,							// layout index
	// 		   2, GL_FLOAT,					// dimension and repr
	// 		   sizeof(vecfield_vert), 				// size of vertex
	// 		   offsetof(vecfield_vert,pos));		// offset of member
	// vao.attrib(1,							// layout index
	// 		   2, GL_FLOAT,					// dimension and repr
	// 		   sizeof(vecfield_vert), 				// size of vertex
	// 		   offsetof(vecfield_vert,vec));		// offset of member
}

vecFieldRenderer::vecFieldRenderer(Field const& f) : FieldRenderer(f) {}
void vecFieldRenderer::init() {
	int n = field.solver.dim();
	buffer = new vecfield_vert[n*n];
	line_shad = Shader::from_source("vector_vert", "color");
	mesh.create();
}
#define max(a,b) (a>b?a:b)
void vecFieldRenderer::prepare() {
	#define COEFF_SCALE (80.0f)
	const float size = 2.0f / (n - 1); // Calculate the size of each step in the grid
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.8, 1.2);
	const float coeff = dis(gen) * (COEFF_SCALE / (n - 1));
    std::vector<uint32_t> indices;
    float* const u = field.solver.x_buffer();
    float* const v = field.solver.y_buffer();
    // Fill vertices buffer
    size_t bufi = 0;
    for (size_t j = 0; j < (n); j+=4) {
        for (size_t i = 0; i < (n); i+=4) {

        	glm::vec2 vec = glm::vec2(u[i + j * n], v[i + j * n]);
        	float l = glm::length(vec);
        	vec *= (1./l);
        	l = sqrt(l+1) - 1;
        	l *= 2;
        	l = max(l-0.1,0.);
        	vec *= l;
            if (l<=0.) continue;

            // Scale positions to [-1, 1]
            buffer[bufi].pos = glm::vec2(-1.0f + size * i, -1.0f + size * j);

            // Calculate and scale vector
            buffer[bufi].vec = buffer[bufi].pos + coeff * vec; // Adjust vector by the given coefficient and add to position

            // Add line indices (each vertex to its corresponding vector tip)
            size_t index = 2*bufi;
            indices.push_back(index); // start of line (position)
            indices.push_back(index+1); // end of line (position + vector)
            bufi++;
        }
    }

    mesh.vao.bind();
	mesh.vbo.bind();
	mesh.vbo.buffer_data(bufi, buffer);
	mesh.ibo.bind();
	mesh.ibo.buffer(indices);
	mesh.vao.unbind();
	mesh.vbo.unbind();
	mesh.ibo.unbind();
	glLineWidth(2.0f);
	glEnable(GL_LINE_SMOOTH); // not implemented on macos :(
}
void vecFieldRenderer::render() {
	line_shad.bind();
	gl.draw_mesh(mesh, GL_LINES);
}

void vecFieldRenderer::destroy() {
	delete [] buffer;
}
