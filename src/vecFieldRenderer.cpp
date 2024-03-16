#include "vecFieldRenderer.h"
#include <flgl/tools.h>
#include <random>

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
void vecFieldRenderer::prepare() {
	const float size = 2.0f / (n - 1); // Calculate the size of each step in the grid
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.8, 1.2);
	const float coeff = dis(gen) * (20.0f / (n - 1));
    std::vector<uint32_t> indices;
    float* const u = field.solver.x_buffer();
    float* const v = field.solver.y_buffer();

    // Fill vertices buffer
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < n; ++i) {

            // Scale positions to [-1, 1]
            buffer[i+j*n].pos = glm::vec2(-1.0f + size * i, -1.0f + size * j);

            // Calculate and scale vector
            glm::vec2 vec = glm::vec2(u[i + j * n], v[i + j * n]);
            buffer[i+j*n].vec = buffer[i+j*n].pos + coeff * vec; // Adjust vector by the given coefficient and add to position

            // Add line indices (each vertex to its corresponding vector tip)
            size_t index = i + j * n;
            indices.push_back(2*index); // start of line (position)
            indices.push_back(2*index+1); // end of line (position + vector)
        }
    }

    mesh.vao.bind();
	mesh.vbo.bind();
	mesh.vbo.buffer_data(n*n, buffer);
	mesh.ibo.bind();
	mesh.ibo.buffer(indices);
	mesh.vao.unbind();
	mesh.vbo.unbind();
	mesh.ibo.unbind();
}
void vecFieldRenderer::render() {
	gl.clear();
	line_shad.bind();
	gl.draw_mesh(mesh, GL_LINES);
}

void vecFieldRenderer::destroy() {
	delete [] buffer;
}
