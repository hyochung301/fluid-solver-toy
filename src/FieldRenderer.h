#ifndef FIELD_RENDERER_H
#define FIELD_RENDERER_H
#include <flgl.h>
#include "Field.h"

struct FieldRenderer {
	Mesh<Vt_classic> quad;
	Shader field_shad;
	Texture field_tex;
	const int n;
	Field const& field;

	FieldRenderer(Field const& f) : field(f), n(f.solver.dim()) {}
	void init() {
		quad = DefaultMeshes::tile<Vt_classic>();
		field_shad = Shader::from_source("passthrough_vert", "tex");
		field_tex.create();
		field_tex.bind();
		field_tex.pixelate(false);
    	field_tex.unbind();
	}
	void buffer_texture() {
		field_tex.bind();
    	field_tex.alloc(GL_RG32F, n, n, GL_RG, GL_FLOAT, field.buff());
		field_tex.unbind();
	}
	void render() {
		gl.clear();
		field_tex.bind();
		field_tex.bind_to_unit(0);
		field_shad.bind();
		field_shad.uInt("uTexslot", 0);
		gl.draw_mesh(quad);
		field_tex.unbind();
		field_shad.unbind();
	}
};

#endif
