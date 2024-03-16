#include "rgFieldRenderer.h"
#include <flgl/tools.h>


rgFieldRenderer::rgFieldRenderer(Field const& f) : FieldRenderer(f) {}
void rgFieldRenderer::init() {
	quad = DefaultMeshes::tile<Vt_classic>();
	field_shad = Shader::from_source("passthrough_vert", "tex");
	field_tex.create();
	field_tex.bind();
	field_tex.pixelate(false);
	field_tex.unbind();
}
void rgFieldRenderer::prepare() {
	field_tex.bind();
	field_tex.alloc(GL_RG32F, n, n, GL_RG, GL_FLOAT, field.buff());
	field_tex.unbind();
}
void rgFieldRenderer::render() {
	gl.clear();
	field_tex.bind();
	field_tex.bind_to_unit(0);
	field_shad.bind();
	field_shad.uInt("uTexslot", 0);
	gl.draw_mesh(quad);
	field_tex.unbind();
	field_shad.unbind();
}

void rgFieldRenderer::destroy() {
}
