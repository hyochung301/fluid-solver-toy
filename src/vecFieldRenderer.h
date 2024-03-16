#ifndef VECTOR_FIELD_RENDERER_H
#define VECTOR_FIELD_RENDERER_H
#include "FieldRenderer.h"
#include <flgl.h>
#include "Field.h"
#include "FieldRenderer.h"

struct vecfield_vert {
	glm::vec2 pos;
	glm::vec2 vec;
};

struct vecFieldRenderer : public FieldRenderer {

	vecfield_vert* buffer;

	Shader line_shad;
	Mesh<vecfield_vert> mesh;

	vecFieldRenderer(Field const& f);
	void init() override final;
	void prepare() override final;
	void render() override final;
	void destroy() override final;
};

#endif
