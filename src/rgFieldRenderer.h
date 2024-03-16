#ifndef RG_FIELD_RENDERER_H
#define RG_FIELD_RENDERER_H
#include "FieldRenderer.h"
#include <flgl.h>
#include "Field.h"
#include "FieldRenderer.h"

struct rgFieldRenderer : public FieldRenderer {
	Mesh<Vt_classic> quad;
	Shader field_shad;
	Texture field_tex;

	rgFieldRenderer(Field const& f);
	void init() override final;
	void prepare() override final;
	void render() override final;
	void destroy() override final;
};

#endif
