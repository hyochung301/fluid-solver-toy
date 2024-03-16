#ifndef FIELD_RENDERER_H
#define FIELD_RENDERER_H
#include <flgl.h>
#include "Field.h"

struct FieldRenderer {
	const int n;
	Field const& field;

	FieldRenderer(Field const& f);
	virtual void init() 	= 0;
	virtual void prepare() 	= 0;
	virtual void render() 	= 0;
	virtual void destroy() 	= 0;
};

#endif
