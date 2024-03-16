#include "FieldRenderer.h"

FieldRenderer::FieldRenderer(Field const& f) : field(f), n(f.solver.dim()) {}
FieldRenderer::~FieldRenderer() {}
