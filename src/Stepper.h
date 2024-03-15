#ifndef STEPPER_H
#define STEPPER_H

#include <flgl/glm.h>

struct Stepper {
	int x1, y1, x2, y2;
	Stepper(int x1, int y1, int x2, int y2);
	class StepIterator {
		int x1, y1, x2, y2;
		glm::vec2 pos;
		float* u, * v;
		float slope, step;
		bool complete;
	public:
		StepIterator(int x1, int y1, int x2, int y2);
		StepIterator(bool);
		glm::ivec2 operator*() const;
		StepIterator& operator++();
		bool operator==(StepIterator const& other) const;
		bool operator!=(StepIterator const& other) const;
	};
	const StepIterator begin() const;
	const StepIterator end() const;
};

#endif
