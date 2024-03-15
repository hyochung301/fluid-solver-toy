#include "Stepper.h"
#include <cmath>

Stepper::Stepper(int a, int b, int c, int d)
:	x1(a), y1(b), x2(c), y2(d) {}

const Stepper::StepIterator Stepper::begin() const {
	return StepIterator(this->x1,this->y1,this->x2,this->y2);
};

const Stepper::StepIterator Stepper::end() const {
	return StepIterator(true);
};

Stepper::StepIterator::StepIterator(int a, int b, int c, int d)
:	x1(a), y1(b), x2(c), y2(d), pos(x1,y1), u(0), v(0), slope(0), step(1.f), complete(false)
{
	bool back = false;
	if (std::fabs(x2 - x1) > std::fabs(y2 - y1)) {
		u = &pos.x; v = &pos.y;
		slope = x2==x1 ? 0.f : ((float)(y2 - y1)) / ((float)(x2 - x1));
		if (x2 < x1) back = true;
	} else {
		u = &pos.y; v = &pos.x;
		slope = y2==y1 ? 0.f : ((float)(x2 - x1)) / ((float)(y2 - y1));
		if (y2 < y1) back = true;
	}
	if (!back) {return;}
	slope *= -1;
	step = -1.f;
}

Stepper::StepIterator::StepIterator(bool c){complete = c;}

glm::ivec2 Stepper::StepIterator::operator*() const {
	return {
		(int)(std::round(pos.x)),
		(int)(std::round(pos.y))
	};
}

Stepper::StepIterator& Stepper::StepIterator::operator++() {
	if ((*(*this)) == glm::ivec2(x2,y2)) {
		this->complete = true; return *this;
	}
	(*u) = (*u) + step;
	(*v) = (*v) + slope;
	return *this;
}

bool Stepper::StepIterator::operator==(Stepper::StepIterator const& other) const {
	return this->complete == other.complete;
}

bool Stepper::StepIterator::operator!=(Stepper::StepIterator const& other) const {
	return this->complete != other.complete;
}
