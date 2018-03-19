#include "range.hpp"

range::bound &range::bound::operator-() {
    value = -value;
    return *this;
}
