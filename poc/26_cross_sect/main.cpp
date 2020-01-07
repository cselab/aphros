#undef NDEBUG
#include <cassert>
#include <iostream>
#include <string>
#include <memory>

using Scal = double;
using Vect = GVect<Scal, 3>;

int main(int argc, const char** argv) {
  std::cout << Vect(5);
}
