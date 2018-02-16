#include <iostream>
#include <cstdlib>

int main(int c, char** v) {
  // parse arg
  std::string a = v[1];
  std::string b = v[2];
  std::cout 
      << "a=" << a << " "
      << "b=" << b << " "
      << std::endl;
  std::cerr 
      << "err a=" << a << " "
      << "err b=" << b << " "
      <<  std::endl;

  if (a != "a") {
    return 11;
  }
  if (b != "b") {
    return 12;
  }
  return 0;
}
