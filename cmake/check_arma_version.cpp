#include <iostream>
#include <armadillo>

int main() {
  std::cout << ARMA_VERSION_MAJOR << "."
            << ARMA_VERSION_MINOR << "."
            << ARMA_VERSION_PATCH;
  return 1;
}
