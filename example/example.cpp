#include <iostream>
#include <algorithm>

#include <triqs/arrays/arpack/arpack_worker.hpp>

// This example shows how to partially diagonalize a large sparse symmetric
// matrix and find a number of its low-lying eigenvalues.

// Size of the matrix
const int N = 10000;

// We are going to use a band matrix with this bandwidth
const int bandwidth = 5;

// The number of low-lying eigenvalues we want to compute
const int N_ev = 10;

using namespace triqs::arrays::arpack;

int main(int argc, char* argv[]) {

 // Construct a worker object for the symmetric case.
 // Other options would be
 // * `arpack_worker<Asymmetric>' for general real matrices.
 // * `arpack_worker<Complex>' for general complex matrices.
 arpack_worker<Symmetric> worker(N);

 // Linear operator representing multiplication of a given vector by our matrix
 // Arguments `from_index' and `to_index` provide an alternative way to access
 // the ARPACK workspace vectors:
 //  worker.workspace_vector(from_index)
 //  worker.workspace_vector(to_index)
 auto matrix = [](vector_const_view<double> from, int from_index,
                  vector_view<double> to, int to_index) {
  to() = 0; // Clear result

  // to_i = \sum_j A_{ij} from_j
  // A_ij = |i-j| / (1 + i + j), |i-j| <= bandwidth, zero otherwise
  for(int i = 0; i < N; ++i) {
   int j_min = std::max(0, i - bandwidth);
   int j_max = std::min(N, i + bandwidth);
   for(int j = j_min; j <= j_max; ++j) {
    to(i) += double(std::abs(i - j)) / (1 + i + j) * from(j);
   }
  }
 };

 // Specify parameters for the worker
 using params_t = arpack_worker<Symmetric>::params_t;
 params_t params(N_ev,               // Number of low-lying eigenvalues
                 params_t::Smallest, // We want the smallest eigenvalues
                 true                // Yes, we want the eigenvectors (Ritz vectors) as well
                );

 // Run diagonalization!
 worker(matrix, params);

 // Print found eigenvalues
 std::cout << "Eigenvalues (Ritz values):" << std::endl;
 std::cout << worker.eigenvalues() << std::endl;

 // Check A*v = \lambda*v
 auto const& lambda = worker.eigenvalues();
 auto const& v = worker.eigenvectors();
 vector<double> lhs(N), rhs(N);

 for(int i = 0; i < N_ev; ++i) {          // For each eigenpair ...
  matrix(v(range(), i), 0, lhs, 0);       // calculate A*v
  rhs = lambda(i) * v(range(), i);        // and \lambda*v

  std::cout << i << ": deviation = " << norm2_sqr(rhs - lhs) / (N*N) << std::endl;
 }

 // Print some computation statistics
 auto stats = worker.stats();

 std::cout << "Number of Arnoldi update iterations: " << stats.n_iter << std::endl;
 std::cout << "Number of 'converged' Ritz values: " << stats.n_converged << std::endl;
 std::cout << "Total number of OP*x operations: " << stats.n_op_x_operations << std::endl;
 std::cout << "Total number of steps of re-orthogonalization: " << stats.n_reorth_steps << std::endl;

 return 0;
}
