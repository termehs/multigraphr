#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector calculateQSeq(double m, NumericVector deg_seq, String model) {
  int n = deg_seq.size();
  int sum_deg = sum(deg_seq);
  // Check if the sum of deg_seq is even
  if (sum_deg % 2 == 1) {
    stop("not a graphical degree sequence, the sum is not equal to 2m.");
  }

  NumericMatrix Q_mat(n, n);
  NumericVector Q_seq;

  if (model == "IEAS") {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i == j) {
          Q_mat(i, j) = deg_seq[i] * (deg_seq[i] - 1) / (2 * m * (2 * m - 1));
        } else {
          Q_mat(i, j) = 2 * deg_seq[i] * deg_seq[j] / (2 * m * (2 * m - 1));
        }
      }
    }
  } else if (model == "ISA") {
    NumericVector p_seq = deg_seq / (2 * m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i == j) {
          Q_mat(i, j) = pow(p_seq[i], 2);
        } else {
          Q_mat(i, j) = 2 * p_seq[i] * p_seq[j];
        }
      }
    }
  } else {
    stop("model must be one of 'IEAS' or 'ISA'.");
  }

  // Extract lower triangle including diagonal as a vector
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      Q_seq.push_back(Q_mat(i, j));
    }
  }

  return Q_seq;
}
