#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <R_ext/BLAS.h>

#ifndef FCONE
#define FCONE
#endif

using namespace Rcpp;

namespace {

inline std::size_t idx_dk(int d, int k, int D) {
  return static_cast<std::size_t>(d) + static_cast<std::size_t>(D) * k;
}

inline std::size_t idx_dkl(int d, int k, int l, int D, int K) {
  return static_cast<std::size_t>(d) +
         static_cast<std::size_t>(D) * (k + K * l);
}

inline std::size_t idx_t(int row, int col, int g, int K) {
  return static_cast<std::size_t>(row) +
         static_cast<std::size_t>(K) * (col + K * g);
}

}  // namespace

// Poisson group updates using donor-level residual and basis summaries.
// [[Rcpp::export]]
List flasso_compressed_rep_cd_fit_cpp(
    NumericMatrix G,
    NumericMatrix B,
    IntegerVector donor_index,
    NumericVector y,
    NumericVector lambda,
    NumericVector feature_off,
    NumericVector Q_array,
    NumericVector group_multiplier,
    double alpha,
    double eps,
    int max_iter,
    NumericVector offset_vec) {
  const int D = G.nrow();
  const int J = G.ncol();
  const int n = B.nrow();
  const int K = B.ncol();
  const int p = J * K;
  const int L = lambda.size();

  if (donor_index.size() != n || y.size() != n || offset_vec.size() != n ||
      feature_off.size() != p || Q_array.size() != static_cast<R_xlen_t>(J) * K * K ||
      group_multiplier.size() != J || L == 0)
    stop("Incompatible dimensions in the fast-TWiST solver inputs.");

  std::vector<int> donor(n);
  for (int i = 0; i < n; ++i) {
    const int d = donor_index[i] - 1;
    if (d < 0 || d >= D) stop("donor_index out of range");
    donor[i] = d;
  }

  const double* g_ptr = G.begin();
  const double* b_ptr = B.begin();
  const double* q_ptr = Q_array.begin();

  std::vector<double> S(static_cast<std::size_t>(D) * K, 0.0);
  std::vector<double> M(static_cast<std::size_t>(D) * K * K, 0.0);
  for (int i = 0; i < n; ++i) {
    const int d = donor[i];
    for (int k = 0; k < K; ++k) {
      const double bik = b_ptr[static_cast<std::size_t>(i) + static_cast<std::size_t>(n) * k];
      S[idx_dk(d, k, D)] += bik;
      for (int l = 0; l < K; ++l) {
        const double bil = b_ptr[static_cast<std::size_t>(i) + static_cast<std::size_t>(n) * l];
        M[idx_dkl(d, k, l, D, K)] += bik * bil;
      }
    }
  }

  NumericMatrix beta(p, L);
  NumericVector intercept(L);
  IntegerVector iter(L);
  LogicalVector converged(L, false);

  double ysum = 0.0;
  double exposure = 0.0;
  for (int i = 0; i < n; ++i) {
    ysum += y[i];
    exposure += std::exp(offset_vec[i]);
  }
  const double initial_intercept =
      std::log(std::max(ysum, 1e-300) / std::max(exposure, 1e-300));

  std::vector<double> a(p, 0.0), bcoef(p, 0.0);
  std::vector<double> b_new(K, 0.0), delta(K, 0.0), z(K, 0.0);
  std::vector<double> q_delta(K, 0.0), gh(K, 0.0);
  std::vector<double> eta_cache(n, 0.0);
  std::vector<double> H(static_cast<std::size_t>(D) * K, 0.0);
  std::vector<double> F(static_cast<std::size_t>(D) * K, 0.0);
  std::vector<int> active(J, 0);

  double a0 = initial_intercept;
  double feature_offset = 0.0;
  int total_iter = 0;
  double rsum = 0.0;

  auto fill_eta_cache = [&]() {
    double max_eta = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
      double eta = offset_vec[i] + a0 + feature_offset;
      for (int k = 0; k < K; ++k) {
        eta += b_ptr[static_cast<std::size_t>(i) + static_cast<std::size_t>(n) * k] *
               F[idx_dk(donor[i], k, D)];
      }
      eta_cache[i] = eta;
      max_eta = std::max(max_eta, eta);
    }
    return max_eta;
  };

  auto rebuild_residual_stats = [&](double v) {
    std::fill(H.begin(), H.end(), 0.0);
    rsum = 0.0;
    for (int i = 0; i < n; ++i) {
      const int d = donor[i];
      const double eta = eta_cache[i];
      const double mu = std::exp(eta);
      const double ri = (y[i] - mu) / v;
      rsum += ri;
      for (int k = 0; k < K; ++k) {
        H[idx_dk(d, k, D)] +=
            b_ptr[static_cast<std::size_t>(i) + static_cast<std::size_t>(n) * k] * ri;
      }
    }
  };

  auto apply_group_delta = [&](int g, const std::vector<double>& del) {
    std::fill(q_delta.begin(), q_delta.end(), 0.0);
    double off_delta = 0.0;
    for (int h = 0; h < K; ++h) {
      off_delta += feature_off[g * K + h] * del[h];
      for (int l = 0; l < K; ++l) {
        q_delta[l] += q_ptr[idx_t(l, h, g, K)] * del[h];
      }
    }

    feature_offset += off_delta;
    double rsum_shift = off_delta * static_cast<double>(n);
    for (int l = 0; l < K; ++l) {
      const double qd = q_delta[l];
      if (qd == 0.0) continue;
      for (int d = 0; d < D; ++d) {
        const double gd = g_ptr[static_cast<std::size_t>(d) + static_cast<std::size_t>(D) * g];
        F[idx_dk(d, l, D)] += gd * q_delta[l];
        rsum_shift += gd * S[idx_dk(d, l, D)] * qd;
      }
    }
    rsum -= rsum_shift;

    if (off_delta != 0.0) {
      for (int k = 0; k < K; ++k) {
        for (int d = 0; d < D; ++d) {
          H[idx_dk(d, k, D)] -= off_delta * S[idx_dk(d, k, D)];
        }
      }
    }
    for (int l = 0; l < K; ++l) {
      const double qd = q_delta[l];
      if (qd == 0.0) continue;
      for (int k = 0; k < K; ++k) {
        for (int d = 0; d < D; ++d) {
          const double gd = g_ptr[static_cast<std::size_t>(d) + static_cast<std::size_t>(D) * g];
          H[idx_dk(d, k, D)] -= gd * M[idx_dkl(d, k, l, D, K)] * qd;
        }
      }
    }
  };

  auto fit_group = [&](int g, int ell, double v, double& max_change) {
    for (int l = 0; l < K; ++l) {
      double raw = 0.0;
      for (int d = 0; d < D; ++d) {
        raw += g_ptr[static_cast<std::size_t>(d) + static_cast<std::size_t>(D) * g] *
               H[idx_dk(d, l, D)];
      }
      gh[l] = raw;
    }

    for (int h = 0; h < K; ++h) {
      double cp = feature_off[g * K + h] * rsum;
      for (int l = 0; l < K; ++l) {
        cp += q_ptr[idx_t(l, h, g, K)] * gh[l];
      }
      z[h] = cp / static_cast<double>(n) + a[g * K + h];
    }

    double z_ss = 0.0;
    for (int h = 0; h < K; ++h) z_ss += z[h] * z[h];
    const double z_norm = std::sqrt(z_ss);
    std::fill(b_new.begin(), b_new.end(), 0.0);
    const double mult = group_multiplier[g];
    const double lam1 = lambda[ell] * mult * alpha;
    const double lam2 = lambda[ell] * mult * (1.0 - alpha);
    if (z_norm > 0.0) {
      const double len = std::max(v * z_norm - lam1, 0.0) / (v * (1.0 + lam2));
      for (int h = 0; h < K; ++h) b_new[h] = len * z[h] / z_norm;
    }

    bool changed = false;
    for (int h = 0; h < K; ++h) {
      const int col = g * K + h;
      delta[h] = b_new[h] - a[col];
      if (delta[h] != 0.0) {
        changed = true;
        max_change = std::max(max_change, std::abs(delta[h]));
      }
      bcoef[col] = b_new[h];
    }
    if (changed) apply_group_delta(g, delta);
    return bcoef[g * K] != 0.0;
  };

  // Carry the fitted coefficients forward along the penalty path.
  for (int ell = 0; ell < L; ++ell) {
    while (total_iter < max_iter) {
      bool active_converged = false;
      double v = 1.0;

      while (total_iter < max_iter) {
        iter[ell]++;
        total_iter++;

        const double max_eta = fill_eta_cache();
        v = std::exp(max_eta);
        rebuild_residual_stats(v);

        const double shift0 = rsum / static_cast<double>(n);
        double max_change = std::abs(shift0);
        a0 += shift0;
        rsum -= shift0 * n;
        for (int k = 0; k < K; ++k) {
          for (int d = 0; d < D; ++d) {
            H[idx_dk(d, k, D)] -= shift0 * S[idx_dk(d, k, D)];
          }
        }

        for (int g = 0; g < J; ++g) {
          if (!active[g]) continue;
          fit_group(g, ell, v, max_change);
        }

        a = bcoef;
        if (max_change < eps) {
          active_converged = true;
          break;
        }
      }

      if (!active_converged) break;

      int violations = 0;
      double max_change_scan = 0.0;
      for (int g = 0; g < J; ++g) {
        if (active[g]) continue;
        const bool first_nonzero = fit_group(g, ell, v, max_change_scan);
        if (first_nonzero) {
          active[g] = 1;
          violations++;
        }
      }

      a = bcoef;
      if (violations == 0) {
        converged[ell] = true;
        break;
      }
    }

    intercept[ell] = a0;
    for (int j = 0; j < p; ++j) beta(j, ell) = bcoef[j];
  }

  return List::create(
      _["beta_ortho"] = beta,
      _["intercept_ortho"] = intercept,
      _["iter"] = iter,
      _["converged"] = converged);
}

// Compute grpreg::maxgrad for the donor-level representation without
// materializing the full design. The accumulation order matches grpreg:
// group by group, column by column, then cell by cell.
// [[Rcpp::export]]
double flasso_rep_maxgrad_cell_cpp(
    NumericMatrix G,
    NumericMatrix B,
    IntegerVector donor_index,
    NumericVector r,
    NumericVector center,
    NumericVector scale,
    NumericVector raw_T_array,
    NumericVector group_multiplier) {
  const int D = G.nrow();
  const int J = G.ncol();
  const int n = B.nrow();
  const int K = B.ncol();
  const int p = J * K;

  if (donor_index.size() != n || r.size() != n || center.size() != p || scale.size() != p ||
      raw_T_array.size() != static_cast<R_xlen_t>(J) * K * K || group_multiplier.size() != J)
    stop("Incompatible dimensions in the fast-TWiST gradient inputs.");

  std::vector<int> donor(n);
  for (int i = 0; i < n; ++i) {
    const int d = donor_index[i] - 1;
    if (d < 0 || d >= D) stop("donor_index out of range");
    donor[i] = d;
  }

  const double* g_ptr = G.begin();
  const double* b_ptr = B.begin();
  const double* center_ptr = center.begin();
  const double* scale_ptr = scale.begin();
  const double* raw_t_ptr = raw_T_array.begin();
  const double* r_ptr = r.begin();

  std::vector<double> std_input_block(static_cast<std::size_t>(n) * K, 0.0);
  std::vector<double> xblock(static_cast<std::size_t>(n) * K, 0.0);
  double zmax = 0.0;

  for (int g = 0; g < J; ++g) {
    if (group_multiplier[g] == 0.0) continue;
    for (int l = 0; l < K; ++l) {
      const int raw_col = g * K + l;
      for (int i = 0; i < n; ++i) {
        const double gd =
            g_ptr[static_cast<std::size_t>(donor[i]) + static_cast<std::size_t>(D) * g];
        const double raw = gd * b_ptr[static_cast<std::size_t>(i) + static_cast<std::size_t>(n) * l];
        std_input_block[static_cast<std::size_t>(i) + static_cast<std::size_t>(n) * l] =
            (raw - center_ptr[raw_col]) / scale_ptr[raw_col];
      }
    }

    const char trans = 'N';
    const double one = 1.0;
    const double zero = 0.0;
    F77_CALL(dgemm)(&trans, &trans, &n, &K, &K, &one,
                    std_input_block.data(), &n,
                    raw_t_ptr + static_cast<std::size_t>(K) * K * g, &K,
                    &zero, xblock.data(), &n FCONE FCONE);

    double z_ss = 0.0;
    for (int h = 0; h < K; ++h) {
      double cp = 0.0;
      const std::size_t col_off = static_cast<std::size_t>(n) * h;
      for (int i = 0; i < n; ++i) cp += xblock[col_off + i] * r_ptr[i];
      z_ss += cp * cp;
    }
    const double z = std::sqrt(z_ss) / group_multiplier[g];
    zmax = std::max(zmax, z);
  }

  return zmax;
}
