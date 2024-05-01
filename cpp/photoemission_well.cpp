// Before execution, please check that you properly prepared the constants at the beginning of main().

#define USE_MATH_DEFINES
#include "photoemission_well.hpp"
#include <cmath>
#include <complex>
#include <cstdio>

using namespace std;
int main() {
  // Constants
  /// System
  double half_well_length = 10;
  double half_system_length = 1000;
  double potential_well = -1.0;
  int numPoints_per_unit = 10;
  double k_convergence_threshold = 1e-10;
  /// Photoemission simulations
  double Ekin_min = .2;
  double Ekin_delta = 0.02;
  int numEkin = 41;
  int numLEStates_sim = 1; // should be < number of bound states
  double deltaH_amplitude = 1;
  double propagation_length = 100;
  double integration_range = 5;
  double deltak_coefficient = 0.2;
  double photoelectron_range_min = 20;
  double photoelectron_range_max = 80;
  /// Export of the results
  const char *filePath_V = "./output/Potential.dat";
  const char *filePath_Psi_bound = "./output/Psi_bound.dat";
  const char *fileFormat_Photoelectron_amplitude = "./output/Photoelectron_amplitude_%d.dat";
  const char *fileFormat_Psi_unbound_cosine = "./output/Psi_unbound_cosine_%d_%d.dat";
  const char *fileFormat_Psi_unbound_cosine_wfn = "./output/Psi_unbound_cosine_wfn_%d_%d.dat";
  const char *fileFormat_Psi_unbound_sine = "./output/Psi_unbound_sine_%d_%d.dat";
  const char *fileFormat_Psi_unbound_sine_wfn = "./output/Psi_unbound_sine_wfn_%d_%d.dat";
  const char *fileFormat_Psi_pw_wfn = "./output/Psi_pw_wfn_%d_%d.dat";
  const char *fileFormat_Psi_ex = "./output/Psi_ex_%d_%d.dat";
  int numUnboundStates_export = 1;
  /// Physical constants
  const double fine_structure_constant = 7.297e-3;
  // utility variables
  // k is skipped because k is used as the wavevector variable
  int i, j, l, m;
  const double avoid_divergence = 1e-5;
  // shortened variables
  double a = half_well_length;
  double b = half_system_length;
  double V0 = potential_well;

  // Step 1: prepare the potential
  // x range: [-vacuum_length, period_slab*num_layers+vacuum_length]
  printf("Step 1: Preparation of the slab potential\n");
  double dx = 1.0 / numPoints_per_unit;
  double system_length = 2.0 * half_system_length;
  int numPoints = round(system_length * numPoints_per_unit) + 1;
  printf("System length: %.1f (Bohr), Number of points: %d\n", system_length, numPoints);

  double *V = alloc_dvector(numPoints);
  double *x = alloc_dvector(numPoints);

  for (i = 0; i < numPoints; i++) {
    x[i] = dx * i - half_system_length;
    if (-half_well_length <= x[i] && x[i] <= half_well_length) {
      V[i] = potential_well;
    } else {
      V[i] = 0.0;
    }
  }

  FILE *export_V = fopen(filePath_V, "w");
  fprintf(export_V, "# Model potential\n");
  fprintf(export_V, "# x V(x)\n");
  for (i = 0; i < numPoints; i++) {
    fprintf(export_V, "%.2f %.3f\n", x[i], V[i]);
  }
  fclose(export_V);
  printf("The slab potential was exported to %s\n", filePath_V);

  // Step 2: Obtain bound eigenstate(s)
  printf("Step 2: Bound eigenstates\n");
  double K_max = sqrt(-2 * V0);
  int K_count_total = ceil(K_max / (0.5 * M_PI / a));
  printf("%d eigenstates will be found in total\n", K_count_total);
  double **Psis_bound = alloc_dmatrix(K_count_total, numPoints);
  double *Ks_bound = alloc_dvector(K_count_total);
  double *kappas_bound = alloc_dvector(K_count_total);
  double *Es_bound = alloc_dvector(K_count_total);
  double *alphas_bound = alloc_dvector(K_count_total);

  printf("Step 2-1: cosine type\n");
  int K_count_cosine = ceil(K_max / (M_PI / a));
  printf("%d eigenstates will be found\n", K_count_cosine);
  for (i = 0; i < K_count_cosine; i++) {
    // K range: [i*PI, (i+1/2)*PI]
    double K_left = i * M_PI / a + avoid_divergence;
    double K_right = (i + 0.5) * M_PI / a - avoid_divergence;
    if (K_right > K_max) {
      K_right = K_max - avoid_divergence;
    }
    double diffFunc_left = diffFunc_bound_cosine(K_left, V0, a);
    double diffFunc_right = diffFunc_bound_cosine(K_right, V0, a);
    double K_center, diffFunc_center;
    if (diffFunc_left * diffFunc_right < 0) {
      while (K_right - K_left > k_convergence_threshold) {
        K_center = 0.5 * (K_left + K_right);
        diffFunc_center = diffFunc_bound_cosine(K_center, V0, a);
        if (diffFunc_left * diffFunc_center < 0) {
          K_right = K_center;
          diffFunc_right = diffFunc_center;
        } else {
          K_left = K_center;
          diffFunc_left = diffFunc_center;
        }
      }
    } else {
      printf("Bisection error\n");
      return -1;
    }
    printf("Eigenstate[%d]: K=%.7f\n", i, K_center);
    double K = K_center;
    double E = 0.5 * K * K + V0;
    double kappa = sqrt(-2.0 * E);
    double alpha = cos(K * a) * exp(kappa * a);
    double norm = 2.0 * (a / 2.0 + 1.0 / 4.0 / K * sin(2.0 * K * a) + alpha * alpha / 2.0 / kappa * exp(-2.0 * kappa * a));
    Ks_bound[i * 2] = K;
    kappas_bound[i * 2] = kappa;
    Es_bound[i * 2] = E;
    alphas_bound[i * 2] = alpha;
    for (j = 0; j < numPoints; j++) {
      if (x[j] < -a) {
        Psis_bound[i * 2][j] = alpha * exp(kappa * x[j]) / sqrt(norm);
      } else if (x[j] < a) {
        Psis_bound[i * 2][j] = cos(K * x[j]) / sqrt(norm);
      } else {
        Psis_bound[i * 2][j] = alpha * exp(-kappa * x[j]) / sqrt(norm);
      }
    }
  }

  printf("Step 2-2: sine type\n");
  int K_count_sine = ceil(K_max / (M_PI / a) - 0.5);
  printf("%d eigenstates will be found\n", K_count_sine);
  for (i = 0; i < K_count_sine; i++) {
    // K range: [(i+1/2)*PI, (i+1)*PI]
    double K_left = (i + 0.5) * M_PI / a + avoid_divergence;
    double K_right = (i + 1.0) * M_PI / a - avoid_divergence;
    if (K_right > K_max) {
      K_right = K_max - avoid_divergence;
    }
    double diffFunc_left = diffFunc_bound_sine(K_left, V0, a);
    double diffFunc_right = diffFunc_bound_sine(K_right, V0, a);
    double K_center, diffFunc_center;
    if (diffFunc_left * diffFunc_right < 0) {
      while (K_right - K_left > k_convergence_threshold) {
        K_center = 0.5 * (K_left + K_right);
        diffFunc_center = diffFunc_bound_sine(K_center, V0, a);
        if (diffFunc_left * diffFunc_center < 0) {
          K_right = K_center;
          diffFunc_right = diffFunc_center;
        } else {
          K_left = K_center;
          diffFunc_left = diffFunc_center;
        }
      }
    } else {
      printf("Bisection error\n");
      return -1;
    }
    printf("Eigenstate[%d]: K=%.7f\n", i, K_center);
    double K = K_center;
    double E = 0.5 * K * K + V0;
    double kappa = sqrt(-2.0 * E);
    double alpha = sin(K * a) * exp(kappa * a);
    double norm = 2.0 * (a / 2.0 - 1.0 / 4.0 / K * sin(2.0 * K * a) + alpha * alpha / kappa / 2.0 * exp(-2.0 * kappa * a));
    Ks_bound[i * 2 + 1] = K;
    kappas_bound[i * 2 + 1] = kappa;
    Es_bound[i * 2 + 1] = E;
    alphas_bound[i * 2 + 1] = alpha;
    for (j = 0; j < numPoints; j++) {
      if (x[j] < -a) {
        Psis_bound[i * 2 + 1][j] = -alpha * exp(kappa * x[j]) / sqrt(norm);
      } else if (x[j] < a) {
        Psis_bound[i * 2 + 1][j] = sin(K * x[j]) / sqrt(norm);
      } else {
        Psis_bound[i * 2 + 1][j] = alpha * exp(-kappa * x[j]) / sqrt(norm);
      }
    }
  }

  FILE *export_Psi_bound = fopen(filePath_Psi_bound, "w");
  fprintf(export_Psi_bound, "# Bound eigenstates\n");
  fprintf(export_Psi_bound, "# Eigenvalues:\n");
  fprintf(export_Psi_bound, "# ");
  for (i = 0; i < K_count_total; i++) {
    fprintf(export_Psi_bound, "%.3f ", Es_bound[i]);
  }
  fprintf(export_Psi_bound, "\n");
  fprintf(export_Psi_bound, "# x ");
  for (i = 0; i < K_count_total; i++) {
    fprintf(export_Psi_bound, "Psi[%d](x) ", i);
  }
  fprintf(export_Psi_bound, "\n");
  for (i = 0; i < numPoints; i++) {
    fprintf(export_Psi_bound, "%.2f ", x[i]);
    for (j = 0; j < K_count_total; j++) {
      fprintf(export_Psi_bound, "%.5f ", Psis_bound[j][i]);
    }
    fprintf(export_Psi_bound, "\n");
  }
  fclose(export_Psi_bound);
  printf("The bound eigenstates were exported to %s\n", filePath_Psi_bound);

  // Step 3: Direct photoemission simulations
  printf("Step 3: Direct photoemission simulations\n");
  complex<double> *Psi_ex = alloc_zvector(numPoints);
  double *Psi_unbound = alloc_dvector(numPoints);
  double *dPsi_bound = alloc_dvector(numPoints);
  for (i = 0; i < numLEStates_sim; i++) {
    printf("Simulation for eigenstate #%d\n", i);
    char filePath_Photoelectron_amplitude[1024];
    sprintf(filePath_Photoelectron_amplitude, fileFormat_Photoelectron_amplitude, i);
    FILE *export_Photoelectron_amplitude = fopen(filePath_Photoelectron_amplitude, "w");
    fprintf(export_Photoelectron_amplitude, "# Kinetic energy, Photoelectron amplitude (average), Photoelectron amplitude (stdev), Photoelectron amplitude (PW), Photoelectron amplitude (TR-LEED)\n");
    // differenciate
    for (j = 0; j < numPoints; j++) {
      if (j == 0) {
        dPsi_bound[j] = (Psis_bound[i][j + 1] - Psis_bound[i][j]) / dx;
      } else if (j == numPoints - 1) {
        dPsi_bound[j] = (Psis_bound[i][j] - Psis_bound[i][j - 1]) / dx;
      } else {
        dPsi_bound[j] = (Psis_bound[i][j + 1] - Psis_bound[i][j - 1]) / dx / 2;
      }
    }
    for (j = 0; j < numEkin; j++) {
      double Ekin = Ekin_min + Ekin_delta * j;
      double k = sqrt(2.0 * Ekin);
      printf("Photoelectron kinetic energy: %.3f\n", Ekin);
      printf("Photoelectron wavevector: %.3f\n", k);
      double hn = Ekin - Es_bound[i];
      double kph = hn * fine_structure_constant;
      printf("Photon energy: %.3f\n", hn);
      printf("Photon momentum: %.3f\n", kph);
      double t = propagation_length / k;
      printf("T_max: %.3f\n", t);
      double half_integration_range_k = integration_range * 2.0 * M_PI / t;
      printf("Half k range: %.6f\n", half_integration_range_k);
      double dk = M_PI / b * deltak_coefficient;
      printf("Delta k: %.6f\n", dk);
      int k_index_offset = ceil(half_integration_range_k / dk);
      int k_array_size = 2 * k_index_offset + 1;
      double *n_float_array = alloc_dvector(k_array_size);

      // cosine-type
      printf("Step A: Unbound state search: cosine-type\n");
      int k_count_cosine = 0;
      for (l = 0; l < k_array_size; l++) {
        double kl = k + (l - k_index_offset) * dk;
        n_float_array[l] = n_float_cosine(kl, V0, a, b);
        // printf("kl=%.7f, n_float=%.7f\n", kl, n_float_array[l]);
        if (l > 0 && floor(n_float_array[l]) - floor(n_float_array[l - 1]) == 1) {
          // printf("kl=%.7f, n_float=%.7f %.7f\n", kl, n_float_array[l - 1], n_float_array[l]);
          k_count_cosine++;
        }
      }
      printf("%d eigenstates will be found\n", k_count_cosine);
      double *ks_unbound_c = alloc_dvector(k_count_cosine);
      double *Ks_unbound_c = alloc_dvector(k_count_cosine);
      double *alphas_unbound_c = alloc_dvector(k_count_cosine);
      double *thetas_unbound_c = alloc_dvector(k_count_cosine);
      double *norms_unbound_c = alloc_dvector(k_count_cosine);
      int k_index = 0;
      for (l = 1; l < k_array_size; l++) {
        if (floor(n_float_array[l]) - floor(n_float_array[l - 1]) == 1) {
          double k_left = k + (l - k_index_offset - 1) * dk;
          double k_right = k + (l - k_index_offset) * dk;
          double n_float_left = n_float_cosine(k_left, V0, a, b);
          double n_float_right = n_float_cosine(k_right, V0, a, b);
          double k_center = 0.5 * (k_left + k_right);
          double n_float_center;
          while (k_right - k_left > k_convergence_threshold) {
            k_center = 0.5 * (k_left + k_right);
            n_float_center = n_float_cosine(k_center, V0, a, b);
            if (floor(n_float_center) - floor(n_float_left) == 1) {
              k_right = k_center;
              n_float_right = n_float_center;
            } else if (floor(n_float_right) - floor(n_float_center) == 1) {
              k_left = k_center;
              n_float_left = n_float_center;
            } else {
              printf("Bisection error\n");
              return -1;
            }
          }
          // printf("Eigenstate[%d]: k=%.7f\n", k_index, k_center);
          double k = k_center;
          double E = 0.5 * k * k;
          double K = sqrt(2.0 * (E - V0));
          double theta = atan(1.0 / tan(K * a) * (-k / K)) - k * a;
          double alpha = cos(K * a) / sin(k * a + theta);
          ks_unbound_c[k_index] = k;
          Ks_unbound_c[k_index] = K;
          thetas_unbound_c[k_index] = theta;
          alphas_unbound_c[k_index] = alpha;
          norms_unbound_c[k_index] = 2.0 * (a / 2.0 + 1.0 / 4.0 / K * sin(2.0 * K * a) +
                                            alpha * alpha * ((b - a) / 2.0 + 1.0 / 4.0 / k * (sin(2.0 * (k * a + theta)) - sin(2.0 * (k * b + theta)))));
          k_index++;
        }
      }
      char filePath_Psi_unbound_cosine[1024];
      sprintf(filePath_Psi_unbound_cosine, fileFormat_Psi_unbound_cosine, i, j);
      FILE *export_Psi_unbound_cosine = fopen(filePath_Psi_unbound_cosine, "w");
      fprintf(export_Psi_unbound_cosine, "# Unound eigenstates (cosine)\n");
      fprintf(export_Psi_unbound_cosine, "# index k K theta alpha\n");
      for (l = 0; l < k_count_cosine; l++) {
        fprintf(export_Psi_unbound_cosine, "%d %.7f %.7f %.7f %.7f\n", l,
                ks_unbound_c[l], Ks_unbound_c[l], thetas_unbound_c[l], alphas_unbound_c[l]);
      }
      fclose(export_Psi_unbound_cosine);
      char filePath_Psi_unbound_cosine_wfn[1024];
      sprintf(filePath_Psi_unbound_cosine_wfn, fileFormat_Psi_unbound_cosine_wfn, i, j);
      FILE *export_Psi_unbound_cosine_wfn = fopen(filePath_Psi_unbound_cosine_wfn, "w");
      fprintf(export_Psi_unbound_cosine_wfn, "# Unound eigenstates (cosine)\n");
      fprintf(export_Psi_unbound_cosine_wfn, "# x ");
      for (l = 0; l < numUnboundStates_export; l++) {
        fprintf(export_Psi_unbound_cosine_wfn, "Psi[%d](x) ", l);
      }
      fprintf(export_Psi_unbound_cosine_wfn, "\n");
      for (l = 0; l < numPoints; l++) {
        fprintf(export_Psi_unbound_cosine_wfn, "%.2f ", x[l]);
        for (m = 0; m < numUnboundStates_export; m++) {
          fprintf(export_Psi_unbound_cosine_wfn, "%.7f ",
                  unbound_wfn_cosine(x[l], a, b, ks_unbound_c[m], Ks_unbound_c[m], alphas_unbound_c[m], thetas_unbound_c[m], norms_unbound_c[m]));
        }
        fprintf(export_Psi_unbound_cosine_wfn, "\n");
      }
      fclose(export_Psi_unbound_cosine_wfn);
      if (k_index == k_count_cosine) {
        printf("All unbound eigenstates (cosine type) are found\n");
      } else {
        printf("Unbound eigenstates (cosine type) are missing\n");
        return -1;
      }

      // sine-type
      printf("Step B: Unbound state search: sine-type\n");
      int k_count_sine = 0;
      for (l = 0; l < k_array_size; l++) {
        double kl = k + (l - k_index_offset) * dk;
        n_float_array[l] = n_float_sine(kl, V0, a, b);
        // printf("kl=%.7f, n_float=%.7f\n", kl, n_float_array[l]);
        if (l > 0 && floor(n_float_array[l]) - floor(n_float_array[l - 1]) == 1) {
          // printf("kl=%.7f, n_float=%.7f %.7f\n", kl, n_float_array[l - 1], n_float_array[l]);
          k_count_sine++;
        }
      }
      printf("%d eigenstates will be found\n", k_count_sine);
      k_index = 0;
      double *ks_unbound_s = alloc_dvector(k_count_sine);
      double *Ks_unbound_s = alloc_dvector(k_count_sine);
      double *alphas_unbound_s = alloc_dvector(k_count_sine);
      double *thetas_unbound_s = alloc_dvector(k_count_sine);
      double *norms_unbound_s = alloc_dvector(k_count_sine);
      for (l = 1; l < k_array_size; l++) {
        if (floor(n_float_array[l]) - floor(n_float_array[l - 1]) == 1) {
          double k_left = k + (l - k_index_offset - 1) * dk;
          double k_right = k + (l - k_index_offset) * dk;
          double n_float_left = n_float_sine(k_left, V0, a, b);
          double n_float_right = n_float_sine(k_right, V0, a, b);
          double k_center = 0.5 * (k_left + k_right);
          double n_float_center;
          while (k_right - k_left > k_convergence_threshold) {
            k_center = 0.5 * (k_left + k_right);
            n_float_center = n_float_sine(k_center, V0, a, b);
            if (floor(n_float_center) - floor(n_float_left) == 1) {
              k_right = k_center;
              n_float_right = n_float_center;
            } else if (floor(n_float_right) - floor(n_float_center) == 1) {
              k_left = k_center;
              n_float_left = n_float_center;
            } else {
              printf("Bisection error\n");
              return -1;
            }
          }
          // printf("Eigenstate[%d]: k=%.7f\n", k_index, k_center);
          double k = k_center;
          double E = 0.5 * k * k;
          double K = sqrt(2.0 * (E - V0));
          double theta = atan(tan(K * a) * (k / K)) - k * a;
          double alpha = sin(K * a) / sin(k * a + theta);
          ks_unbound_s[k_index] = k;
          Ks_unbound_s[k_index] = K;
          thetas_unbound_s[k_index] = theta;
          alphas_unbound_s[k_index] = alpha;
          norms_unbound_s[k_index] = 2.0 * (a / 2.0 - 1.0 / 4.0 / K * sin(2.0 * K * a) +
                                            alpha * alpha * ((b - a) / 2.0 + 1.0 / 4.0 / k * (sin(2.0 * (k * a + theta)) - sin(2.0 * (k * b + theta)))));
          k_index++;
        }
      }
      char filePath_Psi_unbound_sine[1024];
      sprintf(filePath_Psi_unbound_sine, fileFormat_Psi_unbound_sine, i, j);
      FILE *export_Psi_unbound_sine = fopen(filePath_Psi_unbound_sine, "w");
      fprintf(export_Psi_unbound_sine, "# Unound eigenstates (sine)\n");
      fprintf(export_Psi_unbound_sine, "# index k K theta alpha\n");
      for (l = 0; l < k_count_sine; l++) {
        fprintf(export_Psi_unbound_sine, "%d %.7f %.7f %.7f %.7f\n", l,
                ks_unbound_s[l], Ks_unbound_s[l], thetas_unbound_s[l], alphas_unbound_s[l]);
      }
      fclose(export_Psi_unbound_sine);
      char filePath_Psi_unbound_sine_wfn[1024];
      sprintf(filePath_Psi_unbound_sine_wfn, fileFormat_Psi_unbound_sine_wfn, i, j);
      FILE *export_Psi_unbound_sine_wfn = fopen(filePath_Psi_unbound_sine_wfn, "w");
      fprintf(export_Psi_unbound_sine_wfn, "# Unound eigenstates (sine)\n");
      fprintf(export_Psi_unbound_sine_wfn, "# x ");
      for (l = 0; l < numUnboundStates_export; l++) {
        fprintf(export_Psi_unbound_sine_wfn, "Psi[%d](x) ", l);
      }
      fprintf(export_Psi_unbound_sine_wfn, "\n");
      for (l = 0; l < numPoints; l++) {
        fprintf(export_Psi_unbound_sine_wfn, "%.2f ", x[l]);
        for (m = 0; m < numUnboundStates_export; m++) {
          fprintf(export_Psi_unbound_sine_wfn, "%.7f ",
                  unbound_wfn_sine(x[l], a, b, ks_unbound_s[m], Ks_unbound_s[m], alphas_unbound_s[m], thetas_unbound_s[m], norms_unbound_s[m]));
        }
        fprintf(export_Psi_unbound_sine_wfn, "\n");
      }
      fclose(export_Psi_unbound_sine_wfn);
      if (k_index == k_count_sine) {
        printf("All unbound eigenstates (sine type) are found\n");
      } else {
        printf("Unbound eigenstates (sine type) are missing\n");
        return -1;
      }

      // Excited state calculations
      printf("Step C: Excited state calculations\n");
      for (l = 0; l < numPoints; l++) {
        Psi_ex[l] = complex<double>(0.0, 0.0);
      }
      // cosine type
      for (l = 0; l < k_count_cosine; l++) {
        for (m = 0; m < numPoints; m++) {
          Psi_unbound[m] = unbound_wfn_cosine(x[m], a, b, ks_unbound_c[l], Ks_unbound_c[l], alphas_unbound_c[l], thetas_unbound_c[l], norms_unbound_c[l]);
        }
        double Matrix_element = calc_matrix_element(numPoints, dx, x, Psi_unbound, dPsi_bound, deltaH_amplitude);
        // printf("Index %d: Mat = %.10f %.10f\n", l, Matrix_element.real(), Matrix_element.imag());
        double E = 0.5 * ks_unbound_c[l] * ks_unbound_c[l];
        double E_diff = E - Ekin;
        for (m = 0; m < numPoints; m++) {
          Psi_ex[m] += Matrix_element * (1.0 - complex<double>(cos(E_diff * t), sin(E_diff * t))) / E_diff * complex<double>(cos(E * t), -sin(E * t)) * Psi_unbound[m];
        }
      }
      // sine type
      for (l = 0; l < k_count_sine; l++) {
        for (m = 0; m < numPoints; m++) {
          Psi_unbound[m] = unbound_wfn_sine(x[m], a, b, ks_unbound_s[l], Ks_unbound_s[l], alphas_unbound_s[l], thetas_unbound_s[l], norms_unbound_s[l]);
        }
        double Matrix_element = calc_matrix_element(numPoints, dx, x, Psi_unbound, dPsi_bound, deltaH_amplitude);
        // printf("Index %d: Mat = %.10f %.10f\n", l, Matrix_element.real(), Matrix_element.imag());
        double E = 0.5 * ks_unbound_s[l] * ks_unbound_s[l];
        double E_diff = E - Ekin;
        for (m = 0; m < numPoints; m++) {
          Psi_ex[m] += Matrix_element * (1.0 - complex<double>(cos(E_diff * t), sin(E_diff * t))) / E_diff * complex<double>(cos(E * t), -sin(E * t)) * Psi_unbound[m];
        }
      }
      // amplitude calculation
      double sum_amp = 0.0;
      double sum2_amp = 0.0;
      int numPoints_calculation = 0;
      for (l = 0; l < numPoints; l++) {
        if (x[l] < photoelectron_range_min) {
          continue;
        }
        if (x[l] > photoelectron_range_max) {
          break;
        }
        double Psi_abs = abs(Psi_ex[l]);
        sum_amp += Psi_abs;
        sum2_amp += Psi_abs * Psi_abs;
        numPoints_calculation++;
      }
      double ave_amp = sum_amp / numPoints_calculation;
      double stdev_amp = sqrt(sum2_amp / numPoints_calculation - ave_amp * ave_amp);
      // export
      char filePath_Psi_ex[1024];
      sprintf(filePath_Psi_ex, fileFormat_Psi_ex, i, j);
      FILE *export_Psi_ex = fopen(filePath_Psi_ex, "w");
      fprintf(export_Psi_ex, "# Excited states\n");
      fprintf(export_Psi_ex, "# x Re(Psi_ex(x)) Im(Psi_ex(x))\n");
      for (l = 0; l < numPoints; l++) {
        fprintf(export_Psi_ex, "%.2f %.10f %.10f\n", x[l], Psi_ex[l].real(), Psi_ex[l].imag());
      }
      fclose(export_Psi_ex);

      printf("Step D: Plane-wave matrix element calculation\n");
      // k, Ekin are provided at lines 235-236
      double K = sqrt(2.0 * (Ekin - V0));
      complex<double> Amp_J = (K + k) / 2.0 / K * complex<double>(cos((K - k) * a), -sin((K - k) * a));
      complex<double> Amp_S = (K - k) / 2.0 / K * complex<double>(cos((K + k) * a), sin((K + k) * a));
      complex<double> Amp_I = 1.0 / 4.0 / k / K * ((K + k) * (K + k) * complex<double>(cos(2.0 * (K - k) * a), -sin(2.0 * (K - k) * a)) - (K - k) * (K - k) * complex<double>(cos(2.0 * (K + k) * a), sin(2.0 * (K + k) * a)));
      complex<double> Amp_R = (K * K - k * k) / 4 / k / K * complex<double>(0, 2 * sin(2 * K * a));
      // check
      // printf("|I|^2-|R|^2 = %.3f\n", abs(Amp_I) * abs(Amp_I) - abs(Amp_R) * abs(Amp_R));
      double Amp_norm = (b - a) * (1 + abs(Amp_I) * abs(Amp_I) + abs(Amp_R) * abs(Amp_R)) + 2.0 * a * (abs(Amp_J) * abs(Amp_J) + abs(Amp_S) * abs(Amp_S));
      complex<double> *Psi_PW = alloc_zvector(numPoints);
      for (l = 0; l < numPoints; l++) {
        Psi_PW[l] = unbound_wfn_pw(x[l], a, b, k, K, Amp_J, Amp_S, Amp_I, Amp_R, Amp_norm);
      }
      complex<double> Matrix_element = calc_matrix_element(numPoints, dx, x, Psi_PW, dPsi_bound, deltaH_amplitude);
      // printf("Mat: %.7e\n", abs(Matrix_element));
      double dos = (2.0 * b) / 2 / M_PI;
      double amp_PW = 2.0 * M_PI / k / sqrt(Amp_norm) * dos * abs(Matrix_element);
      // TR-LEED (x<0 side)
      for (l = 0; l < numPoints; l++) {
        Psi_PW[l] = conj(unbound_wfn_pw(x[l], a, b, k, K, Amp_J, Amp_S, Amp_I, Amp_R, Amp_norm));
      }
      Matrix_element = calc_matrix_element(numPoints, dx, x, Psi_PW, dPsi_bound, deltaH_amplitude);
      double amp_PW2 = 2.0 * M_PI / k / sqrt(Amp_norm) * dos * abs(Matrix_element) * abs(Amp_I);
      fprintf(export_Photoelectron_amplitude, "%.3f %.7f %.7f %.7f %.7f\n", Ekin, ave_amp, stdev_amp, amp_PW, amp_PW2);
      // export
      char filePath_Psi_pw_wfn[1024];
      sprintf(filePath_Psi_pw_wfn, fileFormat_Psi_pw_wfn, i, j);
      FILE *export_Psi_pw_wfn = fopen(filePath_Psi_pw_wfn, "w");
      fprintf(export_Psi_pw_wfn, "# PW final states\n");
      fprintf(export_Psi_pw_wfn, "# x Re(Psi_pw(x)) Im(Psi_pw(x))\n");
      for (l = 0; l < numPoints; l++) {
        fprintf(export_Psi_pw_wfn, "%.2f %.10f %.10f\n", x[l], Psi_PW[l].real(), Psi_PW[l].imag());
      }
      fclose(export_Psi_pw_wfn);

      free_dvector(n_float_array);

      free_dvector(ks_unbound_c);
      free_dvector(Ks_unbound_c);
      free_dvector(thetas_unbound_c);
      free_dvector(alphas_unbound_c);
      free_dvector(norms_unbound_c);

      free_dvector(ks_unbound_s);
      free_dvector(Ks_unbound_s);
      free_dvector(thetas_unbound_s);
      free_dvector(alphas_unbound_s);
      free_dvector(norms_unbound_s);
    }
    fclose(export_Photoelectron_amplitude);
  }
}