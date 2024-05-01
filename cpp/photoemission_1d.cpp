// Before execution, please check that you properly prepared the constants at the beginning of main().

#define USE_MATH_DEFINES
#include "photoemission_1d.hpp"
#include <cmath>
#include <complex>
#include <cstdio>

using namespace std;
int main() {
  // Constants
  /// System
  double period_slab = 4.0;
  double potential_slab_center = -1.15;
  double potential_slab_amplitude = 0.85;
  double vacuum_length = 100;
  int num_layers = 9;
  int numPoints_per_unit = 10;
  /// Photoemission simulations
  double Ekin_min = .5;
  double Ekin_delta = 0.1;
  double numEkin = 1;
  int numLEStates_sim = 1;
  double deltaH_amplitude = 0.001;
  double propagation_length = 30;
  double lower_integration_range = 0.5;
  double upper_integration_range = 2.0;
  /// Export of the results
  const char *filePath_V = "./output/Potential.dat";
  const char *filePath_PsiLE = "./output/Psi_lowEnergy.dat";
  int numLEStates_export = 10;
  const char *fileFormat_Psi_ex = "./output/Psi_ex_%d_%d.dat";
  /// Physical constants
  const double fine_structure_constant = 7.297e-3;
  // utility variables
  // k is skipped because k is used as the wavevector variable
  int i, j, l, m;

  // Step 1: prepare the potential
  // x range: [-vacuum_length, period_slab*num_layers+vacuum_length]
  printf("Step 1: Preparation of the slab potential\n");
  double slab_length = period_slab * num_layers;
  double system_length = vacuum_length * 2.0 + slab_length;
  double dx = 1.0 / numPoints_per_unit;
  int numPoints = round(system_length * numPoints_per_unit);
  printf("System length: %.1f (Bohr), Number of points: %d\n", system_length, numPoints);

  double *V = alloc_dvector(numPoints);
  double *x = alloc_dvector(numPoints);

  for (i = 0; i < numPoints; i++) {
    x[i] = dx * i - vacuum_length;
    if (-period_slab * 0.25 <= x[i] && x[i] <= slab_length + period_slab * 0.25) {
      V[i] = potential_slab_center - potential_slab_amplitude * cos(2.0 * M_PI * x[i] / period_slab);
    } else {
      double y_scale = 2.0 * (-potential_slab_center);
      double x_scale = 4.0 * M_PI * potential_slab_amplitude / period_slab / (-potential_slab_center);
      if (x[i] < 0.0) {
        double x_offset = x[i] + period_slab * 0.25;
        V[i] = potential_slab_center * 2 + y_scale / (exp(x_offset * x_scale) + 1.0);
      } else {
        double x_offset = x[i] - slab_length - period_slab * 0.25;
        V[i] = potential_slab_center * 2 + y_scale / (exp(-x_offset * x_scale) + 1.0);
      }
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

  // Step 2: Obtain eigenstate(s)
  printf("Step 2: Eigenstates\n");
  // Hmat is inverted
  // Eigenvectors are vertical (i-th eigenvector is in the (*, i) components -> [i][*])
  double **Hmat = alloc_dmatrix(numPoints);
  for (i = 0; i < numPoints; i++) {
    for (j = 0; j < numPoints; j++) {
      if (i == j) {
        Hmat[i][j] = V[i] + 1.0 / dx / dx;
      } else if (i == j - 1 || i == j + 1) {
        Hmat[i][j] = -0.5 / dx / dx;
      } else {
        Hmat[i][j] = 0.0;
      }
    }
  }
  char JOBZ = 'V';
  char UPLO = 'U';
  double *W = alloc_dvector(numPoints);
  double *WORK;
  double WORK_dummy;
  int LWORK = -1;
  int INFO;
  dsyev_(&JOBZ, &UPLO, &numPoints, &Hmat[0][0], &numPoints, &W[0], &WORK_dummy, &LWORK, &INFO);
  if (INFO == 0) {
    LWORK = round(WORK_dummy);
    printf("WORK size: %d\n", LWORK);
    WORK = alloc_dvector(LWORK);
  } else {
    printf("Error on dsyev_ test\n");
    return -1;
  }
  dsyev_(&JOBZ, &UPLO, &numPoints, &Hmat[0][0], &numPoints, &W[0], &WORK[0], &LWORK, &INFO);
  if (INFO == 0) {
    printf("dsyev_ succeeded\n");
  } else {
    printf("Error on dsyev_\n");
    return -1;
  }
  // Psis is NOT inverted
  double **Psis = alloc_dmatrix(numPoints);
  for (i = 0; i < numPoints; i++) {
    // printf("E[%2d] = %8.3f\n", i, W[i]);
    // Normalization
    double norm = 0.0;
    for (j = 0; j < numPoints; j++) {
      norm += Hmat[i][j] * Hmat[i][j] * dx;
    }
    double coef = 1.0 / sqrt(norm);
    // printf("coef: %.3f\n", coef);
    for (j = 0; j < numPoints; j++) {
      Psis[i][j] = Hmat[i][j] * coef;
    }
  }
  FILE *export_PsiLE = fopen(filePath_PsiLE, "w");
  fprintf(export_PsiLE, "# Low-energy eigenstates\n");
  fprintf(export_PsiLE, "# Eigenvalues:\n");
  fprintf(export_PsiLE, "# ");
  for (i = 0; i < numLEStates_export; i++) {
    fprintf(export_PsiLE, "%.3f ", W[i]);
  }
  fprintf(export_PsiLE, "\n");
  fprintf(export_PsiLE, "# x ");
  for (i = 0; i < numLEStates_export; i++) {
    fprintf(export_PsiLE, "Psi[%d](x) ", i);
  }
  fprintf(export_PsiLE, "\n");
  for (i = 0; i < numPoints; i++) {
    fprintf(export_PsiLE, "%.2f ", x[i]);
    for (j = 0; j < numLEStates_export; j++) {
      fprintf(export_PsiLE, "%.5f ", Psis[j][i]);
    }
    fprintf(export_PsiLE, "\n");
  }
  fclose(export_PsiLE);
  printf("The low-energy eigenstates were exported to %s\n", filePath_PsiLE);
  free_dvector(WORK);
  free_dmatrix(Hmat);

  // Step 3: Direct photoemission simulations
  printf("Step 3: Direct photoemission simulations\n");
  complex<double> *Psi_ex = alloc_zvector(numPoints);
  char TRANS = 'N';
  complex<double> ALPHA(1.0, 0.0);
  int INC = 1;
  complex<double> BETA(0.0, 0.0);
  double *dPsi_bound = alloc_dvector(numPoints);
  for (i = 0; i < numLEStates_sim; i++) {
    printf("Simulation for eigenstate #%d\n", i);
    // differenciate
    for (j = 0; j < numPoints; j++) {
      if (j == 0) {
        dPsi_bound[j] = (Psis[i][j + 1] - Psis[i][j]) / dx;
      } else if (j == numPoints - 1) {
        dPsi_bound[j] = (Psis[i][j] - Psis[i][j - 1]) / dx;
      } else {
        dPsi_bound[j] = (Psis[i][j + 1] - Psis[i][j - 1]) / dx / 2;
      }
    }

    for (j = 0; j < numEkin; j++) {
      double Ekin = Ekin_min + Ekin_delta * j;
      double k = sqrt(2.0 * Ekin);
      printf("Photoelectron kinetic energy: %.3f\n", Ekin);
      printf("Photoelectron wavevector: %.3f\n", k);
      double hn = Ekin - W[i];
      double kph = hn * fine_structure_constant;
      printf("Photon energy: %.3f\n", hn);
      printf("Photon momentum: %.3f\n", kph);
      double t = propagation_length / k;
      printf("T_max: %.3f\n", t);
      // Excited part calculation based on the perturbation theory
      for (l = 0; l < numPoints; l++) {
        Psi_ex[l] = complex<double>(0.0, 0.0);
      }
      for (l = 0; l < numPoints; l++) {
        if (W[l] < lower_integration_range * Ekin) {
          continue;
        }
        if (W[l] > upper_integration_range * Ekin) {
          break;
        }
        printf("Eigenstate[%d]: eigenvalue=%.3f\n", l, W[l]);
        complex<double> Matrix_element(0.0, 0.0);
        for (m = 0; m < numPoints; m++) {
          Matrix_element += Psis[l][m] * deltaH_amplitude * dPsi_bound[m] * dx;
        }
        double E_diff = W[l] - Ekin;
        complex<double> cnt = (1.0 - complex<double>(cos(E_diff * t), sin(E_diff * t))) / E_diff;
        for (m = 0; m < numPoints; m++) {
          Psi_ex[m] += cnt * Psis[l][m] * complex<double>(cos(W[l] * t), -sin(W[l] * t));
        }
      }
      // For debug: export the wave function
      char filePath_Psi_ex[1024];
      sprintf(filePath_Psi_ex, fileFormat_Psi_ex, i, j);
      FILE *export_Psi_ex = fopen(filePath_Psi_ex, "w");
      fprintf(export_Psi_ex, "# Excited part of the wavefunction\n");
      fprintf(export_Psi_ex, "# x Re(Psi(t)) Im(Psi(t))\n");
      for (l = 0; l < numPoints; l++) {
        fprintf(export_Psi_ex, "%.2f %.5f %.5f\n", x[l], Psi_ex[l].real(), Psi_ex[l].imag());
      }
    }
  }

  free_dmatrix(Psis);
  free_dvector(W);
}