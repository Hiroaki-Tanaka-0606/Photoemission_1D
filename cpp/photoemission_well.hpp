#define USE_MATH_DEFINES
#include <cmath>
#include <complex>
using namespace std;

double **alloc_dmatrix(int N) {
  double **matrix = new double *[N];
  double *alloc = new double[N * N];
  for (int i = 0; i < N; i++) {
    matrix[i] = &alloc[i * N];
  }
  return matrix;
}
complex<double> **alloc_zmatrix(int N) {
  complex<double> **matrix = new complex<double> *[N];
  complex<double> *alloc = new complex<double>[N * N];
  for (int i = 0; i < N; i++) {
    matrix[i] = &alloc[i * N];
  }
  return matrix;
}
double **alloc_dmatrix(int N, int M) {
  double **matrix = new double *[N];
  double *alloc = new double[N * M];
  for (int i = 0; i < N; i++) {
    matrix[i] = &alloc[i * M];
  }
  return matrix;
}
double *alloc_dvector(int N) {
  return new double[N];
}
complex<double> *alloc_zvector(int N) {
  return new complex<double>[N];
}
void free_dmatrix(double **mat) {
  delete[] mat[0];
  delete[] mat;
}
void free_zmatrix(complex<double> **mat) {
  delete[] mat[0];
  delete[] mat;
}
void free_dvector(double *vec) {
  delete[] vec;
}
void free_zvector(complex<double> *vec) {
  delete[] vec;
}

double diffFunc_bound_cosine(double K, double V0, double a) {
  double E = 0.5 * K * K + V0;
  double kappa = sqrt(-2.0 * E);
  return tan(K * a) - kappa / K;
}

double diffFunc_bound_sine(double K, double V0, double a) {
  double E = 0.5 * K * K + V0;
  double kappa = sqrt(-2.0 * E);
  return tan(K * a) + K / kappa;
}

double n_float_cosine(double k, double V0, double a, double b) {
  double E = 0.5 * k * k;
  double K = sqrt(2.0 * (E - V0));
  double theta = atan(1.0 / tan(K * a) * (-k / K)) - k * a;
  double kbpt = k * b + theta;
  return kbpt / M_PI;
}

double n_float_sine(double k, double V0, double a, double b) {
  double E = 0.5 * k * k;
  double K = sqrt(2.0 * (E - V0));
  double theta = atan(tan(K * a) * (k / K)) - k * a;
  double kbpt = k * b + theta;
  return kbpt / M_PI;
}

double unbound_wfn_cosine(double x, double a, double b, double k, double K, double alpha, double theta, double norm) {
  if (x < -a) {
    return alpha * sin(-k * x + theta) / sqrt(norm);
  } else if (x < a) {
    return cos(K * x) / sqrt(norm);
  } else {
    return alpha * sin(k * x + theta) / sqrt(norm);
  }
}

double unbound_wfn_sine(double x, double a, double b, double k, double K, double alpha, double theta, double norm) {
  if (x < -a) {
    return -alpha * sin(-k * x + theta) / sqrt(norm);
  } else if (x < a) {
    return sin(K * x) / sqrt(norm);
  } else {
    return alpha * sin(k * x + theta) / sqrt(norm);
  }
}

complex<double> unbound_wfn_pw(double x, double a, double b, double k, double K, complex<double> J, complex<double> S, complex<double> I, complex<double> R, double norm) {
  if (x < -a) {
    return (I * complex<double>(cos(k * x), sin(k * x)) + R * complex<double>(cos(k * x), -sin(k * x))) / sqrt(norm);
  } else if (x < a) {
    return (J * complex<double>(cos(K * x), sin(K * x)) + S * complex<double>(cos(K * x), -sin(K * x))) / sqrt(norm);
  } else {
    return complex<double>(cos(k * x), sin(k * x)) / sqrt(norm);
  }
}

// <Psi1| dH d/dx |Psi2>
double calc_matrix_element(int N, double dx, double *x, double *Psi1, double *dPsi2, double deltaH) {
  double sum = 0.0;
  int i;
  for (i = 0; i < N - 1; i++) {
    sum += Psi1[i] * dPsi2[i] * dx;
  }
  return sum * deltaH;
}

// <Psi1| dH d/dx |Psi2>
complex<double> calc_matrix_element(int N, double dx, double *x, complex<double> *Psi1, double *dPsi2, double deltaH) {
  complex<double> sum = 0.0;
  // double norm = 0.0;
  int i;
  for (i = 0; i < N - 1; i++) {
    sum += Psi1[i] * dPsi2[i] * dx;
    // norm += abs(Psi1[i]) * abs(Psi1[i]) * dx;
  }
  // printf("Norm: %.3f\n", norm);
  return sum * deltaH;
}
