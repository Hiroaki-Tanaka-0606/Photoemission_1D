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

extern "C" {
void dsyev_(
    char *JOBZ,
    char *UPLO,
    int *N,
    double *A,
    int *LDA,
    double *W,
    double *WORK,
    int *LWORK,
    int *INFO);
void zgemv_(char *TRANS,
            int *M,
            int *N,
            complex<double> *ALPHA,
            complex<double> *A,
            int *LDA,
            complex<double> *X,
            int *INCX,
            complex<double> *BETA,
            complex<double> *Y,
            int *INCY);
}