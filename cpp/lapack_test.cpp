#include <cstdio>
#include <chrono>

using namespace std;

extern "C" {
    void dgemm_(
        char* transa,
        char* transb,
        int* m,
        int* n,
        int* k,
        double* alpha,
        double* a,
        int* lda,
        double* b,
        int* ldb,
        double* beta,
        double* c,
        int* ldc
    );
}


double** alloc_matrix(int N){
    double** matrix=new double*[N];
    double* alloc=new double[N*N];
    for(int i=0; i<N; i++){
        matrix[i]=&alloc[i*N];
    }
    return matrix;
}
void free_matrix(double** mat){
    delete[] mat[0];
    delete[] mat;
}

int main(){
    int N=100;
    int T=1000;

    double** A=alloc_matrix(N);
    double** B=alloc_matrix(N);
    double** C=alloc_matrix(N);

    int i,j;
    int t;
    
    char notrans='N';
    double alpha=1.0;
    double beta=0.0;

    chrono::system_clock::time_point start=chrono::system_clock::now();
    for(t=0; t<T; t++){
        for(i=0; i<N; i++){
            for(j=0; j<N; j++){
                A[i][j]=i+j;
                B[i][j]=i-j;
            }
        }
        dgemm_(&notrans, &notrans, &N, &N, &N, &alpha, &A[0][0], &N, &B[0][0], &N, &beta, &C[0][0], &N);
    }
    chrono::system_clock::time_point end=chrono::system_clock::now();
    double duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();


    printf("Dgemm %d times: %.0f [ms]\n", T, duration);
    free_matrix(A);
    free_matrix(B);
}