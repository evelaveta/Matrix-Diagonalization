#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "nr.h"
#include "nrutil.h"
#define N 3
void printmatrix(float **M, int size, FILE *F);

void printmatrix(float **M, int size, FILE *F)
{
    for (int i = 1; i <= size; i++)
    {
        for (int j = 1; j <= size; j++)
        {
           fprintf(F, " %10f ",M[i][j]);
        }
        fprintf(F, "\n");
    }
    fprintf(F, "\n");
}
int main(void)
{
    FILE *F;
    float **S,**v, **S_copy;
    float *e, *r, *d;
    int counter_jacobi = 0, counter_QR = 0;
    int n;

    if((F = fopen("output.txt", "wb")) == NULL)
    {
        printf("Can`t find file\n");
        exit(EXIT_FAILURE);
    }
    //printf("Enter matrix size: ");
    //scanf("%d",&n);
    n = N;
    S = dmatrix(1, n, 1, n);
    S_copy = dmatrix(1, n, 1, n);
    v = dmatrix(1, n, 1, n);
    e = dvector(1, n);
    d = dvector(1, n);
    r = dvector(1, n);
    for(int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
           S[i][j] = 0;
           S_copy[i][j] = 0;
        }
    }
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= i; j++)
        {
                S[i][j] = (double)rand()/RAND_MAX;
                S[j][i] = S[i][j];
                S_copy[i][j] = S[i][j];
        }
    }
    fprintf(F, "Symmetric Matrix: \n");
    printmatrix(S, n, F);
    jacobi(S, n, d, v, &counter_jacobi);
    fprintf(F, "Eigenvectors Matrix(Jacobi): \n");
    printmatrix(v, n, F);
    for (int i = 1; i <= n; i++)
    {
        fprintf(F, "Eigenvalue(Jacobi) %d: %f \n", i, d[i]);
    }
    fprintf(F, "\n");
    fprintf(F, "\n");
    tred2(S_copy, n, r, e);
    tqli(r, e, n, S_copy, &counter_QR);
    fprintf(F, "Eigenvectors Matrix(OR): \n");
    printmatrix(S_copy, n, F);
    for (int i = 1; i <= n; i++)
    {
        fprintf(F, "Eigenvalue(QR) %d: %f \n", i, r[i]);
    }
    fprintf(F, "\n");
    fprintf(F, "Jacobi_counter = %d\n",counter_jacobi);
    fprintf(F, "QR_counter = %d\n",counter_QR);
    fclose(F);
    return 0;
}
