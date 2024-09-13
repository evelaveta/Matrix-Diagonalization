#include <stdio.h>
#include <stdlib.h>
#include "nr.h"
#include "nrutil.h"
    
void printm (float **M, int size)
{
    for (int i = 1; i <= size; i++)
    {
        for (int j = 1; j <= size; j++)
        {
           printf(" %10f ",M[i][j]);
        }
        printf( "\n");
    }
    printf("\n");
}
int main(void)
{
    float **S,**v, **S_copy, **res;
    float *e, *r, *d;
    int cj = 0, cqr = 0;
    int n;

    scanf("%d",&n);
    
    float **S2;
    S = dmatrix(1, n, 1, n);
    S_copy = dmatrix(1, n, 1, n);
    S2 = dmatrix(1, n, 1, n);
    v = dmatrix(1, n, 1, n);
    res = dmatrix(1, n, 1, n);
    e = dvector(1, n);
    d = dvector(1, n);
    r = dvector(1, n);
    
    for (int i = 1; i <= n; i++)
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
                S[i][j] = (int)rand() % 100000 /500;
                S[j][i] = S[i][j];
                S_copy[i][j] = S[i][j];
                S2[i][j] = S[i][j];
        }
    }
    
    printf("Матрица: \n");
    printm(S, n);
    
    jacobi(S, n, d, v, &cj);

     for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= i; j++)
        {
                res[i][j] = 0;
                
                if(i == j){
                    res[i][j] = d[i];
                }
        }
    }

    printf("СВ Якоби: \n");
    printm(v, n);
    
    for (int i = 1; i <= n; i++)
    {
        printf("СЗ Якоби %d: %f \n", i, d[i]);
    }
    printf("\n");
    printf("\n");
    
    tred2(S_copy, n, r, e);
    tqli(r, e, n, S_copy, &cqr);
    
    printf("СВ QR: \n");
    printm(S_copy, n);
    
    for (int i = 1; i <= n; i++)
    {
        printf("СЗ QR %d: %f \n", i, r[i]);
    }
    printf("\n");
    
    printf("Якоби счетчик= %d\n", cj);
    printf("QR счетчик= %d\n", cqr);
    
    printm(res, n);
    return 0;
}
