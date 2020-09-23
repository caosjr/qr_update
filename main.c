#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>

#define TAM 47

void print_matrix(double *matrix, int m, int n){
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++) {
            printf("%lf\t", matrix[i*m + j]);
        }
        printf("\n");
    }
}

void rand_matrix(double *matrix, int m, int n){
    double a = 2;
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++) {
            matrix[i*m + j] = ((float)rand()/(float)(RAND_MAX)) * a;
        }
    }
}

void rand_identity_matrix(int *matrix){
    int a = rand()%TAM;
    memset(matrix, 0, TAM*sizeof(int));
    matrix[a] = 1;
}

void rank_one_update(double *A, double *B, double *u, int *v){
    int j = 0;
    while(v[j] < 1) j++;
    for (int i = 0; i < TAM; i++){
        B[i*TAM + j] = A[i*TAM + j] + u[i]; //B = A + u*v'
    }
}

double sign(double x){
    return (x > 0) - (x < 0);
}

void identity(double *A) {
    memset(A, 0, TAM * TAM * sizeof(double));
    for (int i = 0; i < TAM; i++) {
        A[i * TAM + i] = 1;
    }
}

void givens_rotation(double a, double b, double *c, double *s, double *r) {
    double t = 0, u = 0;
    *c = 0;
    *s = 0;
    *r = 0;
    if (b == 0) {
        *c = sign(a);
        *s = 0;
        *r = fabs(a);
    } else if (a == 0) {
        *c = 0;
        *s = -sign(b);
        *r = fabs(b);
    } else if (fabs(a) > fabs(b)) {
        t = b / a;
        u = sign(a) * fabs(sqrt(1 + pow(t, 2)));
        *c = 1 / u;
        *s = -*c * t;
        *r = a * u;
    } else {
        t = a / b;
        u = sign(b) * fabs(sqrt(1 + pow(t, 2)));
        *s = -1 / u;
        *c = -*s * t;
        *r = b * u;
    }
}

void hessenberg_qr(double *Q, double *R) {
    int n = TAM;
    double c = 0, s = 0, r = 0, new_row = 0, new_col = 0;
    identity(Q);

    for (int k = 0; k <= n - 2; k++) {
        givens_rotation(R[k * n + k], R[(k + 1) * n + k], &c, &s, &r);
        for (int j = 0; j < n; j++) {
            new_row = c * R[k * n + j] - s * R[(k + 1) * n + j];
            R[(k + 1) * n + j] = s * R[k * n + j] + c * R[(k + 1) * n + j];
            R[k * n + j] = new_row;
            new_col = c * Q[j * n + k] - s * Q[j * n + (k + 1)];
            Q[j * n + (k + 1)] = s * Q[j * n + k] + c * Q[j * n + (k + 1)];
            Q[j * n + k] = new_col;
        }
    }
}

void w_calc(double *Q, double *u, double *w) {
    int n = TAM;
    double sum;

    for (int i = 0; i < n; i++) {
        sum = 0;
        for (int k = 0; k < n; k++) {
            //transposed Q
            sum += Q[k * n + i] * u[k]; //first[i][k]*second[k][j];
        }
        w[i] = sum;
    }
}

void mult_mat(double *A, double *B) {
    int n = TAM;
    double sum = 0;
    double C[n][n];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum = 0;
            for (int k = 0; k < n; k++) {
                sum = sum + A[i * TAM + k] * B[k * TAM + j];
            }
            C[i][j] = sum;
        }
    }

    //MELHORAR ESSA COPIA
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = C[i][j];
        }
    }
}

void qr_update(double *Q, double *R, double *u, int *v) {
    int n = TAM, cont = 0;
    double c = 0, s = 0, r = 0, new_row = 0, new_col = 0;
    double w[n], Q1[n][n];

    w_calc(Q, u, w);

    for (int k = n - 2; k >= 0; k--) {
        givens_rotation(w[k], w[k + 1], &c, &s, &r);
        w[k + 1] = 0;
        w[k] = r;

        for (int j = 0; j < n; j++) {
            new_row = c * R[k * n + j] - s * R[(k + 1) * n + j];
            R[(k + 1) * n + j] = s * R[k * n + j] + c * R[(k + 1) * n + j];
            R[k * n + j] = new_row;
            new_col = c * Q[j * n + k] - s * Q[j * n + (k + 1)];
            Q[j * n + (k + 1)] = s * Q[j * n + k] + c * Q[j * n + (k + 1)];
            Q[j * n + k] = new_col;
        }
    }
    cont = 0;
    while (v[cont] != 1) cont++; //gambiarra
    R[0 * n + cont] += w[0];
    hessenberg_qr(&Q1[0][0], R);
    mult_mat(Q, &Q1[0][0]);
}

void qr_method(double *dataA, double *dataQ, double *dataR) {
    int i, j, k, m;
    double sum = 0;

    for (i = 0; i < TAM; i++) {
        sum = 0;
        for (m = 0; m < TAM; m++) {
            sum += dataA[m * TAM + i] * dataA[m * TAM + i];
        }
        dataR[i * TAM + i] = sqrt(sum);

        for (k = 0; k < TAM; k++) {
            dataQ[k * TAM + i] = dataA[k * TAM + i] / dataR[i * TAM + i];
        }

        for (j = i + 1; j < TAM; j++) {
            sum = 0;
            for (k = 0; k < TAM; k++) {
                sum += dataQ[k * TAM + i] * dataA[k * TAM + j];
            }
            dataR[i * TAM + j] = sum;

            for (k = 0; k < TAM; k++) {
                dataA[k * TAM + j] = dataA[k * TAM + j] - (dataR[i * TAM + j] * dataQ[k * TAM + i]);
            }
        }
    }
}

int main() {
    double dataA[TAM * TAM];// = {2, 4, -2, 4, 9, -3, -2, -3, 7};
    double dataB[TAM * TAM];
    double u[TAM] = {0, 0, 2};
    int v[TAM] = {0, 1, 0};
    double dataQ[TAM * TAM], dataQ1[TAM * TAM];
    double dataR[TAM * TAM], dataR1[TAM * TAM];
    struct timeval before, after;
    long int accumTimeUpdate = 0, accumTimeQr = 0;
    int interval = 10000;

    rand_matrix(dataA, TAM, TAM);
//    printf("\nMatriz A:\n");
//    print_matrix(dataA, TAM, TAM);
    memcpy(dataQ, dataA, sizeof(dataA));
    memset(dataR, 0, sizeof(dataR));
    memcpy(dataB, dataA, sizeof(dataA));

    qr_method(dataA, dataQ, dataR);

    for (int n = 0; n < interval; n++) {
        memcpy(dataQ1, dataB, sizeof(dataB));
        memset(dataR1, 0, sizeof(dataR1));

        rand_matrix(u, 1, TAM);
//        printf("\nVetor u:\n");
//        print_matrix(u, 1, TAM);

//        printf("\nVetor v:\n");
        rand_identity_matrix(v);

        rank_one_update(dataA, dataB, u, v);
//        printf("\nMatriz B:\n");
//        print_matrix(dataB, TAM, TAM);


//        printf("\n\nMatriz Q\n");
//        print_matrix(dataQ, TAM, TAM);
//        printf("\n\nMatriz R\n");
//        print_matrix(dataR, TAM, TAM);

        //update
        gettimeofday(&before, NULL);
        qr_update(dataQ, dataR, u, v);
        gettimeofday(&after, NULL);
        accumTimeUpdate += ((after.tv_sec * 1000000 + after.tv_usec) - (before.tv_sec * 1000000 + before.tv_usec));

        gettimeofday(&before, NULL);
        qr_method(dataB, dataQ1, dataR1);
        gettimeofday(&after, NULL);
        accumTimeQr += ((after.tv_sec * 1000000 + after.tv_usec) - (before.tv_sec * 1000000 + before.tv_usec));

//        printf("\n\nUpdate Matriz Q\n");
//        print_matrix(dataQ, TAM, TAM);
//        printf("\n\nUpdate Matriz R\n");
//        print_matrix(dataR, TAM, TAM);
    }
    printf("\nQR update decomposition B: %ld", accumTimeUpdate);
    printf("\nQR decomposition B: %ld", accumTimeQr);

    return 0;
}