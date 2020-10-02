#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>

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

void rand_identity_matrix(int *matrix, int n){
    //int a = rand()%n;
    int a = 1;
    memset(matrix, 0, n*sizeof(int));
    matrix[a] = 1;
}

void rank_one_update(double *A, double *B, double *u, int *v, int n){
    int j = 0;
    while(v[j] < 1) j++;
    for (int i = 0; i < n; i++){
        B[i*n + j] = A[i*n + j] + u[i];
    }
}

double sign(double x){
    return (x > 0) - (x < 0);
}

void identity(double *A, int n) {
    memset(A, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        A[i * n + i] = 1;
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

void hessenberg_qr(double *Q, double *R, int n) {
    double c = 0, s = 0, r = 0, new_row = 0, new_col = 0;
    identity(Q, n);

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

void w_calc(double *Q, double *u, double *w, int n) {
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

void mult_mat(double *A, double *B, int n) {
    double sum = 0;
    //double C[n][n];
    double C[n];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum = 0;
            for (int k = 0; k < n; k++) {
                sum = sum + A[i * n + k] * B[k * n + j];
            }
            C[j] = sum;
            //C[i][j] = sum;
        }
        memcpy(&A[i*n], C, n * sizeof(double));
    }

    //MELHORAR ESSA COPIA
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            A[i * n + j] = C[i][j];
//        }
//    }
}

void qr_update(double *Q, double *R, double *u, int *v, int n) {
    int cont = 0;
    double c = 0, s = 0, r = 0, new_row = 0, new_col = 0;
    double w[n], Q1[n][n];

    w_calc(Q, u, w, n);

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
    hessenberg_qr(&Q1[0][0], R, n);
    mult_mat(Q, &Q1[0][0], n);
}

void qr_method(double *dataA, double *dataQ, double *dataR, int n) {
    int i, j, k, m;
    double sum = 0;

    for (i = 0; i < n; i++) {
        sum = 0;
        for (m = 0; m < n; m++) {
            sum += dataA[m * n + i] * dataA[m * n + i];
        }
        dataR[i * n + i] = sqrt(sum);

        for (k = 0; k < n; k++) {
            dataQ[k * n + i] = dataA[k * n + i] / dataR[i * n + i];
        }

        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = 0; k < n; k++) {
                sum += dataQ[k * n + i] * dataA[k * n + j];
            }
            dataR[i * n + j] = sum;

            for (k = 0; k < n; k++) {
                dataA[k * n + j] = dataA[k * n + j] - (dataR[i * n + j] * dataQ[k * n + i]);
            }
        }
    }
}

int main() {
//    double dataA[n * n];// = {2, 4, -2, 4, 9, -3, -2, -3, 7};
//    double dataB[n * n];
//    double u[n] = {0, 0, 2};
//    int v[n] = {0, 1, 0};
//    double dataQ[n * n], dataQ1[n * n];
//    double dataR[n * n], dataR1[n * n];

    double *dataA, *dataB, *dataQ, *dataR, *dataQ1, *dataR1, *u;
    struct timeval before, after;
    long int accumTimeUpdate = 0, accumTimeQr = 0;
    int *v, interval = 100, n = 2, num_commands = 8;

    FILE * temp = fopen("dataUpdate.temp", "w");
    FILE * temp2 = fopen("dataQR.temp", "w");
    FILE * gnuplotPipe = popen ("gnuplot -persist", "w");
    char * commandsForGnuplot[] = {
            "set xrange[2:100]",
            "set yrange[0:2000]",
            "set style line 1 linecolor rgb \'#0060ad\' linetype 1 linewidth 2",
            "set style line 2 linecolor rgb \'#dd181f\' linetype 1 linewidth 2",
            "set title \"QR Decomposition vs QR Decomposition Update\"",
            "set xlabel \"Matrix size\"",
            "set ylabel \"Computation time\"",
            "plot 'dataUpdate.temp' with lines linestyle 1, 'dataQR.temp' with lines linestyle 2"
    };
//    "plot '-' with lines \n"
//    for (int c = 0; c < num_commands; c++){
//        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[c]);
//    }

    dataA = (double *) malloc(n*n *sizeof(double));
    dataB = (double *) malloc(n*n *sizeof(double));
    dataQ = (double *) malloc(n*n *sizeof(double));
    dataQ1 = (double *) malloc(n*n *sizeof(double));
    dataR = (double *) malloc(n*n *sizeof(double));
    dataR1 = (double *) malloc(n*n *sizeof(double));
    u = (double *) malloc(n *sizeof(double));
    v = (int *) malloc(n *sizeof(int));

//    dataA[0] = 2;
//    dataA[1] = 4;
//    dataA[2] = -2;
//    dataA[3] = 4;
//    dataA[4] = 9;
//    dataA[5] = -3;
//    dataA[6] = -2;
//    dataA[7] = -3;
//    dataA[8] = 7;
//
//    u[0] = 0;
//    u[1] = 0;
//    u[2] = 2;
//
//    v[0] = 0;
//    v[1] = 1;
//    v[2] = 0;

    for (int x = 0; x < interval; x++) {
        dataA = (double *) realloc(dataA, n*n *sizeof(double));
        dataB = (double *) realloc(dataB, n*n *sizeof(double));
        dataQ = (double *) realloc(dataQ, n*n *sizeof(double));
        dataQ1 = (double *) realloc(dataQ1, n*n *sizeof(double));
        dataR = (double *) realloc(dataR, n*n *sizeof(double));
        dataR1 = (double *) realloc(dataR1, n*n *sizeof(double));
        u = (double *) realloc(u, n *sizeof(double));
        v = (int *) realloc(v, n *sizeof(int));

        rand_matrix(dataA, n, n);

        memcpy(dataQ, dataA, n*n* sizeof(double));
        memset(dataR, 0, n*n* sizeof(double));
        memcpy(dataB, dataA, n*n* sizeof(double));

        qr_method(dataA, dataQ, dataR, n);
        rand_matrix(u, 1, n);
        rand_identity_matrix(v, n);

        rank_one_update(dataA, dataB, u, v, n);
        memcpy(dataQ1, dataB, n*n* sizeof(double));
        memset(dataR1, 0, n*n* sizeof(double));

        //update
        gettimeofday(&before, NULL);
        qr_update(dataQ, dataR, u, v, n);
        gettimeofday(&after, NULL);
        accumTimeUpdate = ((after.tv_sec * 1000000 + after.tv_usec) - (before.tv_sec * 1000000 + before.tv_usec));
        fprintf(temp, "%d %g\n", x, (double) accumTimeUpdate);

        gettimeofday(&before, NULL);
        qr_method(dataB, dataQ1, dataR1, n);
        gettimeofday(&after, NULL);
        accumTimeQr = ((after.tv_sec * 1000000 + after.tv_usec) - (before.tv_sec * 1000000 + before.tv_usec));
        fprintf(temp2, "%d %g\n", x, (double) accumTimeQr);
        n++;
    }
    fclose(temp);
    fclose(temp2);

    for (int c = 0; c < num_commands; c++){
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[c]);
    }
    fflush(gnuplotPipe);
    free(dataA);
    free(dataB);
    free(dataQ);
    free(dataQ1);
    free(dataR);
    free(dataR1);
    free(u);
    free(v);

    pclose(gnuplotPipe);

    return 0;
}