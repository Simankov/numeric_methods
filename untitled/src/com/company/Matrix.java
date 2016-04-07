package com.company;

/**
 * Created by admin on 15/11/15.
 */
public class Matrix {


    // return x^T y
    public static Double dot(Double[] x, Double[] y) {
        Double sum = 0.0;
        for (int i = 0; i < x.length; i++)
            sum += x[i] * y[i];
        return sum;
    }

    // return C = A^T
    public static Double[][] transpose(Double[][] A) {
        int m = A.length;
        int n = A[0].length;
        Double[][] C = new Double[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[j][i] = A[i][j];
        return C;
    }

    // return C = A + B
    public static Double[][] add(Double[][] A, Double[][] B) {
        int m = A.length;
        int n = A[0].length;
        Double[][] C = new Double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] = A[i][j] + B[i][j];
        return C;
    }

    // return C = A + B
    public static Double[] add(Double[] A, Double[] B) {
        int n = A.length;
        Double[] C = new Double[n];
        for (int i = 0; i < n; i++)
                C[i] = A[i] + B[i];
        return C;
    }

    // return C = A - B
    public static Double[][] subtract(Double[][] A, Double[][] B) {
        int m = A.length;
        int n = A[0].length;
        Double[][] C = new Double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] = A[i][j] - B[i][j];
        return C;
    }

    public static Double[] subtract(Double[] A, Double[] B) {
        int n = A.length;
        Double[] C = new Double[n];
            for (int i = 0; i < n; i++)
                C[i] = A[i] - B[i];
        return C;
    }

    public static Double euclidNorm(Double [] x){
        int n = x.length;
        Double result = 0.0;
        for (int i = 0; i<n; i++){
            result += Math.pow(x[i],2);
        }
        return Math.sqrt(result);
    }

    // return C = A * B
    public static Double[][] multiply(Double[][] A, Double[][] B) {
        int mA = A.length;
        int nA = A[0].length;
        int mB = B.length;
        int nB = B[0].length;
        Double[][] C = new Double[mA][nB];
        for (int i = 0; i < mA; i++)
            for (int j = 0; j < nB; j++)
                for (int k = 0; k < nA; k++)
                    C[i][j] += A[i][k] * B[k][j];
        return C;
    }

    // matrix-vector multiplication (y = A * x)
    public static Double[] multiply(Double[][] A, Double[] x) {
        int m = A.length;
        int n = A[0].length;
        Double[] y = new Double[m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                y[i] += A[i][j] * x[j];
        return y;
    }

    public static Double[] multiply(Double[] x, Double a) {
        int n = x.length;
        Double[] C = new Double[n];
        for (int i = 0; i < n; i++) {
            C[i] = x[i] * a;
        }
        return C;
    }

    // vector-matrix multiplication (y = x^T A)
    public static Double[] multiply(Double[] x, Double[][] A) {
        int m = A.length;
        int n = A[0].length;
        Double[] y = new Double[n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                y[j] += A[i][j] * x[i];
        return y;
    }

    // Gaussian elimination with partial pivoting
    public static Double[] eliminateWithGauss(Double[][] A, Double[] b) {
        int N = b.length;

        for (int p = 0; p < N; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            Double[] temp = A[p];
            A[p] = A[max];
            A[max] = temp;
            Double t = b[p];
            b[p] = b[max];
            b[max] = t;

            // pivot within A and b
            for (int i = p + 1; i < N; i++) {
                Double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < N; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // back substitution
        Double[] x = new Double[N];
        for (int i = N - 1; i >= 0; i--) {
            Double sum = 0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }

}