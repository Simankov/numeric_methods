package com.company;

import java.util.ArrayList;
import java.util.Vector;

/**
 * Created by admin on 27/11/15.
 */
public class LeastSquares {
    static double A = -1;
    static double B = 1;
    static int POWER = 8;
    static int SIZE_OF_PARTITION = 10000;
    static int COUNT_SOURCE_POINTS = 99;
    static double [] testPartition = new double[SIZE_OF_PARTITION];
    static double [][] Q = new double[COUNT_SOURCE_POINTS][POWER];
    static double [] y = new double[COUNT_SOURCE_POINTS];
    static double [] x = new double[COUNT_SOURCE_POINTS];
    static double [] leastSquaresCoefficients = new double[COUNT_SOURCE_POINTS];
    static double [] legandreCoefficients = new double[COUNT_SOURCE_POINTS];

    static ArrayList<Vector<Point>> data = new ArrayList<Vector<Point>>();
    static Vector<Double> partition = new Vector<Double>();
    static Vector<Point> value_func = new Vector<Point>();
    static Vector<Point> value_leastSquares = new Vector<Point>();
    static Vector<Point> value_legandre = new Vector<Point>();


    public static void start() {
        double startPoint = A;
        for (int i = 0; i < COUNT_SOURCE_POINTS; i++) {
            x[i] = startPoint;
            y[i] = function(startPoint);
            startPoint += (B - A) / 4;
        }

        for (int i = 0; i < Q.length; i++) {
            for (int j = 0; j < Q[i].length; j++) {
                Q[i][j] = basis(x[i], j);
            }
        }

        double[][] QTQ = Matrix.multiply(Matrix.transpose(Q), Q);
        double[] QTy = Matrix.multiply(Matrix.transpose(Q), y);
        leastSquaresCoefficients = Matrix.eliminateWithGauss(QTQ, QTy);
        for (int i = 0; i < POWER; i++) {
            legandreCoefficients[i] = integralFQk(A, B, i) / integralQk2(A, B, i);
        }
    }

    public static double basis(Double x, Integer i){
        return Math.pow(x,i);
    }

    public static double function(double x){
        return x*Math.tan(x);
    }

    public static double leastSquareFunction(double x){
        double sum = 0;
        for (int i =0; i<POWER; i++){
            sum+=leastSquaresCoefficients[i]*basis(x,i);
        }
        return sum;
    }

    public static double legandreApproximatedFunction(double x){
        double sum = 0;
        for (int i =0; i<POWER; i++){
            sum+=legandreCoefficients[i]*LegandrePolinomial(x,i);
        }
        return sum;
    }

    public static double LegandrePolinomial(double x, double i){
        if (i==0.0) return 1.0;
        if (i==1.0) return  x;
        double n = i - 1 ;
        double a = ( 2*n + 1 )/(n + 1);
        double b = n/(n+1);

        return a * x * LegandrePolinomial(x,n) - b * LegandrePolinomial(x,n-1);
    }

    public static double integralFQk(double a, double b, int i){
        double dx = 1e-6;
        double xn_1 = a;
        double sum = 0;
        double xn = a + dx;
        while (xn<=b){
            sum += LegandrePolinomial(xn,i)*function(xn) * dx;
            xn_1 = xn;
            xn = xn + dx;
        }
        return sum;
    }

    public static double integralF2(double a, double b){
        double dx = 1e-6;
        double xn_1 = a;
        double sum = 0;
        double xn = a + dx;
        while (xn<=b){
            sum += Math.pow(function(xn),2) * dx;
            xn_1 = xn;
            xn = xn + dx;
        }
        return sum;
    }


    public static double integralQk2(double a, double b, double i){
        double dx = 1e-6;
        double xn_1 = a;
        double sum = 0;
        double xn = a + dx;
        while (xn<=b){
            sum += Math.pow(LegandrePolinomial(xn,i),2) * dx;
            xn_1 = xn;
            xn = xn + dx;
        }
        return sum;
    }


    static void compute(){
        fillPartition();
        for (int i=0; i<SIZE_OF_PARTITION; i++) {
            double x = partition.get(i);
            double value = function(x);
            double value_leastSquares_y = leastSquareFunction(x);
            double value_legandre_y = legandreApproximatedFunction(x);
            value_func.add(new Point(x, value));
            value_leastSquares.add(new Point(x, value_leastSquares_y));
            value_legandre.add(new Point(x,value_legandre_y));
        }
        data.add(value_leastSquares);
        data.add(value_func);
        data.add(value_legandre);
    }

    static public ArrayList<Vector<Point>> getData() {
        compute();
        return data ;
    }

    static void fillPartition(){
        double start = A;
        for (int i =0; i<SIZE_OF_PARTITION+1; i++){
            partition.add(start);
            start += (B-A)/SIZE_OF_PARTITION;
        }
    }
}
