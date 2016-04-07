package com.company;

import java.util.ArrayList;
import java.util.Vector;

/**
 * Created by admin on 27/11/15.
 */
public class LeastSquares {
    static Double A = -1.0;
    static Double B = 1.0;
    static int POWER = 8;
    static int SIZE_OF_PARTITION = 10000;
    static int COUNT_SOURCE_POINTS = 99;
    static Double [] testPartition = new Double[SIZE_OF_PARTITION];
    static Double [][] Q = new Double[COUNT_SOURCE_POINTS][POWER];
    static Double [] y = new Double[COUNT_SOURCE_POINTS];
    static Double [] x = new Double[COUNT_SOURCE_POINTS];
    static Double [] leastSquaresCoefficients = new Double[COUNT_SOURCE_POINTS];
    static Double [] legandreCoefficients = new Double[COUNT_SOURCE_POINTS];

    static ArrayList<Vector<Point>> data = new ArrayList<Vector<Point>>();
    static Vector<Double> partition = new Vector<Double>();
    static Vector<Point> value_func = new Vector<Point>();
    static Vector<Point> value_leastSquares = new Vector<Point>();
    static Vector<Point> value_legandre = new Vector<Point>();


    public static void start() {
        Double startPoint = A;
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

       Double[][] QTQ = Matrix.multiply(Matrix.transpose(Q), Q);
        Double[] QTy = Matrix.multiply(Matrix.transpose(Q), y);
        leastSquaresCoefficients = Matrix.eliminateWithGauss(QTQ, QTy);
        for (int i = 0; i < POWER; i++) {
            legandreCoefficients[i] = integralFQk(A, B, i) / integralQk2(A, B, new Double(i));
        }
    }

    public static Double basis(Double x, Integer i){
        return Math.pow(x,i);
    }

    public static Double function(Double x){
        return x*Math.tan(x);
    }

    public static Double leastSquareFunction(Double x){
        Double sum = 0.0;
        for (int i =0; i<POWER; i++){
            sum+=leastSquaresCoefficients[i]*basis(x,i);
        }
        return sum;
    }

    public static Double legandreApproximatedFunction(Double x){
        Double sum = 0.0;
        for (int i =0; i<POWER; i++){
            sum+=legandreCoefficients[i]*LegandrePolinomial(x,new Double(i));
        }
        return sum;
    }

    public static Double LegandrePolinomial(Double x, Double i){
        if (i==0.0) return 1.0;
        if (i==1.0) return  x;
        Double n = i - 1 ;
        Double a = ( 2*n + 1 )/(n + 1);
        Double b = n/(n+1);

        return a * x * LegandrePolinomial(x,n) - b * LegandrePolinomial(x,n-1);
    }

    public static Double integralFQk(Double a, Double b, int i){
        Double dx = 1e-6;
        Double xn_1 = a;
        Double sum = 0.0;
        Double xn = a + dx;
        while (xn<=b){
            sum += LegandrePolinomial(xn,new Double(i))*function(xn) * dx;
            xn_1 = xn;
            xn = xn + dx;
        }
        return sum;
    }

    public static Double integralF2(Double a, Double b){
        Double dx = 1e-6;
        Double xn_1 = a;
        Double sum = 0.0;
        Double xn = a + dx;
        while (xn<=b){
            sum += Math.pow(function(xn),2) * dx;
            xn_1 = xn;
            xn = xn + dx;
        }
        return sum;
    }


    public static Double integralQk2(Double a, Double b, Double i){
        Double dx = 1e-6;
        Double xn_1 = a;
        Double sum = 0.0;
        Double xn = a + dx;
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
            Double x = partition.get(i);
            Double value = function(x);
            Double value_leastSquares_y = leastSquareFunction(x);
            Double value_legandre_y = legandreApproximatedFunction(x);
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
        Double start = A;
        for (int i =0; i<SIZE_OF_PARTITION+1; i++){
            partition.add(start);
            start += (B-A)/SIZE_OF_PARTITION;
        }
    }
}
