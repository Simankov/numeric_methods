package com.company;

import java.io.Externalizable;
import java.io.PrintWriter;

/**
 * Created by admin on 15/11/15.
 */
public class Minimization {
    static Integer i = 0;
    static Integer a = 5;
    static Integer b = 5;
    static Integer c = 13;
    static double m = 3.62772;
    static Double eps = 1.0E-6;
    static double[][] A = {{4,1,1},{1,2*a,-1},{1,-1,2*b}};
    static double[] B = {1,-2,3};
    static double[] minimum = {-4.0/17.0,37.0/187.0,-48.0/187.0};
    static PrintWriter writer;

    Minimization(){
        try {
            writer = new PrintWriter("output.txt","UTF-8");
        } catch (Exception e){

        }
    }
    public static void gradientMethod(){
        writer.print("Градиентный метод начался:");
        double[] xprev = {0,0,0};
        double[] next = Matrix.add( xprev ,  Matrix.multiply(grad(xprev),u(xprev)) );
        double [] diff = new double [3];
        while (Matrix.euclidNorm(grad(next))/m >= eps){
            double[] temp = next;
            next = Matrix.add( xprev ,  Matrix.multiply(grad(xprev),u(xprev)) );
            xprev = temp;
            i++;
            writer.print("Шаг номер:");
            writer.println(i);
            writer.printf("Различие величины нормы градиента функции и eps: %.8f\n", Matrix.euclidNorm(grad(next)) - eps);
            diff = Matrix.subtract(minimum, next);
            writer.printf("Следующая точка: (%.8f,%.8f,%.8f)\n", next[0], next[1], next[2]);
            writer.printf("Разница между настоящим минимумом и текущей точкой: (%.8f,%.8f,%.8f)\n", diff[0], diff[1], diff[2]);


        }
        writer.printf("Минимум найден в точке: (%.8f,%.8f,%.8f)\n", next[0], next[1], next[2]);
        writer.printf("Разница между настоящим минимумом и найденным: (%.8f,%.8f,%.8f)\n", diff[0], diff[1], diff[2]);
        writer.printf("Значение нормы градиента функции: %.8f\n", Matrix.euclidNorm(grad(next)));
    }

    public static void coordinateMethod() {
        writer.print("Координатный метод начался:");
        double[] xprev = {0, 0, 0};
        double[] e1 = {1, 0, 0};
        int j = 0;
        double[] diff = new double[3];
        double[] next = Matrix.add(xprev, Matrix.multiply(e1, u_coord(xprev, ei(0))));
        while (Matrix.euclidNorm(grad(next))/m >= eps) {
            j++;
            double[] temp = next;
            next = Matrix.add(xprev, Matrix.multiply(ei(j), u_coord(xprev, ei(j))));
            xprev = temp;
            writer.print("Шаг номер:");
            writer.println(j);
            writer.printf("Различие величины нормы градиента функции и eps: %.8f\n", Matrix.euclidNorm(grad(next)) - eps);
            diff = Matrix.subtract(minimum, next);
            writer.printf("Следующая точка: (%.8f,%.8f,%.8f)\n", next[0], next[1], next[2]);
            writer.printf("Разница между настоящим минимумом и текущей точкой: (%.8f,%.8f,%.8f)\n", diff[0], diff[1], diff[2]);

        }

        writer.printf("Минимум найден в точке: (%.8f,%.8f,%.8f)\n", next[0], next[1], next[2]);
        writer.printf("Разница между настоящим минимумом и найденным: (%.8f,%.8f,%.8f)\n", diff[0], diff[1], diff[2]);
        writer.printf("Значение нормы градиента функции: %.8f\n\n\n\n\n\n\n\n"
, Matrix.euclidNorm(grad(next)));
    }

    static double[] ei(Integer i){
        Integer remainder = i%3;
        if (remainder == 0) return new double[] {1,0,0};
        if (remainder == 1) return new double[] {0,1,0};
        if (remainder == 2) return new double[] {0,0,1};
        return new double[] {1,0,0};
    }

    static double function(double[] x) {
        double[] xTA = Matrix.multiply(x, A);
        double xTAx = Matrix.dot(xTA, x);
        double xTb = Matrix.dot(x, B);
        return 0.5 * xTAx + xTb + c;
    }
    static void stop(){
        writer.close();
    }

    static double u(double [] x){
        double[] grad = grad(x);
        double gradientNorm = Matrix.euclidNorm(grad);
        double [] gradTA = Matrix.multiply(grad, A);
        double gradTAgrad = Matrix.dot(gradTA,grad);
        double uk =  -Math.pow(gradientNorm,2) / (gradTAgrad);
        return uk;
    }


    static double u_coord(double [] x, double[] ei){
        double[] grad = grad(x);
        double eiTgrad = Matrix.dot(ei,grad);
        double[] eiTA = Matrix.multiply(ei, A);
        double eiTAei = Matrix.dot(eiTA, ei);
        double uk =  -eiTgrad/(eiTAei);
        return uk;
    }

    static double [] grad(double [] x){
        double [] Ax = Matrix.multiply(A,x);
        return  Matrix.add(Ax,B);
    }
}
