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
    static Double m = 3.62772;
    static Double eps = 1.0E-6;
    static Double[][] A = {{4.0,1.0,1.0},{1.0,new Double(2*a),-1.0},{1.0,-1.0,new Double(2*b)}};
    static Double[] B = {1.0,-2.0,3.0};
    static Double[] minimum = {-4.0/17.0,37.0/187.0,-48.0/187.0};
    static PrintWriter writer;

    Minimization(){
        try {
            writer = new PrintWriter("output.txt","UTF-8");
        } catch (Exception e){

        }
    }
    public static void gradientMethod(){
        writer.print("Градиентный метод начался:");
        Double[] xprev = {0.0,0.0,0.0};
        Double[] next = Matrix.add( xprev ,  Matrix.multiply(grad(xprev),u(xprev)) );
        Double [] diff = new Double [3];
        while (Matrix.euclidNorm(grad(next))/m >= eps){
            Double[] temp = next;
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
        Double[] xprev = {0.0, 0.0, 0.0};
        Double[] e1 = {1.0, 0.0, 0.0};
        int j = 0;
        Double[] diff = new Double[3];
        Double[] next = Matrix.add(xprev, Matrix.multiply(e1, u_coord(xprev, ei(0))));
        while (Matrix.euclidNorm(grad(next))/m >= eps) {
            j++;
            Double[] temp = next;
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

    static Double[] ei(Integer i){
        Integer remainder = i%3;
        if (remainder == 0) return new Double[] {1.0,0.0,0.0};
        if (remainder == 1) return new Double[] {0.0,1.0,0.0};
        if (remainder == 2) return new Double[] {0.0,0.0,1.0};
        return new Double[] {1.0,0.0,0.0};
    }

    static Double function(Double[] x) {
        Double[] xTA = Matrix.multiply(x, A);
        Double xTAx = Matrix.dot(xTA, x);
        Double xTb = Matrix.dot(x, B);
        return 0.5 * xTAx + xTb + c;
    }
    static void stop(){
        writer.close();
    }

    static Double u(Double [] x){
        Double[] grad = grad(x);
        Double gradientNorm = Matrix.euclidNorm(grad);
        Double [] gradTA = Matrix.multiply(grad, A);
        Double gradTAgrad = Matrix.dot(gradTA,grad);
        Double uk =  -Math.pow(gradientNorm,2) / (gradTAgrad);
        return uk;
    }


    static Double u_coord(Double [] x, Double[] ei){
        Double[] grad = grad(x);
        Double eiTgrad = Matrix.dot(ei,grad);
        Double[] eiTA = Matrix.multiply(ei, A);
        Double eiTAei = Matrix.dot(eiTA, ei);
        Double uk =  -eiTgrad/(eiTAei);
        return uk;
    }

    static Double [] grad(Double [] x){
        Double [] Ax = Matrix.multiply(A,x);
        return  Matrix.add(Ax,B);
    }
}
