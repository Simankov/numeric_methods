package com.company;

/**
 * Created by sergey on 21.04.16.
 */
public class DifferenialEquations {
    public double A = 3;
    public double B = 4;
    public double max_x_k = 3.14;
    public double delta = Math.pow(1 / ( this.x_0 > this.max_x_k ? this.x_0 : this.x_k ), this.s+1) +
            Math.pow ( norm(f(0.0, new double[]{this.B*this.x_k,this.A*this.x_k} )) ,this.s+1 );

    public double x_0 = 0;
    public double x_k = 3.14;
    public double s = 2;
    public double epsilon = 0.0001;
    public double norm(double[] x){
        double sum = 0;
        for (int i = 0; i<x.length; i++){
            sum += x[i]*x[i];
        }
        sum = Math.sqrt(sum);
        return sum;
    }
    public double h = Math.pow(epsilon / delta, 1/ (s+1));
    public double[] f(double x, double[] y){
        return new double[] {A*y[1], -B*y[0]};
    }





}
