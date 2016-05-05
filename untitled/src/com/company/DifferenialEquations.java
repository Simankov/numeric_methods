package com.company;
import java.util.ArrayList;
import java.util.Vector;
/**
 * Created by sergey on 21.04.16.
 */
public class DifferenialEquations {
    public double A = 2.0/15.0;
    public double B = 3.0/14.0;
    public Double pi = 3.14;
    public double max_x_k = pi;
    public double delta(){
        return Math.pow(1 / (max_x_k), this.s+1) +
                Math.pow ( Math.sqrt(Math.pow(this.A*this.y_0()[1],2) + Math.pow(-this.B *this.y_0()[1],2)) ,this.s+1 );
    }


    public double x_0 = 0.0;
    public double x_k = pi;
    public double s = 2;
    public Double[] y_0() {return new Double[] {this.B*pi, this.A*pi };}; //todo
    public double epsilon = 0.0000001;
    public double a_epsilon = 0.00001;

    public Vector<Point> y1data = new Vector<Point>();
    public Vector<Point> y2data = new Vector<Point>();
    public Vector<Point> y1exactData = new Vector<Point>();
    public Vector<Point> y2exactData = new Vector<Point>();
    public Vector<Point> y1diff = new Vector<Point>();
    public Vector<Point> y2diff = new Vector<Point>();
    public Vector<Point> ydiff_est = new Vector<Point>();
    public Vector<Point> hs = new Vector<Point>();

     public Double[] y(Double arg){
        return new Double[]{
                (
                        Math.pow(this.B, 3.0 / 2.0) * pi * Math.cos(Math.sqrt(this.A * this.B) * arg) +
                                Math.pow(this.A, 3.0 / 2.0) * pi * Math.sin(Math.sqrt(this.A * this.B) * arg)
                ) / Math.sqrt(B),

                (
                        Math.pow(this.A, 3.0 / 2.0) * pi * Math.cos(Math.sqrt(this.A * this.B) * arg) -
                                Math.pow(this.B, 3.0 / 2.0) * pi * Math.sin(Math.sqrt(this.A * this.B) * arg)
                ) / Math.sqrt(A)
        };

    }

    static public double norm(Double[] x){
        double sum = 0;
        for (int i = 0; i<x.length; i++){
            sum += x[i]*x[i];
        }
        sum = Math.sqrt(sum);
        return sum;
    }
    public double h() {
        return Math.pow(epsilon / delta(), 1.0 / (s + 1));
    }
    public Double[] f(Double x, Double[] y){
        return new Double[] {A*y[1], -B*y[0]};
    }
    public double xsi = 2.0/9.0;
    public Double a21 = xsi;
    public double b2 = 1/(2*xsi);
    public double b1 = 1 - 1/(2*xsi);

    public Double[] k1(double h,double x_i, Double[] y_x_i){
        Double res[] = f(x_i, y_x_i);
        return new Double[] {res[0]*h, res[1]*h};
    }

    public Double[] k2(double h,double x_i, Double[] y_x_i){
        Double arg[] = Matrix.add(y_x_i, Matrix.multiply(k1(h, x_i, y_x_i), a21));
        return Matrix.multiply(f(x_i + xsi*h,arg),h);
    }

    public Double y_x_i(){
        Double estimate = -1.0;
        Double h = h();
        Double x_i = this.x_0;
        Double[] y_i = this.y_0();
        Double[] y_i_h_2 = this.y_0();
        Double[] y_i_h  = this.y_0();

       do {

           for (int i = 1; i <= 2; i++) {
               h = h / i;
               x_i = x_0;
               y_i = y_0();
               y1data = new Vector<Point>();
               y2data = new Vector<Point>();
               y1exactData = new Vector<Point>();
               y2exactData = new Vector<Point>();
               y1diff = new Vector<Point>();
               y2diff = new Vector<Point>();
               while (x_i < x_k) {
                   h = getH(x_i,h);
                   ydiff_est.add(new Point(x_i, Math.abs(getEstimateFail(y_i, x_i, h))));
                   y_i = y_x_i(x_i, y_i, h);
                   x_i = x_i + h;

                   y1data.add(new Point(x_i, y_i[0]));
                   y2data.add(new Point(x_i, y_i[1]));

//                   System.out.println(x_i + "|value|" + (y_i[0] - y(x_i)[0]) + " " + (y_i[1] - y(x_i)[1]));
                   y1exactData.add(new Point(x_i, y(x_i)[0]));
                   y2exactData.add(new Point(x_i, y(x_i)[1]));

                   y1diff.add(new Point(x_i, Math.abs(y_i[0] - y(x_i)[0])));
                   y2diff.add(new Point(x_i, Math.abs(y_i[1] - y(x_i)[1])));
               }
               if (i == 1) {
                   y_i_h = y_i;
               } else {
                   y_i_h_2 = y_i;
               }
           }
       } while (norm( Matrix.subtract(y_i_h,y_i_h_2))/(1.0-Math.pow(2,-this.s)) >= epsilon);
        System.out.println("value in pi:" + norm (Matrix.subtract(y_i_h_2,y(pi))));
        return h;
    }

    public Double[] y_x_i(Double x_i, Double[] y_i, Double h){
        Double[] f1 = Matrix.multiply(k1(h, x_i, y_i), b1);
        Double[] f2 = Matrix.multiply(k2(h, x_i, y_i), b2);
        y_i = Matrix.add(y_i,f1);
        y_i = Matrix.add(y_i,f2);
        return y_i;
    }

    public ArrayList<Vector<Point>> getData(){
        ArrayList<Vector<Point>> data = new ArrayList<Vector<Point>>();
        data.add(y1data);
        data.add(y2data);
        data.add(y1exactData);
        data.add(y2exactData);
        return data;
    }

    public ArrayList<Vector<Point>> getDifference(){
        ArrayList<Vector<Point>> data = new ArrayList<Vector<Point>>();
        data.add(y1diff);
        data.add(y2diff);
        data.add(ydiff_est);
        return data;
    }

    public ArrayList<Vector<Point>> getAutomaticH(){
        ArrayList<Vector<Point>> data = new ArrayList<Vector<Point>>();
        data.add(hs);
        return data;
    }

    public Double[] automaticRunge(){
        double h = y_x_i();
        Double[] y_i = this.y_0();
        double x_i = x_0;
        while (x_i < this.x_k){
            hs.add(new Point(x_i,h));
            Double estimateFail = getEstimateFail(y_i,x_i,h);
            if (estimateFail > a_epsilon * Math.pow(2,s)){
                h = getH(x_i,h/2.0);
                y_i = y_x_i(x_i,y_i,h);
                x_i = x_i + h;
            } else if (estimateFail > a_epsilon && estimateFail <= a_epsilon * Math.pow(2,s)){
                //todo
                h = getH(x_i,h);
                y_i = y_x_i(x_i,y_i,h/2);
                x_i = x_i + h/2; // or h*2
            } else if (estimateFail >= a_epsilon / Math.pow(2,s+1) && estimateFail <= a_epsilon){
                h = getH(x_i,h);
                y_i = y_x_i(x_i,y_i,h);
                x_i = x_i + h;
            } else if (estimateFail < a_epsilon/Math.pow(2,s+1)){
                h = getH(x_i,h*2.0);
                y_i = y_x_i(x_i,y_i,h);
                x_i = x_i + h;
            }
        }

        System.out.println("Automatic:" + norm(Matrix.subtract(y(x_i), y_i)));
        return y_i;
    }

    public Double getEstimateFail(Double[] y_i, Double x_i, Double h){
        Double[] y_h_2 = y_x_i(x_i,y_i,h/2.0);
        Double[] y_h = y_x_i(x_i,y_i,h);
        Double est_fail = norm( Matrix.subtract(y_h,y_h_2))/(1.0-Math.pow(2,-this.s));
        return est_fail;
    }

    Double getH(Double x_i,Double h){
        if (x_i + h > x_k) {
           return x_k - x_i;
        } else {
            return h;
        }
    }


}
