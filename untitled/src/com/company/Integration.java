package com.company;

import java.util.HashMap;

public class Integration {
    public static Double aa = 1.5;
    public static Double bb = 2.3;
    public static Double alpha = 1.0/5.0;
    public static Double realValue = 25.0102;
    public static Double realValueWithP = 32.2195;
    public Double[] difference = new Double[5];
    public HashMap<String, Double> SH_i = new HashMap<String, Double>();

    public static Double f(Double arg){
        return 2*Math.cos(3.5*arg)*Math.exp(5*arg/3)+3*Math.sin(1.5*arg)*Math.exp(-4*arg)+3;
    }

    public Double p(Double arg){
        return 1/Math.pow(arg - aa, alpha);
    }

    public Double firstQuadrature(Double a, Double b){
        return (b-a)*f(a);
    }

    public Double secondQuadrature(Double a, Double b){
        return (b-a)*f((a + b) / 2.0);
    }

    public Double thirdQuadrature(Double a, Double b){
        return (b-a)/2.0*(f(a)+f(b));
    }

    public Double fourthQuadrature(Double a, Double b){
        return (b-a)/6.0*(f(a)+4*f((a+b)/2.0)+f(b));
    }

    public double riemannIntegral(Double n){
        Double sum = 0.0;
        Double x_k = aa;
        Double x_k_1 = aa;
        for (int j = 0; j<n; j++){
            x_k_1 = x_k;
            x_k += (bb-aa)/n;
            sum += f((x_k+x_k_1)/2)*(bb-aa)/n;
        }
        System.out.println(realValue - sum);
        return sum;
    }

    public double NewtonKotes(){
        return NewtonKotes(1,aa,bb);
    }

    public double smartIntegration(Double n){
        for (int type = 0; type<=4; type++) {
            Double a_k = aa;
            Double b_k = aa;
            Double sum = 0.0;
            for (int i = 0; i < n; i++) {
                a_k = b_k;
                b_k += (bb - aa) / n;
                switch (type) {
                    case 1:
                        sum += firstQuadrature(a_k, b_k);
                        break;
                    case 2:
                        sum += secondQuadrature(a_k, b_k);
                        break;
                    case 3:
                        sum += thirdQuadrature(a_k, b_k);
                        break;
                    case 4:
                        sum += fourthQuadrature(a_k, b_k);
                        break;
                }
            }
            difference[type] = Math.abs(realValue - sum);
        }
        return 0.0;
    }




    public double NewtonKotes(int n,double a,double b) {
            double a_k = a;
            double b_k = a;
            double sum = 0;
            for (int i = 0; i < n; i++) {
                a_k = b_k;
                b_k += (b - a) / n;
                sum += NewtonKotes(a_k, b_k);
            }
        return sum;
    }

    public double NewtonKotes(double a_k, double b_k){
        double x1 = a_k;
        double x2 = (a_k+b_k)/2;
        double x3 = b_k;
        double [][] A = new double[][] {{1,1,1},{x1,x2,x3},{x1*x1,x2*x2,x3*x3}};
        double [] b = new double[] {nu(a_k,b_k,0),nu(a_k,b_k,1),nu(a_k,b_k,2)};
        double [] Ai = Matrix.eliminateWithGauss(A,b);
        double sum = 0;
        sum += f(x1)*Ai[0] +f(x2)*Ai[1] + f(x3)*Ai[2] ;
        return sum;
    }

    public double[] solveEquation(double a, double b, double c, double d){
        double q = (2*b*b*b / (27*a*a*a) - b*c / (3*a*a) + d / a ) ;
        double p = (3*a*c - b*b)/(3*a*a);
        double D = q*q/4 + p*p*p/9;
        assert (D>0);
        assert (p<0);
        double r = q>0 ? Math.sqrt(Math.abs(p)) : -Math.sqrt(Math.abs(p));
        double arg = q / (r*r*r);
        double phi = Math.acos(arg);
        double y1 = -2*r*Math.cos(phi/3.0);
        double y2 = 2*r*Math.cos(3.1415926/3.0 - phi/3.0);
        double y3 = 2*r*Math.cos(3.1415926/3.0 + phi/3.0);

       return new double[] {y1 - b/(3*a),y2 - b/(3*a),y3 - b/(3*a)};
    }

    public double Gauss(int n,int nsplit, double a,double b){

        double sum = 0;
        double a_k = a;
        double b_k = a;
        for (int k=0; k<nsplit; k++) {
            a_k = b_k;
            b_k += (b - a) / nsplit;
            double[][] Q = new double[n][n];
            double[] b1 = new double[n];
            for (int i = 0; i < n; i++) {
                b1[i] = (-1) * nu(aa,bb, n + i);
                for (int j = 0; j < n; j++) {
                    Q[i][j] = nu(aa,bb, i + (n - 1 - j));
                }
            }

            double[] ai = Matrix.eliminateWithGauss(Q, b1);
            double[] xi = solveEquation(1, ai[0], ai[1], ai[2]);
            double[][] A = new double[][]{{1, 1, 1}, {xi[0], xi[1], xi[2]}, {xi[0] * xi[0], xi[1] * xi[1], xi[2] * xi[2]}};
            double[] bb = new double[]{nu(aa,this.bb,0), nu(aa,this.bb,1), nu(aa,this.bb,2)};
            double[] Ai = Matrix.eliminateWithGauss(A, bb);
            double result = 0;
            for (int i = 0; i < n; i++) {
                result += Ai[i] * f(xi[i]);
            }
            sum += result;
        }
        return sum;
    }
    public double nu(double a_k,double b_k,int i){
        return nu(i,b_k) - nu(i,a_k);
    }


    public double nu(int i,double arg){

        switch(i){
            case 0: return Math.pow(arg - aa,1-alpha) / (1-alpha);
            case 1: return Math.pow(arg - aa,1-alpha) * (aa-alpha*arg + arg) / ((alpha-2)*(alpha-1));
            case 2: return -Math.pow(arg - aa,1-alpha) * (2*aa*aa - 2*aa*(alpha - 1)*arg + (alpha*alpha - 3*alpha + 2)*arg*arg) / ((alpha - 3)*(alpha - 2)*(alpha - 1));
            case 3: return Math.pow(arg - aa,1-alpha) * (6*Math.pow(aa,3) - 6*Math.pow(aa,2)*(alpha-1)*arg +
                            3*aa*(Math.pow(alpha,2) - 3*alpha + 2)*Math.pow(arg,2) -
                            (Math.pow(alpha,3) - 6*Math.pow(alpha,2) + 11*alpha - 6)*Math.pow(arg,3)) /
                            ((alpha-4)*(alpha - 3)*(alpha - 2)*(alpha - 1));
            case 4:
                return -Math.pow(arg - aa,1-alpha) *
            (24 * Math.pow(aa,4)-24*Math.pow(aa,3) *  (alpha-1)*arg+12*Math.pow(aa,2)*(Math.pow(alpha,2)-3*alpha+2)*Math.pow(arg,2)
                    -4*aa*(Math.pow(alpha,3)-6*Math.pow(alpha,2)+11*alpha-6)*Math.pow(arg,3)+(Math.pow(alpha,4)
                    -10*Math.pow(alpha,3)+35*Math.pow(alpha,2)-50*alpha+24)*Math.pow(arg,4))
                        /
                        ((alpha - 5)*(alpha-4)*(alpha - 3)*(alpha - 2)*(alpha - 1));
            case 5:
                return Math.pow(arg - aa,1-alpha) * (120*Math.pow(aa,5)-120*Math.pow(aa,4)*(alpha-1)*arg+60*Math.pow(aa,3)* (Math.pow(alpha,2)-3*alpha+2)*Math.pow(arg,2)-20*aa*aa*(alpha*alpha*alpha-6*alpha*alpha+11*alpha-6)*Math.pow(arg,3)+5*aa*(Math.pow(alpha,4)-10*Math.pow(alpha,3)+35*Math.pow(alpha,2)-50*alpha+24)*Math.pow(arg,4)
                        -(Math.pow(alpha,5)-15*Math.pow(alpha,4)+85*Math.pow(alpha,3)-225*Math.pow(alpha,2)+274*alpha-120)*Math.pow(arg,5))
                / ((alpha - 6)*(alpha - 5)*(alpha-4)*(alpha - 3)*(alpha - 2)*(alpha - 1));
        }
        return 0;
    }

    public double AitkenProcess(double eps){
        double q = 2.0;
        double i = 0;
        double i_1 = 0;
        double m_i_1 = 0;
        double m_i = 0;
        double i_2 = 0;
        do {
            i = i_1;
            i_1 = i+1;
            i_2 = i_1+1;
            m_i = m_i_1;
            m_i_1 = Math.log ( Math.abs ( (  S(i,q) - S(i_1,q) ) / (S(i_1,q) - S(i_2,q)) ) ) / Math.log(q);
        } while (Math.abs(m_i_1 - m_i)>=eps);
        return m_i_1;
 }


 public double S(Double pow,Double q){
     String key = pow.toString() + " " + q.toString();
     if (SH_i.containsKey(key)){
         return SH_i.get(key);
     }
     double i = Math.pow(q,pow);
     double H = (bb-aa)/i;
     double a_k = aa;
    double b_k = aa;
    double sum = 0;
    for (int j = 0; j<i; j++) {
        a_k = b_k;
        b_k += H;
//       sum+=f((a_k+b_k)/2.0)*H;
          sum+=NewtonKotes(a_k,b_k);
    }
     SH_i.put(key,sum);
    return sum;
 }


  public double RungeMethod(double eps1,double eps2){
      double i=-1;
      double R = 0;
      double q=2;
      double result = 0;
      double m = AitkenProcess(eps2);
      do {
          double i_1 = i+1;
          double i_2 = i+2;
          double SH_1 = S(i_1,q);
          double SH_2 = S(i_2,q);
          double Hi_1_m = Math.pow((bb-aa)/Math.pow(q,i_1),m);
          double Hi_2_m = Math.pow((bb-aa)/Math.pow(q,i_2),m);
          double Y = (SH_1/Hi_1_m - SH_2/Math.pow((i_2*q),m)) / (1/Hi_1_m - 1/Hi_2_m);
          double Cm = ( SH_2 - SH_1 ) /  ( Hi_1_m - Hi_2_m );
          R = Cm * Hi_2_m;
          result = SH_2;
          if (Math.abs(R)<eps1) {
              System.out.println("Runge:"+i_2);
          }
          i++;
      }
      while (Math.abs(R)>=eps1);
      return result;
  }

    public double RichardsonMethod(int r,double eps1,double eps2){
        double q = 2;
        double m = AitkenProcess(eps2);
        double Hi[] = new double[r];
        for (int j=0; j<r; j++){
            Hi[j] = (bb-aa)/Math.pow(q,j);
        }
        double [][] Q = new double[r][r];
        double[] b = new double[r];
        for (int i= 0; i<r; i++){
            for (int j=0;j<r; j++){
                Q[i][j] = (j==0)?1:- Math.pow(Hi[i],m+j);
            }
            b[i]=S((double)i,q);
        }
        double[] x = Matrix.eliminateWithGauss(Q,b);
        double R = 0;
        for (int i = 1; i<r; i++){
            R+=x[i]*Math.pow(Hi[r-1],m+i);
        }
        if (Math.abs(R)<eps1){
            System.out.println("Richardson:"+(r-1));
            return S((double)r-1,q);
        } else {
            return RichardsonMethod(r+1,eps1,eps2);
        }
    }





    public double RichardsonMethod(double eps1,double eps2){
        return RichardsonMethod(2,eps1,eps2);
    }

}
