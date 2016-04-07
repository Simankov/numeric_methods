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

    public Double riemannIntegral(Double n){
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

    public Double NewtonKotes(){
        return NewtonKotes(1,aa,bb);
    }

    public Double smartIntegration(Double n){
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




    public Double NewtonKotes(int n,Double a,Double b) {
            Double a_k = a;
            Double b_k = a;
            Double sum = 0.0;
            for (int i = 0; i < n; i++) {
                a_k = b_k;
                b_k += (b - a) / n;
                sum += NewtonKotes(a_k, b_k);
            }
        return sum;
    }

    public Double NewtonKotes(Double a_k, Double b_k){
        Double x1 = a_k;
        Double x2 = (a_k+b_k)/2;
        Double x3 = b_k;
        Double [][] A = new Double[][] {{1.0,1.0,1.0},{x1,x2,x3},{x1*x1,x2*x2,x3*x3}};
        Double [] b = new Double[] {nu(a_k,b_k,0),nu(a_k,b_k,1),nu(a_k,b_k,2)};
        Double [] Ai = Matrix.eliminateWithGauss(A,b);
        Double sum = 0.0;
        sum += f(x1)*Ai[0] +f(x2)*Ai[1] + f(x3)*Ai[2] ;
        return sum;
    }

    public Double[] solveEquation(Double a, Double b, Double c, Double d){
        Double q = (2*b*b*b / (27*a*a*a) - b*c / (3*a*a) + d / a ) ;
        Double p = (3*a*c - b*b)/(3*a*a);
        Double D = q*q/4 + p*p*p/9;
        assert (D>0);
        assert (p<0);
        Double r = q>0 ? Math.sqrt(Math.abs(p)) : -Math.sqrt(Math.abs(p));
        Double arg = q / (r*r*r);
        Double phi = Math.acos(arg);
        Double y1 = -2*r*Math.cos(phi/3.0);
        Double y2 = 2*r*Math.cos(3.1415926/3.0 - phi/3.0);
        Double y3 = 2*r*Math.cos(3.1415926/3.0 + phi/3.0);

       return new Double[] {y1 - b/(3*a),y2 - b/(3*a),y3 - b/(3*a)};
    }

    public Double Gauss(int n,int nsplit, Double a,Double b){

        Double sum = 0.0;
        Double a_k = a;
        Double b_k = a;
        for (int k=0; k<nsplit; k++) {
            a_k = b_k;
            b_k += (b - a) / nsplit;
            Double[][] Q = new Double[n][n];
            Double[] b1 = new Double[n];
            for (int i = 0; i < n; i++) {
                b1[i] = (-1) * nu(a,b, n + i);
                for (int j = 0; j < n; j++) {
                    Q[i][j] = nu(a,b, i + (n - 1 - j));
                }
            }

            Double[] ai = Matrix.eliminateWithGauss(Q, b1);
            Double[] xi = solveEquation(1.0, ai[0], ai[1], ai[2]);
            Double[][] A = new Double[][]{{1.0, 1.0, 1.0}, {xi[0], xi[1], xi[2]}, {xi[0] * xi[0], xi[1] * xi[1], xi[2] * xi[2]}};
            Double[] bb = new Double[]{nu(a,b,0), nu(a,b,1), nu(a,b,2)};
            Double[] Ai = Matrix.eliminateWithGauss(A, bb);
            Double result = 0.0;
            for (int i = 0; i < n; i++) {
                result += Ai[i] * f(xi[i]);
            }
            sum += result;
        }
        return sum;
    }
    public Double nu(Double a_k,Double b_k,int i){
        return nu(i,b_k) - nu(i,a_k);
    }


    public Double nu(int i,Double arg){

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
        return 0.0;
    }

    public Double AitkenProcess(Double eps){
        Double q = 2.0;
        Double i = 0.0;
        Double i_1 = 0.0;
        Double m_i_1 = 0.0;
        Double m_i = 0.0;
        Double i_2 = 0.0;
        do {
            i = i_1;
            i_1 = i+1;
            i_2 = i_1+1;
            m_i = m_i_1;
            m_i_1 = Math.log ( Math.abs ( (  S(i,q) - S(i_1,q) ) / (S(i_1,q) - S(i_2,q)) ) ) / Math.log(q);
        } while (Math.abs(m_i_1 - m_i)>=eps);
        return m_i_1;
 }


 public Double S(Double pow,Double q){
     String key = pow.toString() + " " + q.toString();
     if (SH_i.containsKey(key)){
         return SH_i.get(key);
     }
     Double i = Math.pow(q,pow);
     Double H = (bb-aa)/i;
     Double a_k = aa;
    Double b_k = aa;
    Double sum = 0.0;
    for (int j = 0; j<i; j++) {
        a_k = b_k;
        b_k += H;
        sum += Gauss(3,1,a_k,b_k);
    }
     SH_i.put(key,sum);
    return sum;
 }


  public Double RungeMethod(Double eps1,Double eps2){
      Double i=-1.0;
      Double R = 0.0;
      Double q=2.0;
      Double result = 0.0;
      Double m = AitkenProcess(eps2);
      do {
          Double i_1 = i+1;
          Double i_2 = i+2;
          Double SH_1 = S(i_1,q);
          Double SH_2 = S(i_2,q);
          Double Hi_1_m = Math.pow((bb-aa)/Math.pow(q,i_1),m);
          Double Hi_2_m = Math.pow((bb-aa)/Math.pow(q,i_2),m);
          Double Y = (SH_1/Hi_1_m - SH_2/Math.pow((i_2*q),m)) / (1/Hi_1_m - 1/Hi_2_m);
          Double Cm = ( SH_2 - SH_1 ) /  ( Hi_1_m - Hi_2_m );
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

    public Double RichardsonMethod(int r, Double eps1, Double eps2){
        Double q = 2.0;
        Double m = AitkenProcess(eps2);
        Double Hi[] = new Double[r];
        for (int j=0; j<r; j++){
            Hi[j] = (bb-aa)/Math.pow(q,j);
        }
        Double [][] Q = new Double[r][r];
        Double[] b = new Double[r];
        for (int i= 0; i<r; i++){
            for (int j=0;j<r; j++){
                Q[i][j] = (j==0)?1:- Math.pow(Hi[i],m+j);
            }
            b[i]=S(new Double(i),q);
        }
        Double[] x = Matrix.eliminateWithGauss(Q,b);
        Double R = 0.0;
        for (int i = 1; i<r; i++){
            R+=x[i]*Math.pow(Hi[r-1],m+i);
        }
        if (Math.abs(R)<eps1){
            System.out.println("Richardson:"+(r-1));
            return S(new Double(r-1),q);
        } else {
            return RichardsonMethod(r+1,eps1,eps2);
        }
    }





    public Double RichardsonMethod(Double eps1,Double eps2){
        return RichardsonMethod(2,eps1,eps2);
    }

}
