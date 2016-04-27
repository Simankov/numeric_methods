package com.company;


public class Main {
    public static void main(String[] args) {
        double epsilon = 0.00001;
        double epsilonForM = 0.01;

        for (Integration.Type type :
             Integration.Type.values()) {
            Integration integration1 = new Integration();
            integration1.type = type;
            System.out.println(integration1.AitkenProcess(epsilonForM));
            System.out.println(integration1.RungeMethod(epsilon, epsilonForM));
            System.out.println(integration1.RichardsonMethod(epsilon, epsilonForM));
        }
    }
};