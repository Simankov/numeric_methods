package com.company;


public class Main {
    public static void main(String[] args) {
        double epsilon = 0.001;
        double epsilonForM = 0.01;

        DifferenialEquations differenialEquations = new DifferenialEquations();
        System.out.println("finded" + differenialEquations.y_x_i());
        differenialEquations.automaticRunge();

        Plotter plotter = new Plotter("differentialEq",differenialEquations.getData(),new String[]{"y1","y2","exy1","exy2"});
        plotter.plot();

        Plotter plotter2 = new Plotter("differentialEq",differenialEquations.getDifference(),new String[]{"y1diff","y2diff","est"});
        plotter2.plot();
        Plotter plotter3 = new Plotter("differentialEq",differenialEquations.getAutomaticH(),new String[]{"hs"});
        plotter3.plot();


    }
};