package com.company;

import java.util.ArrayList;
import java.util.Vector;



public class Main {
    public static void main(String[] args) {
        Integration integration = new Integration();
        Vector<Point> data = new Vector<Point>();
        ArrayList<Vector<Point>> plotData = new ArrayList<Vector<Point>>();


//        Vector<Point> firstPoints = new Vector<Point>();
//        Vector<Point> secondPoints = new Vector<Point>();
//        Vector<Point> thirdPoints = new Vector<Point>();
//        Vector<Point> simpsonPoints = new Vector<Point>();
//        String[] names = new String[]{"firstQ","secondQ","thirdQ","simpson"};
//
//        for (double i=5; i<=200; i++) {
//            integration.smartIntegration(i);
//            firstPoints.add(new Point(i, integration.difference[1]));
//            secondPoints.add(new Point(i,integration.difference[2]));
//            thirdPoints.add(new Point(i,integration.difference[3]));
//            simpsonPoints.add(new Point(i,integration.difference[4]));
//        }
//
//        plotData.add(firstPoints);
//        plotData.add(secondPoints);
//        plotData.add(thirdPoints);
//        plotData.add(simpsonPoints);

        for (int i=1; i<=10; i++) {
            Double diff = Integration.realValueWithP - integration.NewtonKotes(i);
            data.add(new Point(((double) i), diff));
        }

        integration.Gauss(3,1);




    }
};