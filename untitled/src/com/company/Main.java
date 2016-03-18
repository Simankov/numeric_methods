package com.company;

import java.util.ArrayList;
import java.util.Vector;



public class Main {
    public static void main(String[] args) {
        Integration integration = new Integration();
        System.out.println(integration.AitkenProcess(0.00001));
        System.out.println(integration.RungeMethod(0.001,0.00001));
        System.out.println(integration.RichardsonMethod(0.001,0.00001));
    }
};