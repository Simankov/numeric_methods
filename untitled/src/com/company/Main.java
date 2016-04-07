package com.company;


public class Main {
    public static void main(String[] args) {
        Integration integration = new Integration();
        System.out.println(integration.AitkenProcess(0.01));
        System.out.println(integration.RungeMethod(0.00001,0.01));
        System.out.println(integration.RichardsonMethod(0.00001,0.01));

    }
};