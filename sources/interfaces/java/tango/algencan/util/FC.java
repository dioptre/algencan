package tango.algencan.util;

import tango.algencan.base.*;

public class FC implements ObjectiveAndConstraints {

    private double f;
    private double[] c;

    public FC(double f,double[] c) {

        this.f = f;
        this.c = c;

    }

    public double getF()   { return f; }
    public double[] getC() { return c; }

}
