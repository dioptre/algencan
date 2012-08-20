import tango.algencan.base.*;
import tango.algencan.util.*;

public class ToyProbObjFunc extends ObjectiveFunction {

    public ToyProbObjFunc() {

        super();

    }

    public double evalf(int n,double[] x) throws ALGENCANException {

        return x[1];

    }

    public double[] evalg(int n,double[] x) throws ALGENCANException {

        double[] g;

        g    = new double[2];

        g[0] = 0.0;
        g[1] = 1.0;

        return g;

    }

    public Hessian evalh(int n,double[] x) throws ALGENCANException {

        int      hnnz;
        int[]    hlin;
        int[]    hcol;
        double[] hval;
        
        hnnz = 0;

        hlin = new    int[hnnz];
        hcol = new    int[hnnz];
        hval = new double[hnnz];

        return new ALGENCANHessian(hnnz,hlin,hcol,hval);

    }
}
