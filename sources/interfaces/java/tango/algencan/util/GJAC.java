package tango.algencan.util;

import tango.algencan.base.*;

public class GJAC extends SparseStructure
    implements GradientAndJacobian {

    private double[] grad;

    public GJAC(double[] grad,int jcnnz,int[] jcfun,int[] jcvar,
                double[] jcval) {

        super(jcnnz,jcfun,jcvar,jcval);

        this.grad = grad;

    }

    public double[] getG()     { return grad; }

    public int      getJCnnz() { return  nnz; }
    public int[]    getJCfun() { return  lin; }
    public int[]    getJCvar() { return  col; }
    public double[] getJCval() { return  val; }

}