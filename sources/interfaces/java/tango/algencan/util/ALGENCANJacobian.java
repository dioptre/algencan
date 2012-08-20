package tango.algencan.util;

import tango.algencan.base.*;

public class ALGENCANJacobian extends SparseStructure implements Jacobian {

    public ALGENCANJacobian(int jcnnz,int[] jcvar,double[] jcval) {

        super(jcnnz,jcvar,jcval);

    }

    public int getJCnnz()      { return  nnz; }
    public int[] getJCvar()    { return  lin; }
    public double[] getJCval() { return  val; }

}