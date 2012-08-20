package tango.algencan.util;

import tango.algencan.base.*;

public class ALGENCANHessian extends SparseStructure implements Hessian {

    public ALGENCANHessian(int hnnz,int[] hlin,int[] hcol,double[] hval) {

        super(hnnz,hlin,hcol,hval);

    }

    public int      getHnnz() { return  nnz; }
    public int[]    getHlin() { return  lin; }
    public int[]    getHcol() { return  col; }
    public double[] getHval() { return  val; }

}