package tango.algencan.base;

public interface Hessian {

    public int      getHnnz();
    public int[]    getHlin();
    public int[]    getHcol();
    public double[] getHval();

}