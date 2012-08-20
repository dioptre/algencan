package tango.algencan.base;

public interface Jacobian {

    public int      getJCnnz();
    public int[]    getJCvar();
    public double[] getJCval();

}