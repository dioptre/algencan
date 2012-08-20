package tango.algencan.base;

public interface GradientAndJacobian {

    public double[] getG();

    public int      getJCnnz();
    public int[]    getJCfun();
    public int[]    getJCvar();
    public double[] getJCval();

}