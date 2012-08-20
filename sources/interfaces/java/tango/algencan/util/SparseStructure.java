package tango.algencan.util;

public class SparseStructure {

    protected int nnz;
    protected int[] lin,col;
    protected double[] val;

    public SparseStructure() {

        nnz = 0;

    }

    public SparseStructure(int nnz, double[] v) {

        this.nnz = nnz;
        this.val = v;

    }

    public SparseStructure(int nnz, int[] lin, double[] v) {

        this(nnz,v);
        this.lin = lin;

    }

    public SparseStructure(int nnz, int[] lin, int[] col, double[] v) {

        this(nnz,lin,v);
        this.col = col;

    }

}