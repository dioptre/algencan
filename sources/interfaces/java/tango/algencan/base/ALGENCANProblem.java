package tango.algencan.base;

public class ALGENCANProblem {

    private Characteristics probChar;
    private ObjectiveFunction objF;
    private Constraints constr;
    private Combined combined;

    public double epsfeas,epsopt;
    public int iprint,ncomp;

    public ALGENCANProblem(Characteristics probChar,
                           ObjectiveFunction objF,
                           Constraints constr,
                           Combined combined) {

        this.probChar = probChar;
        this.objF = objF;
        this.constr = constr;
        this.combined = combined;

    }

    public ALGENCANProblem(Characteristics probChar,
                           ObjectiveFunction objF,
                           Constraints constr) {

        this(probChar,objF,constr,new Combined());

    }

    public ALGENCANProblem(Characteristics probChar,Combined combined) {

        this(probChar,new ObjectiveFunction(),new Constraints(),combined);

    }

    public Characteristics   getCharacteristics()   { return probChar; }
    public ObjectiveFunction getObjectiveFunction() { return     objF; }
    public Constraints       getConstraints()       { return   constr; }
    public Combined          getCombined()          { return combined; }

    public int getN()            { return probChar.getN(); }

    public int getM()            { return probChar.getM(); }

    public double[] getLBounds() { return probChar.getLBounds(); }

    public double[] getUBounds() { return probChar.getUBounds(); }

    public double[] getX()       { return probChar.getX(); }

    public boolean getCheckder() { return probChar.getCheckder(); }

    public boolean[] getCoded()  { return probChar.getCoded(); }

    public double[] getLambda()  { return probChar.getLambda(); }
    
    public boolean[] getEquatn() { return probChar.getEquatn(); }

    public boolean[] getLinear() { return probChar.getLinear(); }

    public void setAlgencanParam(double epsfeas,double epsopt,
                                 int iprint,int ncomp) {

        this.epsfeas = epsfeas;
        this.epsopt  =  epsopt;
        this.iprint  =  iprint;
        this.ncomp   =   ncomp;

    }

    /*
     * Solver
     */

    public Status algencan() {

        return Algencan(epsfeas,epsopt,iprint,ncomp,getN(),getX(),
                        getLBounds(),getUBounds(),getM(),getLambda(),
                        getEquatn(),getLinear(),getCoded(),getCheckder());
    }


    /*
     * Native method. Uses C - Fortran interface.
     */

    public native Status Algencan(double epsfeas,double epsopt,int iprint,
                         int ncomp,int n,double[] x,double[] l,double[] u,
                         int m,double[] lambda,boolean[] equatn,
                         boolean[] linear,boolean[] coded,boolean checkder);

    /*
     * Methods
     */
    

    public double evalf(double[] x)
        throws ALGENCANException {

        return objF.evalf(getN(),x);

    }


    public double[] evalg(double[] x)
        throws ALGENCANException {

        return objF.evalg(getN(),x);

    }

    public Hessian evalh(double[] x)
        throws ALGENCANException {

        return objF.evalh(getN(),x);

    }


    public double evalc(double[] x,int ind)
        throws ALGENCANException {

        return constr.evalc(getN(),x,ind);

    }

    public Jacobian evaljac(double[] x,int ind)
        throws ALGENCANException {

        return constr.evaljac(getN(),x,ind);

    }


    public Hessian evalhc(double[] x,int ind)
        throws ALGENCANException {

        return constr.evalhc(getN(),x,ind);

    }

    public ObjectiveAndConstraints evalfc(double[] x)
        throws ALGENCANException {

        return combined.evalfc(getN(),x,getM());

    }

    public GradientAndJacobian evalgjac(double[] x)
        throws ALGENCANException {

        return combined.evalgjac(getN(),x,getM());

    }

    public Hessian evalhl(double[] x,double[] lambda,double sf,
                          double[] sc)
        throws ALGENCANException {

        return combined.evalhl(getN(),x,getM(),lambda,sf,sc);
        
    }

    public HLTimesVector evalhlp(double[] x,double[] lambda,double sf,
                                 double[] sc,double[] p,boolean goth)
        throws ALGENCANException {
    
        return combined.evalhlp(getN(),x,getM(),lambda,sf,sc,p,goth);

    }

    public void endp() throws ALGENCANException {

        probChar.endp(getN(),getX(),getM(),getLambda());

    }

}
