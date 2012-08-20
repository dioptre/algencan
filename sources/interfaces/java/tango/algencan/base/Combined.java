package tango.algencan.base;

public class Combined {

    public ObjectiveAndConstraints evalfc(int n,double[] x,int m)
        throws ALGENCANException {

        throw new NotImplementedException( "evalfc" );

    }

    public GradientAndJacobian evalgjac(int n,double[] x,int m)
        throws ALGENCANException {

        throw new NotImplementedException( "evalgjac" );

    }

    public Hessian evalhl(int n,double[] x,int m,double[] lambda,
                          double sf,double[] sc)
        throws ALGENCANException {

        throw new NotImplementedException( "evalhl" );

    }

    public HLTimesVector evalhlp(int n,double[] x,int m,double[] lambda,
                                 double sf,double[] sc,double[] p,
                                 boolean goth)
        throws ALGENCANException {

        throw new NotImplementedException( "evalhlp" );

    }

}