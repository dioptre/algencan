package tango.algencan.base;

public final class Status {

    public final int inform,m,n;
    public final double cnorm,f,nlpsupn,snorm;
    public final double[] lambda,x;

    public Status(int n,double[] x,int m,double[] lambda,double f,
                  double cnorm,double snorm,double nlpsupn,int inform) {

        this.n = n;
        this.m = m;
        this.x = x;
        this.f = f;
        this.lambda  =  lambda;
        this.cnorm   =   cnorm;
        this.nlpsupn = nlpsupn;
        this.snorm   =   snorm;
        this.inform  =  inform;

    }

    public String toString() {

        String s;
        int i;

        s = String.format(
            "\nPrinting Status information\n\n" +
            "Flag of ALGENCAN                                %19c%5d\n" +
            "Functional Value                              = %24.16E\n" +
            "Sup-norm of constraints                       = %17c%7.1E\n" +
            "Sup-norm of complementarity-feasibility       = %17c%7.1E\n" +
            "Sup-norm of the Lagrangian projected gradient = %17c%7.1E\n" +
            "\nFinal point (first %7d components):\n",
            ' ',inform,f,' ',cnorm,' ',snorm,' ',nlpsupn,Math.min(n,6));

        for ( i = 0; i < Math.min(n,6); i++ )
            s = s.concat( String.format( "%c%11.4E", ' ',x[i] ) );

        s = s.concat( String.format(
            "\n\nFinal Lagrange multipliers (first %7d components):\n",
            Math.min(m,6) ) );
                      
        for ( i = 0; i < Math.min(m,6); i++ )
            s = s.concat( String.format( "%c%11.4E", ' ',lambda[i] ) );

        s = s.concat( String.format( "\n" ) );

        return s;

    }

    public int getInform()      { return inform;  }
    public double getCnorm()    { return cnorm;   }
    public double getF()        { return f;       }
    public double getNlpsupn()  { return nlpsupn; }
    public double getSnorm()    { return snorm;   }
    public double[] getLambda() { return lambda;  }
    public double[] getX()      { return x;       }

}
