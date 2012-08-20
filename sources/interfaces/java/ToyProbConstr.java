import tango.algencan.base.*;
import tango.algencan.util.*;

public class ToyProbConstr extends Constraints {

    ToyProbConstr() {

        super();

    }


    public double evalc(int n,double[] x,int ind)
        throws ALGENCANException {

        double c;

        switch( ind ) {

        case 1:
            c = Math.pow(x[0],2) + 1.0 - x[n - 1];
            break;

        case 2:
            c = 2.0 - x[0] - x[n - 1];
            break;

        default:
            throw new ALGENCANException(1, "Constraint " + ind +
                                        " unknown." );

        }

        return c;

    }

    public Jacobian evaljac(int n,double[] x,int ind)
        throws ALGENCANException {

        int      jcnnz;
        int[]    jcvar;
        double[] jcval;

        switch( ind ) {

        case 1:

            jcnnz = 2;

            jcvar = new    int[jcnnz];
            jcval = new double[jcnnz];

            jcvar[0] = 0;
            jcval[0] = 2.0 * x[0];

            jcvar[1] = 1;
            jcval[1] = - 1.0;

            break;

        case 2:

            jcnnz = 2;

            jcvar = new    int[jcnnz];
            jcval = new double[jcnnz];

            jcvar[0] = 0;
            jcval[0] = - 1.0;

            jcvar[1] = 1;
            jcval[1] = - 1.0;

            break;

        default:
            throw new ALGENCANException(1, "Constraint " + ind +
                                        " unknown." );
            
        }

        return new ALGENCANJacobian(jcnnz,jcvar,jcval);

    }

    public Hessian evalhc(int n,double[] x,int ind)
        throws ALGENCANException {

        int         hcnnz;
        int[] hclin,hccol;
        double[]    hcval;

        switch(ind) {

        case 1:
            
            hcnnz = 1;
            
            hclin = new    int[hcnnz];
            hccol = new    int[hcnnz];
            hcval = new double[hcnnz];

            hclin[0] = 0;
            hccol[0] = 0;
            hcval[0] = 2.0;

            break;

        case 2:

            hcnnz = 0;
            
            hclin = new    int[hcnnz];
            hccol = new    int[hcnnz];
            hcval = new double[hcnnz];

            break;

        default:
            throw new ALGENCANException(1,"Constraint " + ind +
                                        " unknown.");

        }

        return new ALGENCANHessian(hcnnz,hclin,hccol,hcval);

    }

}