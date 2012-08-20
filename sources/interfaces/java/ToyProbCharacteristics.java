import tango.algencan.base.*;
import tango.algencan.util.*;

public class ToyProbCharacteristics extends Characteristics {

    ToyProbCharacteristics() {

        super();

        int i;

        n = 2;
        m = 2;

        x = new double[n];

        x[0] = 0.0e0;
        x[1] = 0.0e0;

        l = new double[n];

        l[0] = - 10.0;
        l[1] = - Double.MAX_VALUE;

        u = new double[n];

        u[0] = 10.0;
        u[1] = Double.MAX_VALUE;


        coded = new boolean[10];

        coded[0] =  true; /* evalf    */
        coded[1] =  true; /* evalg    */
        coded[2] =  true; /* evalh    */
        coded[3] =  true; /* evalc    */
        coded[4] =  true; /* evaljac  */
        coded[5] =  true; /* evalhc   */
        coded[6] = false; /* evalfc   */
        coded[7] = false; /* evalgjac */
        coded[8] = false; /* evalhl   */
        coded[9] = false; /* evalhlp  */

        lambda = new double[m];

        for ( i = 0; i < m; i++ )
            lambda[i] = 0.0;

        equatn = new boolean[m];

        for ( i = 0; i < m; i++ )
            equatn[i] = false;

        linear = new boolean[m];

        linear[0] = false;
        linear[1] =  true;

        checkder  =  true;

    }

    public void endp(int n,double[] x,int m,double[] lambda) {}

}
