import tango.algencan.base.*;

public class Algencanma {

    /*
     * Load DLL (or shared library) which contains
     * implementation of native methods.
     *
     * On Linux, it looks for libAlgencan.so
     */

    static  {
        
        System.load(System.getProperty("user.dir") +
                    System.getProperty("file.separator") +
                    System.mapLibraryName("Algencan"));

    }

    public static void main( String[] args ) {

        ALGENCANProblem     p;
        Algencanma algencanma;
        Status         status;

        Characteristics    ch;
        ObjectiveFunction  ob;
        Constraints        co;

        double epsfeas,epsopt;
        int      iprint,ncomp;

        try {

            /* Change here to solve your own problem */

            epsfeas = 1.0e-08;
            epsopt  = 1.0e-08;
            
            iprint  = 10;
            ncomp   = 6;
            
            ch = new ToyProbCharacteristics();
            ob = new ToyProbObjFunc();
            co = new ToyProbConstr();

            p = new ALGENCANProblem(ch,ob,co);
            
            /* Changes end here */
            
            p.setAlgencanParam(epsfeas,epsopt,iprint,ncomp);
            
            status = p.algencan();
            
            p.endp();

            System.out.println(status);

            System.exit(0);

        }
        catch( NotImplementedException e ) {
            
            System.out.println( e );
            System.exit(1);
            
        }
        catch( Exception e ) {
            
            System.out.println( e );
            System.exit(1);
        }

    }

}
