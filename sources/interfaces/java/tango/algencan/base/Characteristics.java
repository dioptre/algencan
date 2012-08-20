package tango.algencan.base;

/**
 * This class is responsible for providing the main characteristics
 * of a NLP problem, such as number of variables, initial point, lower
 * and upper bounds, etc. Every NLP problem HAS to have its own extension
 * of this class.
 **/

public class Characteristics {
    
    
    protected int n,m;
    protected boolean checkder;
    protected boolean[] coded,equatn,linear;
    protected double[] x,l,lambda,u;

    /**
     * Returns the dimension the problem.
     *
     * @return the number 'n' of variables of the NLP problem.
     * @throws NotImplementedException in case this method was not
     *         implemented for an especific problem.
     **/

    public int getN() { return n; }

    public int getM() { return m; }

    /**
     * Returns the choice of testing the derivatives or not.
     *
     * @return true if the derivatives will be tested or false if
     *         if they will not.
     * @throws NotImplementedException in case this method was not
     *         implemented for an especific problem.     
     **/

    public boolean getCheckder() { return checkder; }

    /**
     * Returns the initial point. ALGENCAN also uses this alocated
     * vector as working vector inside its iterations.
     *
     * @return the starting point for the solver.
     * @throws NotImplementedException in case this method was not
     *         implemented for an especific problem.
     **/

    public double[] getX() { return x; }

    /**
     * Returns the lower bounds of the variables in the NLP problem.
     * Note that even if the variable is unbounded, its lower
     * bound must exists (e.g. a positive number with very large
     * module).
     *
     * @return the lower bounds for each variable.
     * @throws NotImplementedException in case this method was not
     *         implemented for an especific problem.
     **/
    
    public double[] getLBounds() { return l; }
    
    /**
     * Returns the upper bounds of the variables in the NLP problem.
     * Note that even if the variable is unbounded, its upper
     * bound must exists (e.g. a negative number with very large
     * module).
     *
     * @return the upper bounds for each varible.
     * @throws NotImplementedException in case this method was not
     *         implemented for an especific problem.
     **/
    
    public double[] getUBounds() { return u; }
    
    public double[] getLambda() { return lambda; }
    
    /**
     * Returns a boolean array with length 10 where position 'i'
     * indicates to ALGENCAN that subroutine 'i' was implemented
     * (true) or not (false) for the specific NLP problem. The order
     * is:
     * <P> 
     * 0 = evalf    <BR>
     * 1 = evalg    <BR>
     * 2 = evalh    <BR>
     * 3 = evalc    <BR>
     * 4 = evaljac  <BR>
     * 5 = evalhc   <BR>
     * 6 = evalfc   <BR>
     * 7 = evalgjac <BR>
     * 8 = evalhl   <BR>
     * 9 = evalhlp  
     *
     * @return a boolean array indicating which subrotine was
     *         implemented.
     * @throws NotImplementedException in case this method was not
     *         implemented for an especific problem.
     **/
    
    public boolean[] getCoded() { return coded; }
    
    public boolean[] getEquatn() { return equatn; }

    public boolean[] getLinear() { return linear; }

    /**
     * This method allows the user to work with the solution returned by
     * ALGENCAN and other data.
     **/

    public void endp(int n,double[] x,int m,double[] lambda)
        throws NotImplementedException {
        
        throw new NotImplementedException( "endp" );
        
    }
    
}
