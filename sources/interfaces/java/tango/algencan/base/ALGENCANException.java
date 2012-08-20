package tango.algencan.base;

public class ALGENCANException extends Exception {

    private int flag;
    private String error;

    public ALGENCANException(int flag) {
        
        super();
        this.flag  = flag;
        this.error = "";

    }

    public ALGENCANException(int flag,String err) {

        this(flag);
        this.error = err;

    }

    public int getFlag() { return this.flag; }
    public String toString() {

        return " *** Java Interface - Flag " + this.flag + ": " +
            error + " ***\n";

    }

}