package tango.algencan.base;

public class NotImplementedException extends ALGENCANException {

    public NotImplementedException() {

        this("");

    }

    public NotImplementedException(String fname) {

        super(1,"Function " + fname + " not implemented.");

    }

}