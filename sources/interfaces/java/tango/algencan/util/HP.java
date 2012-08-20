package tango.algencan.util;

import tango.algencan.base.*;

public class HP extends SparseStructure
    implements HLTimesVector {

    private boolean goth;

    public HP(double[] hp,boolean goth) {

        super(hp.length,hp);
        this.goth = goth;

    }

    public double[] getHP()   {  return val; }
    public boolean  getGoth() { return goth; }

}