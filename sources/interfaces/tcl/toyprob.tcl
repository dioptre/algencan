# =================================================================
# File: toyprob.tcl
# =================================================================

# =================================================================
# Module: Subroutines that define the problem
# =================================================================

# Last update of any of the component of this module: 

# April 28th, 2010.

# Users are encouraged to download periodically updated versions of 
# this code at the TANGO home page:

# www.ime.usp.br/~egbirgin/tango/ 

# ******************************************************************
# ******************************************************************

proc inip {pn px pl pu pm plambda pequatn plinear pcoded pcheckder} {

    upvar $pn               n
    upvar $px               x
    upvar $pl               l
    upvar $pu               u
    upvar $pm               m
    upvar $plambda     lambda
    upvar $pequatn     equatn
    upvar $plinear     linear
    upvar $pcoded       coded
    upvar $pcheckder checkder

    # Number of variables

    set n 2

    # Initial point

    for {set i 0} {$i < [expr {$n - 1}]} {incr i} {
        set x($i) 0.0
    }

    set x([expr {$n - 1}]) 0.0

    # Lower and upper bounds

    for {set i 0} {$i < [expr {$n - 1}]} {incr i} {
        set l($i) -10.0
        set u($i)  10.0
    }

    set l([expr {$n - 1}]) -1.0e20
    set u([expr {$n - 1}])  1.0e20

    # Number of constraints (equalities plus inequalities)

    set m 2

    for {set i 0} {$i < $m} {incr i} {
        set lambda($i) 0.0
    }

    for {set i 0} {$i < $m} {incr i} {
        set equatn($i) 0
    }

    set linear(0) 0
    set linear(1) 1

    # Coded:
    #
    # 0 - evalf
    # 1 - evalg
    # 2 - evalh
    # 3 - evalc
    # 4 - evaljac
    # 5 - evalhc
    # 6 - evalfc
    # 7 - evalgjac
    # 8 - evalhl
    # 9 - evalhlp

    set coded(0) 1
    set coded(1) 1
    set coded(2) 1
    set coded(3) 1
    set coded(4) 1
    set coded(5) 1
    set coded(6) 0
    set coded(7) 0
    set coded(8) 0
    set coded(9) 0

    # Checkder

    set checkder 1

}

# **********************************************************************
# **********************************************************************

proc evalf {n px pf pflag} {

    upvar $px       x
    upvar $pf       f
    upvar $pflag flag

    set flag 0

    set f $x(1)

}

# **********************************************************************
# **********************************************************************

proc evalg {n px pg pflag} {

    upvar $px       x
    upvar $pg       g
    upvar $pflag flag

    set flag 0

    set g(0) 0.0
    set g(1) 1.0

}

# **********************************************************************
# **********************************************************************

proc evalh {n px phlin phcol phval phnnz pflag} {

    upvar $px       x
    upvar $phlin hlin
    upvar $phcol hcol
    upvar $phval hval
    upvar $phnnz hnnz
    upvar $pflag flag

    set flag 0
    set hnnz 0

}

# **********************************************************************
# **********************************************************************

proc evalc {n px ind pcind pflag} {

    upvar $px       x
    upvar $pcind cind
    upvar $pflag flag

    set flag 0

    if {$ind == 1} {
        set cind [expr {pow($x(0), 2) + 1.0 - $x(1)}]
    } elseif {$ind == 2} {
        set cind [expr {2.0 - $x(0) - $x(1)}]
    } else {
        set flag -1;
    }

}

# **********************************************************************
# **********************************************************************

proc evaljac {n px ind pjcvar pjcval pjcnnz pflag} {

    upvar $px         x
    upvar $pjcvar jcvar
    upvar $pjcval jcval
    upvar $pjcnnz jcnnz
    upvar $pflag   flag
    
    set flag 0
    
    if {$ind == 1} {
        set jcnnz 2
        
        set jcvar(0) 1
        set jcval(0) [expr {2.0 * $x(0)}]
        
        set jcvar(1) 2
        set jcval(1) -1.0
    } elseif {$ind == 2} {
        set jcnnz 2
        
        set jcvar(0) 1
        set jcval(0) -1.0
        
        set jcvar(1) 2
        set jcval(1) -1.0
    } else {
        set flag -1;
    }
    
}

# **********************************************************************
# **********************************************************************

proc evalhc {n px ind phclin phccol phcval phcnnz pflag} {

    upvar $px         x
    upvar $phcnnz hcnnz
    upvar $phclin hclin
    upvar $phccol hccol
    upvar $phcval hcval
    upvar $pflag   flag

    set flag 0

    if {$ind == 1} {
        set hcnnz 1

        set hclin(0) 1
        set hccol(0) 1
        set hcval(0) 2.0
    } elseif {$ind == 2} {
        set hcnnz 0
    } else {
        set flag -1
    }

}

# **********************************************************************
# **********************************************************************

proc evalfc {n px pf m pc pflag} {

    upvar $px       x
    upvar $pf       f
    upvar $pc       c
    upvar $pflag flag

    set flag -1

}

# **********************************************************************
# **********************************************************************

proc evalgjac {n px pg m pjcfun pjcvar pjcval pjcnnz pflag} {

    upvar $px         x   
    upvar $pg         g
    upvar $pjcfun jcfun
    upvar $pjcvar jcvar
    upvar $pjcval jcval
    upvar $pjcnnz jcnnz
    upvar $pflag   flag

    set flag -1

}

# **********************************************************************
# **********************************************************************

proc evalhl {n px m plambda sf psc phllin phlcol phlval phlnnz pflag} {

    upvar $px           x
    upvar $plambda lambda
    upvar $psc         sc
    upvar $phllin   hllin
    upvar $phlcol   hlcol
    upvar $phlval   hlval
    upvar $phlnnz   hlnnz
    upvar $pflag     flag

    set flag -1

}

# **********************************************************************
# **********************************************************************

proc evalhlp {n px m plambda sf psc pp php gothl pflag} {

    upvar $px           x
    upvar $plambda lambda
    upvar $psc         sc
    upvar $pp           p
    upvar $php         hp
    upvar $pflag     flag

    set flag -1

}

# **********************************************************************
# **********************************************************************

proc endp {n px pl pu m plambda pequatn plinear} {
    
    upvar $px           x
    upvar $pl           l
    upvar $pu           u
    upvar $plambda lambda
    upvar $pequatn equatn
    upvar $plinear linear

}
