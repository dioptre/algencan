param xini {i in 0..1};

var x {i in 0..1} >= -10, <= 10, := xini[i];

minimize Obj:
    x[1];

subject to R1:
    x[0]^2 + 1 - x[1] <= 0;

subject to R2:
    2 - x[0] - x[1] <= 0;

data;

param xini :=
    0   0
    1   0
;
