/* =================================================================
   File: algencanma.c
   =================================================================

   =================================================================
   Module: Main program
   =================================================================

   Last update of any of the component of this module:

   March 5, 2008.

   Users are encouraged to download periodically updated versions of
   this code at the TANGO home page:

   www.ime.usp.br/~egbirgin/tango/

   *****************************************************************
   *****************************************************************

   TANGO Project

   License:

   All the TANGO Project components are free software; you can
   redistribute it and/or modify it under the terms of the GNU General
   Public License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. You can also find the GPL on the GNU web site.

   Non-free versions of TANGO are available under terms different from
   those of the General Public License. Professors J. M. Martínez
   (martinez@ime.unicamp.br, martinezimecc@gmail.com) or E. G. Birgin
   (egbirgin@ime.usp.br, egbirgin@gmail.com) should be contacted for
   more information related to such a license, future developments
   and/or technical support.

   In addition, we kindly ask you to acknowledge the TANGO Project and
   its authors in any program or publication in which you use of the
   TANGO Project components. For published works that use ALGENCAN we
   suggest referencing:

   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, "On
   Augmented Lagrangian methods with general lower-level constraints",
   to appear in SIAM Journal on Optimization.

   and

   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt,
   "Augmented Lagrangian methods under the Constant Positive Linear
   Dependence constraint qualification", Mathematical Programming 111,
   pp. 5-32, 2008

   For published works that use GENCAN we suggest referencing:

   E. G. Birgin and J. M. Martínez, "Large-scale active-set
   box-constrained optimization method with spectral projected
   gradients", Computational Optimization and Applications 23,
   pp. 101-125, 2002.

   For published works that use SPG we suggest referencing:

   E. G. Birgin, J. M. Martínez and M. Raydan, "Nonmonotone spectral
   projected gradient methods on convex sets", SIAM Journal on
   Optimization 10, pp. 1196-1211, 2000,

   and

   E. G. Birgin, J. M. Martínez and M. Raydan, "Algorithm 813: SPG -
   software for convex-constrained optimization", ACM Transactions on
   Mathematical Software 27, pp. 340-349, 2001.

   (See also other related works in the TANGO Project home page:

   http://www.ime.usp.br/~egbirgin/tango/)

   ***************************************************************** */

/* ******************************************************************
   ****************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "algencanma.h"

int main() {

   int checkder,inform,iprint,m,n,ncomp;
   double cnorm,f,nlpsupn,epsfeas,epsopt,snorm;

   int coded[10];
   int *equatn,*linear;
   double *l,*lambda,*u,*x;

/* SET SOME SOLVER ARGUMENTS */
   param(&epsfeas,&epsopt,&iprint,&ncomp);

/* SET UP PROBLEM DATA */
   inip(&n,&x,&l,&u,&m,&lambda,&equatn,&linear,coded,&checkder);

/* CALL OPTIMIZATION SOLVER */
   algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,
   linear,coded,checkder,&f,&cnorm,&snorm,&nlpsupn,&inform);

/* WRITE ADDITIONAL OUTPUT INFORMATION CODED BY THE USER */
   endp(n,x,l,u,m,lambda,equatn,linear);

   return 0;
}

/* ******************************************************************
   ****************************************************************** */

void param(double *epsfeas,double *epsopt,int *iprint,int *ncomp)
{
  *epsfeas = 1.0e-08;
  *epsopt  = 1.0e-08;

  *iprint = 10;
  *ncomp  = 6;
}
