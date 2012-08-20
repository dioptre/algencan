C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac
      integer nsteps,maxfrt,latop
      double precision ops

C     COMMON ARRAYS
      integer icntl(30),info(20),iw(nnzmax),iw1(2*nsysmax),
     +        ikeep(3*nsysmax),sta(nsysmax),next(nnzmax),curr(nsysmax),
     +        pos(nsysmax),invp(nsysmax)
      double precision s(nsysmax),sdiag(nsysmax),cntl(5),w(nsysmax),
     +        fact(nnzmax)

C     COMMON BLOCKS
      common /lssdat/ s,sdiag,cntl,w,fact,icntl,info,iw,iw1,ikeep,sta,
     +                next,curr,pos,invp,ops,nsteps,maxfrt,latop,
     +                lacneig,lsclsys,lusefac
      save   /lssdat/
