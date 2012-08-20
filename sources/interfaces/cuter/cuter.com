C     COMMON SCALARS
      logical useslacks
      integer nranges,mcuter,ncuter

C     COMMON ARRAYS
      integer ccor(mmax),cmap(mmax),slaind(mmax)
      double precision ca(mmax),cb(mmax)

C     COMMON BLOCKS
      common /probdata/ ca,cb,ccor,cmap,slaind,mcuter,ncuter,nranges,
     +                  useslacks
      save   /probdata/

