c-----------------------------------------------------------------------------
      integer idm,jdm,kdm,nrelax
      parameter (idm=53,jdm=322,kdm=43,nrelax=5)
c
      integer ii,jj,kk,ii1,jj1,nlato,nlongo,nlatn,nlongn
      parameter (ii=idm,jj=jdm,kk=kdm,ii1=ii-1,jj1=jj-1)
      parameter (nlato=1,nlongo=1,nlatn=idm,nlongn=jdm)
c
      integer i,j,k,l,m,n,mm,nn,km,kn,k1m,k1n,lp
c
      common /imau/ lpout, bin_restart, idl_out
      logical lpout, bin_restart, idl_out

      common /linepr/ lp
c-----------------------------------------------------------------------------
