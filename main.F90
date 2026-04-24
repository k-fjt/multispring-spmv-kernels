!c_______________________________________________________________________
      program multispringspmvkernels
      implicit none
      integer n,ne,kd,ncb,n_size,ne_size,np,im
      
#ifdef CRS10
      write(6,*) "using CRS with NVEC",NVEC
#else
      write(6,*) "using EBE with NVEC",NVEC
#endif

!=== READ INPUT DATA SIZES ===
      im=0
      np=1
      
      open(10,file='./data/setting.dat',status='unknown')
      read(10,*) !'num of node' 
      read(10,*)  n
      read(10,*) !'num of tet elem' 
      read(10,*) ne
      read(10,*) !'num of material properties' 
      read(10,*) kd
      close(10)
      
      n_size=n
      ne_size=ne
      ncb=30*n ! maximum number of nonzero components in CRS
      call main(n,ne,kd,ncb,n_size,ne_size,np,im)
      end
!c_______________________________________________________________________
      subroutine main(n,ne,kd,ncb,n_size,ne_size,np,im)
      use omp_lib
      implicit none
      
      ! in scalars
      integer n,ne,kd,n_size,ne_size,np,im,ncb
      
      ! local scalars
      integer i,j,ivec,ie
      real*8 t1,t2
      real*8 dt
      real*8 fmin,fmax,alp,bet
      
      ! local arrays
      real*8 dupli(n_size)
      real*8 coor(3,n_size)
      real*8 up(NVEC,3*n_size),b(NVEC,3*n_size)
      integer num(ne_size),cny(10,ne_size)
      real*8 coe_merge(4,ne_size)
      real*8 rmat(kd,10)
      real*8 younglst(kd,2)
      real*8 d1(ne_size),d2(ne_size),d4(ne_size),d5(ne_size),d6(ne_size)
      real*8 tmp1(ne_size),tmp2(ne_size),tmpone(ne_size),tmpzero(ne_size)
      real*8 c1(ne_size),c2(ne_size),c3(ne_size),c4(ne_size)
      real*8 c5(ne_size),c6(ne_size),c7(ne_size),c8(ne_size)
      real*8 c9(ne_size),coe(ne_size)
#define ms_ng 5
      real*8 ms_Dmat(NVEC,6,6,ms_ng,ne)
#ifdef CRS10
      integer crsptr(n_size+1)
      integer crsind(ncb)
      real*8 crsval(9,ncb)
#endif

!=== READ INPUT DATA ===
      open(42,file='./data/num.bin',form='unformatted',status='old')
      read(42) num
      close(42)

      open(42,file='./data/conn.bin',form='unformatted',status='old')
      read(42) cny
      close(42)

      open(42,file='./data/coor.bin',form='unformatted',status='old')
      read(42) coor
      close(42)

      open(42,file='./data/material.dat',status='old')
      do j=1,kd
      read(42,*)
      do i=1,10
      read(42,*) rmat(j,i)
      enddo
      enddo
      close(42)

#ifdef CRS10
      open(42,file='./data/crsptr.bin',form='unformatted',status='old')
      read(42) crsptr
      close(42)
      open(42,file='./data/crsind.bin',form='unformatted',status='old')
      read(42) crsind
      close(42)
      open(42,file='./data/crsval.bin',form='unformatted',status='old')
      read(42) crsval
      close(42)
#endif

!=== SET MODEL DATA ===
      dt=0.005
      fmin=0.1
      fmax=5.0
      dupli=1.0
      coe_merge=1.0
      call makmatlist(kd,rmat,younglst)
      call calc_alp_bet(fmin,fmax,alp,bet)
      call setelementinfo(n,ne,n_size,ne_size,kd,num,rmat,younglst,dt, &
       coor,cny,alp,bet,coe_merge, &
       d1,d2,d4,d5,d6,tmp1,tmp2,c1,c2,c3,c4,c5,c6,c7,c8,c9,coe)
      
      ms_Dmat=0.0
      do ie=1,ne
      ms_Dmat(:,1,1,:,ie)=d1(ie)
      ms_Dmat(:,1,2,:,ie)=d2(ie)
      ms_Dmat(:,1,3,:,ie)=d2(ie)
      ms_Dmat(:,2,1,:,ie)=d2(ie)
      ms_Dmat(:,2,2,:,ie)=d1(ie)
      ms_Dmat(:,2,3,:,ie)=d2(ie)
      ms_Dmat(:,3,1,:,ie)=d2(ie)
      ms_Dmat(:,3,2,:,ie)=d2(ie)
      ms_Dmat(:,3,3,:,ie)=d1(ie)
      ms_Dmat(:,4,4,:,ie)=d4(ie)
      ms_Dmat(:,5,5,:,ie)=d5(ie)
      ms_Dmat(:,6,6,:,ie)=d6(ie)
      enddo

!=== SET MODEL INPUT ===
      do i=1,3*n
      do ivec=1,NVEC
      b(ivec,i)=cos(mod(i,100)*0.2)+sin(0.3*mod(i,100)*0.2)
      enddo
      enddo
      up=0.0d0

!=== ALLOCATE GPU MEMORY ===
#ifdef _OPENACC
      if(im.eq.0)then
      write(*,*) 'using OpenACC'
      endif
!$acc data copyin(b,up,cny,num,coor,rmat,younglst,dupli, &
#ifdef CRS10
!$acc& crsptr,crsind,crsval, &
#endif
!$acc& ms_Dmat,c1,c2,c3,c4,c5,c6,c7,c8,c9,tmp1,tmp2,coe)
#else
      if(im.eq.0)then
      write(*,*) 'using OpenMP'
      endif
#endif

!=== MEASURE MATRIX-VECTOR PRODUCT PERFORMANCE ===
      do i=1,10
      t1=omp_get_wtime()
!c-------------------------------------------------------------------
#ifdef CRS10
      call kernelCRS(n,ncb,crsptr,crsind,crsval,b,up)
#else
      call kernelEBE(n_size,ne_size,n,ne,b,up,&
       cny,c1,c2,c3,c4,c5,c6,c7,c8,c9,tmp1,tmp2,coe,ms_Dmat)
#endif
!c-------------------------------------------------------------------
      t2=omp_get_wtime()
      write(*,*) 'matvec took',t2-t1
      enddo

!=== PRINT MATRIX-VECTOR PRODUCT RESULTS ===
      call print_norm(0,im,n,n_size,b,dupli)
      call print_norm(1,im,n,n_size,up,dupli)

!=== DEALLOCATE GPU MEMORY ===
#ifdef _OPENACC
!$acc update host(b,up)
!$acc end data
#endif
      end
!c_______________________________________________________________________
      subroutine print_norm(inttype,im,n,n_size,b,dupli)
      implicit none
      integer n,i,im,ierr,inttype,n_size,ivec
      real*8 b(NVEC,3*n_size),dupli(n_size),tmp1(NVEC),tmp2(NVEC)
      real*8 zq_dum(1)

      call inner_product_d(n_size,n,b,b,tmp2,dupli)

      do ivec=1,NVEC
      tmp1(ivec)=sqrt(tmp2(ivec))
      enddo
      if(im.eq.0)then
      write(*,*) 'norm',inttype,tmp1
      endif
      end
!c____________________________________________________________________
      subroutine setelementinfo(n,ne,n_size,ne_size,kd,num,rmat,younglst,dt, &
       coor,cny,alp,bet,coe_merge, &
       d1,d2,d4,d5,d6,tmp1,tmp2,c1,c2,c3,c4,c5,c6,c7,c8,c9,coe)
      implicit none
      integer n,n_size,ne_size,ne,kd
      real*8 coor(3,n_size)
      integer*4 num(ne_size)
      real*8 rmat(kd,10),dt,alp,bet,coe_merge(4,ne_size)
      real*8 younglst(kd,2)
      integer*4 cny(10,ne_size)
      real*8 d1(ne_size),d2(ne_size),d4(ne_size),d5(ne_size),d6(ne_size)
      real*8 tmp1(ne_size),tmp2(ne_size)
      real*8 c1(ne_size),c2(ne_size),c3(ne_size),c4(ne_size)
      real*8 c5(ne_size),c6(ne_size),c7(ne_size),c8(ne_size)
      real*8 c9(ne_size),coe(ne_size)
      real*8 rho,young,rnyu,alpha,beta
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,detjmat
      integer ie,in

      do ie=1,ne
      in=num(ie)
      rho=rmat(in,3)
      young=younglst(in,1)
      rnyu=younglst(in,2)
      alpha=rmat(in,4)*alp *coe_merge(1,ie)
      beta=rmat(in,4)*bet *coe_merge(1,ie)
      tmp1(ie) = 4.0/dt/dt + 2.0*alpha/dt
      tmp2(ie) = 1.0 + 2.0*beta/dt
      x1 = coor(1,cny(1,ie)) - coor(1,cny(4,ie))
      y1 = coor(2,cny(1,ie)) - coor(2,cny(4,ie))
      z1 = coor(3,cny(1,ie)) - coor(3,cny(4,ie))
      x2 = coor(1,cny(2,ie)) - coor(1,cny(4,ie))
      y2 = coor(2,cny(2,ie)) - coor(2,cny(4,ie))
      z2 = coor(3,cny(2,ie)) - coor(3,cny(4,ie))
      x3 = coor(1,cny(3,ie)) - coor(1,cny(4,ie))
      y3 = coor(2,cny(3,ie)) - coor(2,cny(4,ie))
      z3 = coor(3,cny(3,ie)) - coor(3,cny(4,ie))
      c1(ie) = (y2*z1 - y1*z2)
      c2(ie) = (y1*z3 - y3*z1)
      c3(ie) = (y3*z2 - y2*z3)
      c4(ie) = (x1*z2 - x2*z1)
      c5(ie) = (x3*z1 - x1*z3)
      c6(ie) = (x2*z3 - x3*z2)
      c7(ie) = (x2*y1 - x1*y2)
      c8(ie) = (x1*y3 - x3*y1)
      c9(ie) = (x3*y2 - x2*y3)
      detjmat=x3*y2*z1+y3*z2*x1+z3*x2*y1-x3*z2*y1-z3*y2*x1-y3*x2*z1
      coe(ie)=1.0/detjmat
      tmp1(ie)=tmp1(ie)*rho*detjmat
      d1(ie)=((1.-rnyu)*young)/((1.-2.*rnyu)*(1.+rnyu))
      d2(ie)=(rnyu*young)/((1.-2.*rnyu)*(1. + rnyu))
      d4(ie)=young/(2.0*(1+rnyu)) *coe_merge(2,ie)
      d5(ie)=young/(2.0*(1+rnyu)) *coe_merge(3,ie)
      d6(ie)=young/(2.0*(1+rnyu)) *coe_merge(4,ie)
      enddo
      end
!c_______________________________________________________________________
      subroutine makmatlist(kd,rmat,younglst)
      implicit none
      integer kd,it
      real*8 rmat(kd,10)
      real*8 younglst(kd,2)
      real*8 c1tmp,c2tmp,rhtmp,ygtmp,rntmp
      
      do it=1,kd
      c1tmp=rmat(it,1)
      c2tmp=rmat(it,2)
      rhtmp=rmat(it,3)
      call mak_E_nyu(c1tmp,c2tmp,rhtmp,ygtmp,rntmp)
      younglst(it,1)=ygtmp
      younglst(it,2)=rntmp
      enddo
      end
!_______________________________________________________________________
      subroutine mak_E_nyu(c1,c2,rho,E,nyu)
      implicit none
      real*8 c1,c2,rho
      real*8 E,nyu
      E=((3*c1**2*c2**2-4*c2**4)*rho)/(c1**2-c2**2)
      nyu=(c1**2-2*c2**2)/(2*c1**2-2*c2**2)
      end
!_______________________________________________________________________
      subroutine calc_alp_bet(fmin,fmax,alp,bet)
      implicit none
      real*8 alp,bet,fmin,fmax,pi
      pi=atan(1.0)*4.0
      alp=(-2*fmax*fmin*Pi*(3*(fmax**2 - fmin**2) +  &
           2*(fmax**2 + fmax*fmin + fmin**2)*Log(fmin/fmax)))/ &
       (fmax - fmin)**3
      bet=(3*(fmax**2 - fmin**2 + 2*fmax*fmin*Log(fmin/fmax)))/ &
       (2.*(fmax - fmin)**3*Pi) 
      end
!c____________________________________________________________________
      subroutine inner_product_d(n_size,n,z,q,zq,dupli)
      implicit none
      ! in
      integer n,n_size
      real*8 z(NVEC,3,n_size),q(NVEC,3,n_size)
      real*8 dupli(n_size)
      ! out
      real*8 zq(NVEC)
      ! local
      integer i,ivec
      real*8 zq1,zq2,zq3,zq4

      zq=0.0d0
#ifdef _OPENACC

      zq1=0.0d0
#if NVEC>=2
      zq2=0.0d0
#endif
#if NVEC>=4
      zq3=0.0d0
      zq4=0.0d0
#endif
#if NVEC>4
      error NVEC must be 1,2,4,8 or 16
#endif
!$acc parallel loop present(dupli,z,q),reduction(+:zq1 &
#if NVEC>=2
!$acc ,zq2 &
#endif
#if NVEC>=4
!$acc ,zq3,zq4 &
#endif
!$acc )
      do i=1,n
      zq1=zq1+dupli(i)*(z(1,1,i)*q(1,1,i)+z(1,2,i)*q(1,2,i)+z(1,3,i)*q(1,3,i))
#if NVEC>=2
      zq2=zq2+dupli(i)*(z(2,1,i)*q(2,1,i)+z(2,2,i)*q(2,2,i)+z(2,3,i)*q(2,3,i))
#endif
#if NVEC>=4
      zq3=zq3+dupli(i)*(z(3,1,i)*q(3,1,i)+z(3,2,i)*q(3,2,i)+z(3,3,i)*q(3,3,i))
      zq4=zq4+dupli(i)*(z(4,1,i)*q(4,1,i)+z(4,2,i)*q(4,2,i)+z(4,3,i)*q(4,3,i))
#endif
      enddo
!$acc end parallel loop
      zq(1)=zq1
#if NVEC>=2
      zq(2)=zq2
#endif
#if NVEC>=4
      zq(3)=zq3
      zq(4)=zq4
#endif

#else

!$OMP PARALLEL DO default(none) shared (n,z,q,dupli),reduction(+:zq) private (i,ivec)
      do i=1,n
      do ivec=1,NVEC
      zq(ivec)=zq(ivec)+ &
           (z(ivec,1,i)*q(ivec,1,i) &
           +z(ivec,2,i)*q(ivec,2,i) &
           +z(ivec,3,i)*q(ivec,3,i))*dupli(i)
      enddo
      enddo
!$OMP END PARALLEL DO
#endif
#ifdef PRINTLOG
      write(6,*) 'inner_product',zq
#endif
      end
!c____________________________________________________________________

