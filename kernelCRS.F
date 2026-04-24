      subroutine kernelCRS(n,nval,aptr,aind,aval,x,Ax)
      implicit none
      integer, intent(in) :: n,nval
      integer, intent(in) :: aptr(n+1),aind(nval)
      real*8, intent(in) :: aval(3,3,nval)
      real*8, intent(in) :: x(NVEC,3,n)
      real*8, intent(out) :: Ax(NVEC,3,n)
      integer :: i1, i2 ,j1, j2, k1, k2, ivec
      real*8 :: ax_tmp(3)

#ifdef _OPENACC
!$acc parallel loop collapse(2) gang vector
!$acc& present(aptr,aind,aval,x,Ax) private(ax_tmp)
#else
!$OMP PARALLEL DO default(none),shared(n,aptr,aind,aval,x,Ax),
!$OMP& private(i1,ivec,ax_tmp,j2,i2,k1,k2)
#endif
      do i1 = 1, n
        do ivec =1,NVEC
          ax_tmp(:) = 0.0
!$acc loop seq
          do j2 = aptr(i1)+1, aptr(i1+1)
            i2 = aind(j2)
!$acc loop seq
            do k1 = 1, 3
!$acc loop seq
              do k2 = 1, 3
                ax_tmp(k1) = ax_tmp(k1) + aval(k2,k1,j2) * x(ivec,k2,i2)
              enddo
            enddo
          enddo
!$acc loop seq
          do k1 = 1, 3
            Ax(ivec,k1,i1) = ax_tmp(k1)
          enddo
        enddo
      enddo
#ifdef _OPENACC
!$acc end parallel loop
#else
!$OMP END PARALLEL DO
#endif
      end
