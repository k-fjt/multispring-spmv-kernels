#define m00238 (0.002380952380952381)
#define mm0238 (-0.002380952380952381)
#define m00039 (0.0003968253968253968)
#define mm0158 (-0.0015873015873015873)
#define m01269 (0.012698412698412698)
#define m00634 (0.006349206349206349)
#define m00317 (0.0031746031746031746)
#define NIP 5

      subroutine kernelEBE
     & (n_size,ne_size,n,ne,ug,rg,
     &  cny,c1,c2,c3,c4,c5,c6,c7,c8,c9,tmp1,tmp2,coe,ms_Dmat)
      implicit none
      ! in
      integer n,ne,n_size,ne_size
      real*8 ug(NVEC,3,n_size)
      integer cny(10,ne_size),num(ne_size)
      real*8 c1(ne_size),c2(ne_size),c3(ne_size),c4(ne_size)
      real*8 c5(ne_size),c6(ne_size),c7(ne_size),c8(ne_size)
      real*8 c9(ne_size),coe(ne_size)
      real*8 tmp1(NVEC,ne_size),tmp2(NVEC,ne_size)
#define ms_ng 5
      real*8 ms_Dmat(NVEC,6,6,ms_ng,ne)
      ! out
      real*8 rg(NVEC,3,n_size)
      ! local
      real*8 alpha_iip(3,NIP),weight(NIP)
      real*8  u1, u2, u3, u4, u5
      real*8  u6, u7, u8, u9,u10
      real*8 u11,u12,u13,u14,u15
      real*8 u16,u17,u18,u19,u20
      real*8 u21,u22,u23,u24,u25
      real*8 u26,u27,u28,u29,u30
      real*8  BDBu1, BDBu2, BDBu3, BDBu4, BDBu5
      real*8  BDBu6, BDBu7, BDBu8, BDBu9, BDBu10
      real*8 BDBu11,BDBu12,BDBu13,BDBu14,BDBu15
      real*8 BDBu16,BDBu17,BDBu18,BDBu19,BDBu20
      real*8 BDBu21,BDBu22,BDBu23,BDBu24,BDBu25
      real*8 BDBu26,BDBu27,BDBu28,BDBu29,BDBu30
      real*8 Bu1,Bu2,Bu3,Bu4,Bu5,Bu6
      real*8 DBu1,DBu2,DBu3,DBu4,DBu5,DBu6
      integer iu,i,ie,ivec
      real*8 detj
      integer nabs,nabe,ieo
      integer iip
      real*8 dadx11,dadx12,dadx13,dadx21,dadx22
      real*8 dadx23,dadx31,dadx32,dadx33
      real*8 dNda11,dNda12,dNda13,dNda14,dNda15
      real*8 dNda16,dNda17,dNda18,dNda19,dNda1a
      real*8 dNda21,dNda22,dNda23,dNda24,dNda25
      real*8 dNda26,dNda27,dNda28,dNda29,dNda2a
      real*8 dNda31,dNda32,dNda33,dNda34,dNda35
      real*8 dNda36,dNda37,dNda38,dNda39,dNda3a
      real*8 dNdx11,dNdx12,dNdx13,dNdx14,dNdx15
      real*8 dNdx16,dNdx17,dNdx18,dNdx19,dNdx1a
      real*8 dNdx21,dNdx22,dNdx23,dNdx24,dNdx25
      real*8 dNdx26,dNdx27,dNdx28,dNdx29,dNdx2a
      real*8 dNdx31,dNdx32,dNdx33,dNdx34,dNdx35
      real*8 dNdx36,dNdx37,dNdx38,dNdx39,dNdx3a
      real*8 BDBu_iip1, BDBu_iip2, BDBu_iip3, BDBu_iip4, BDBu_iip5
      real*8 BDBu_iip6, BDBu_iip7, BDBu_iip8, BDBu_iip9,BDBu_iip10
      real*8 BDBu_iip11,BDBu_iip12,BDBu_iip13,BDBu_iip14,BDBu_iip15
      real*8 BDBu_iip16,BDBu_iip17,BDBu_iip18,BDBu_iip19,BDBu_iip20
      real*8 BDBu_iip21,BDBu_iip22,BDBu_iip23,BDBu_iip24,BDBu_iip25
      real*8 BDBu_iip26,BDBu_iip27,BDBu_iip28,BDBu_iip29,BDBu_iip30
      real*8 D_11,D_12,D_13,D_22,D_23,D_33,D_44,D_55,D_66
      real*8 tmp1i,tmp2i

      weight(5)=-4.0d0/30.0d0
      weight(1:4)=3.0d0/40.0d0
      alpha_iip(1:3,5) = 0.250
      alpha_iip(1:3,4) = (/1.0/6.0,1.0/6.0,1.0/6.0/)
      alpha_iip(1:3,3) = (/0.5,1.0/6.0,1.0/6.0/)
      alpha_iip(1:3,2) = (/1.0/6.0,0.5,1.0/6.0/)
      alpha_iip(1:3,1) = (/1.0/6.0,1.0/6.0,0.5/)
      
#ifdef _OPENACC
!$acc parallel loop present(rg) collapse(2)
#endif
      do i=1,n
      do ivec=1,NVEC
      rg(ivec,1,i)=0.0d0
      rg(ivec,2,i)=0.0d0
      rg(ivec,3,i)=0.0d0
      enddo
      enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif

#ifdef _OPENACC
!$acc parallel loop collapse(2) present(c1,c2,c3,c4,c5,c6,c7,c8,c9,
!$acc& ug,rg,tmp1,tmp2,coe,cny,ms_Dmat)
#endif
      do ie=1,ne
      do ivec=1,NVEC
      u1=ug(ivec,1,cny(1,ie))
      u2=ug(ivec,2,cny(1,ie))
      u3=ug(ivec,3,cny(1,ie))
      u4=ug(ivec,1,cny(2,ie))
      u5=ug(ivec,2,cny(2,ie))
      u6=ug(ivec,3,cny(2,ie))
      u7=ug(ivec,1,cny(3,ie))
      u8=ug(ivec,2,cny(3,ie))
      u9=ug(ivec,3,cny(3,ie))
      u10=ug(ivec,1,cny(4,ie))
      u11=ug(ivec,2,cny(4,ie))
      u12=ug(ivec,3,cny(4,ie))
      u13=ug(ivec,1,cny(5,ie))
      u14=ug(ivec,2,cny(5,ie))
      u15=ug(ivec,3,cny(5,ie))
      u16=ug(ivec,1,cny(6,ie))
      u17=ug(ivec,2,cny(6,ie))
      u18=ug(ivec,3,cny(6,ie))
      u19=ug(ivec,1,cny(7,ie))
      u20=ug(ivec,2,cny(7,ie))
      u21=ug(ivec,3,cny(7,ie))
      u22=ug(ivec,1,cny(8,ie))
      u23=ug(ivec,2,cny(8,ie))
      u24=ug(ivec,3,cny(8,ie))
      u25=ug(ivec,1,cny(9,ie))
      u26=ug(ivec,2,cny(9,ie))
      u27=ug(ivec,3,cny(9,ie))
      u28=ug(ivec,1,cny(10,ie))
      u29=ug(ivec,2,cny(10,ie))
      u30=ug(ivec,3,cny(10,ie))

      detj=coe(ie)
      dadx11 = c1(ie)
      dadx21 = c2(ie)
      dadx31 = c3(ie)
      dadx12 = c4(ie)
      dadx22 = c5(ie)
      dadx32 = c6(ie)
      dadx13 = c7(ie)
      dadx23 = c8(ie)
      dadx33 = c9(ie)
      
      BDBu1  = 0.0d0
      BDBu2  = 0.0d0
      BDBu3  = 0.0d0
      BDBu4  = 0.0d0
      BDBu5  = 0.0d0
      BDBu6  = 0.0d0
      BDBu7  = 0.0d0
      BDBu8  = 0.0d0
      BDBu9  = 0.0d0
      BDBu10 = 0.0d0
      BDBu11 = 0.0d0
      BDBu12 = 0.0d0
      BDBu13 = 0.0d0
      BDBu14 = 0.0d0
      BDBu15 = 0.0d0
      BDBu16 = 0.0d0
      BDBu17 = 0.0d0
      BDBu18 = 0.0d0
      BDBu19 = 0.0d0
      BDBu20 = 0.0d0
      BDBu21 = 0.0d0
      BDBu22 = 0.0d0
      BDBu23 = 0.0d0
      BDBu24 = 0.0d0
      BDBu25 = 0.0d0
      BDBu26 = 0.0d0
      BDBu27 = 0.0d0
      BDBu28 = 0.0d0
      BDBu29 = 0.0d0
      BDBu30 = 0.0d0

!$acc loop seq
      do iip = 1, NIP
        
        dNda31 = 4.0d0*alpha_iip(3,iip) - 1.0d0
        dNda22 = 4.0d0*alpha_iip(2,iip) - 1.0d0
        dNda13 = 4.0d0*alpha_iip(1,iip) - 1.0d0
        dNda14 = 4.0d0*(alpha_iip(1,iip)+alpha_iip(2,iip)
     &                 +alpha_iip(3,iip)) - 3.0d0
        dNda24 = 4.0d0*(alpha_iip(1,iip)+alpha_iip(2,iip)
     &                 +alpha_iip(3,iip)) - 3.0d0
        dNda34 = 4.0d0*(alpha_iip(1,iip)+alpha_iip(2,iip)
     &                 +alpha_iip(3,iip)) - 3.0d0
        dNda25 = 4.0d0*alpha_iip(3,iip)
        dNda35 = 4.0d0*alpha_iip(2,iip)
        dNda16 = 4.0d0*alpha_iip(2,iip)
        dNda26 = 4.0d0*alpha_iip(1,iip)
        dNda17 = 4.0d0*alpha_iip(3,iip)
        dNda37 = 4.0d0*alpha_iip(1,iip)
        dNda18 = -4.0d0*alpha_iip(3,iip)
        dNda28 = -4.0d0*alpha_iip(3,iip)
        dNda38 = -4.0d0*alpha_iip(1,iip) -4.0d0*alpha_iip(2,iip)
     &            -8.0d0*alpha_iip(3,iip) + 4.0d0
        dNda19 = -4.0d0*alpha_iip(2,iip)
        dNda29 = -4.0d0*alpha_iip(1,iip) -8.0d0*alpha_iip(2,iip)
     &            -4.0d0*alpha_iip(3,iip) + 4.0d0
        dNda39 = -4.0d0*alpha_iip(2,iip)
        dNda1a = -8.0d0*alpha_iip(1,iip) -4.0d0*alpha_iip(2,iip)
     &            -4.0d0*alpha_iip(3,iip) + 4.0d0
        dNda2a = -4.0d0*alpha_iip(1,iip)
        dNda3a = -4.0d0*alpha_iip(1,iip)

!dNdx
        dNdx11 = dNda31*dadx31
        dNdx12 = dNda22*dadx21
        dNdx13 = dNda13*dadx11
        dNdx14 = dNda14*dadx11+dNda24*dadx21+dNda34*dadx31
        dNdx15 = dNda25*dadx21+dNda35*dadx31
        dNdx16 = dNda16*dadx11+dNda26*dadx21
        dNdx17 = dNda17*dadx11+dNda37*dadx31
        dNdx18 = dNda18*dadx11+dNda28*dadx21+dNda38*dadx31
        dNdx19 = dNda19*dadx11+dNda29*dadx21+dNda39*dadx31
        dNdx1a = dNda1a*dadx11+dNda2a*dadx21+dNda3a*dadx31
        dNdx21 = dNda31*dadx32
        dNdx22 = dNda22*dadx22
        dNdx23 = dNda13*dadx12
        dNdx24 = dNda14*dadx12+dNda24*dadx22+dNda34*dadx32
        dNdx25 = dNda25*dadx22+dNda35*dadx32
        dNdx26 = dNda16*dadx12+dNda26*dadx22
        dNdx27 = dNda17*dadx12+dNda37*dadx32
        dNdx28 = dNda18*dadx12+dNda28*dadx22+dNda38*dadx32
        dNdx29 = dNda19*dadx12+dNda29*dadx22+dNda39*dadx32
        dNdx2a = dNda1a*dadx12+dNda2a*dadx22+dNda3a*dadx32
        dNdx31 = dNda31*dadx33
        dNdx32 = dNda22*dadx23
        dNdx33 = dNda13*dadx13
        dNdx34 = dNda14*dadx13+dNda24*dadx23+dNda34*dadx33
        dNdx35 = dNda25*dadx23+dNda35*dadx33
        dNdx36 = dNda16*dadx13+dNda26*dadx23
        dNdx37 = dNda17*dadx13+dNda37*dadx33
        dNdx38 = dNda18*dadx13+dNda28*dadx23+dNda38*dadx33
        dNdx39 = dNda19*dadx13+dNda29*dadx23+dNda39*dadx33
        dNdx3a = dNda1a*dadx13+dNda2a*dadx23+dNda3a*dadx33

        Bu1 = u1*dNdx11
        Bu2 = u2*dNdx21
        Bu3 = u3*dNdx31
        Bu4 = u1*dNdx21 + u2*dNdx11
        Bu5 = u2*dNdx31 + u3*dNdx21
        Bu6 = u3*dNdx11 + u1*dNdx31
        Bu1 = Bu1 + u4*dNdx12
        Bu2 = Bu2 + u5*dNdx22
        Bu3 = Bu3 + u6*dNdx32
        Bu4 = Bu4 + u4*dNdx22 + u5*dNdx12
        Bu5 = Bu5 + u5*dNdx32 + u6*dNdx22
        Bu6 = Bu6 + u6*dNdx12 + u4*dNdx32
        Bu1 = Bu1 + u7*dNdx13
        Bu2 = Bu2 + u8*dNdx23
        Bu3 = Bu3 + u9*dNdx33
        Bu4 = Bu4 + u7*dNdx23 + u8*dNdx13
        Bu5 = Bu5 + u8*dNdx33 + u9*dNdx23
        Bu6 = Bu6 + u9*dNdx13 + u7*dNdx33
        Bu1 = Bu1 + u10*dNdx14
        Bu2 = Bu2 + u11*dNdx24
        Bu3 = Bu3 + u12*dNdx34
        Bu4 = Bu4 + u10*dNdx24 + u11*dNdx14
        Bu5 = Bu5 + u11*dNdx34 + u12*dNdx24
        Bu6 = Bu6 + u12*dNdx14 + u10*dNdx34
        Bu1 = Bu1 + u13*dNdx15
        Bu2 = Bu2 + u14*dNdx25
        Bu3 = Bu3 + u15*dNdx35
        Bu4 = Bu4 + u13*dNdx25 + u14*dNdx15
        Bu5 = Bu5 + u14*dNdx35 + u15*dNdx25
        Bu6 = Bu6 + u15*dNdx15 + u13*dNdx35
        Bu1 = Bu1 + u16*dNdx16
        Bu2 = Bu2 + u17*dNdx26
        Bu3 = Bu3 + u18*dNdx36
        Bu4 = Bu4 + u16*dNdx26 + u17*dNdx16
        Bu5 = Bu5 + u17*dNdx36 + u18*dNdx26
        Bu6 = Bu6 + u18*dNdx16 + u16*dNdx36
        Bu1 = Bu1 + u19*dNdx17
        Bu2 = Bu2 + u20*dNdx27
        Bu3 = Bu3 + u21*dNdx37
        Bu4 = Bu4 + u19*dNdx27 + u20*dNdx17
        Bu5 = Bu5 + u20*dNdx37 + u21*dNdx27
        Bu6 = Bu6 + u21*dNdx17 + u19*dNdx37
        Bu1 = Bu1 + u22*dNdx18
        Bu2 = Bu2 + u23*dNdx28
        Bu3 = Bu3 + u24*dNdx38
        Bu4 = Bu4 + u22*dNdx28 + u23*dNdx18
        Bu5 = Bu5 + u23*dNdx38 + u24*dNdx28
        Bu6 = Bu6 + u24*dNdx18 + u22*dNdx38
        Bu1 = Bu1 + u25*dNdx19
        Bu2 = Bu2 + u26*dNdx29
        Bu3 = Bu3 + u27*dNdx39
        Bu4 = Bu4 + u25*dNdx29 + u26*dNdx19
        Bu5 = Bu5 + u26*dNdx39 + u27*dNdx29
        Bu6 = Bu6 + u27*dNdx19 + u25*dNdx39
        Bu1 = Bu1 + u28*dNdx1a
        Bu2 = Bu2 + u29*dNdx2a
        Bu3 = Bu3 + u30*dNdx3a
        Bu4 = Bu4 + u28*dNdx2a + u29*dNdx1a
        Bu5 = Bu5 + u29*dNdx3a + u30*dNdx2a
        Bu6 = Bu6 + u30*dNdx1a + u28*dNdx3a
        
        DBu1 = ms_Dmat(ivec,1,1,iip,ie)*Bu1+ms_Dmat(ivec,1,2,iip,ie)*Bu2
     &        +ms_Dmat(ivec,1,3,iip,ie)*Bu3+ms_Dmat(ivec,1,4,iip,ie)*Bu4
     &        +ms_Dmat(ivec,1,5,iip,ie)*Bu5+ms_Dmat(ivec,1,6,iip,ie)*Bu6
        DBu2 = ms_Dmat(ivec,2,1,iip,ie)*Bu1+ms_Dmat(ivec,2,2,iip,ie)*Bu2
     &        +ms_Dmat(ivec,2,3,iip,ie)*Bu3+ms_Dmat(ivec,2,4,iip,ie)*Bu4
     &        +ms_Dmat(ivec,2,5,iip,ie)*Bu5+ms_Dmat(ivec,2,6,iip,ie)*Bu6
        DBu3 = ms_Dmat(ivec,3,1,iip,ie)*Bu1+ms_Dmat(ivec,3,2,iip,ie)*Bu2
     &        +ms_Dmat(ivec,3,3,iip,ie)*Bu3+ms_Dmat(ivec,3,4,iip,ie)*Bu4
     &        +ms_Dmat(ivec,3,5,iip,ie)*Bu5+ms_Dmat(ivec,3,6,iip,ie)*Bu6
        DBu4 = ms_Dmat(ivec,4,1,iip,ie)*Bu1+ms_Dmat(ivec,4,2,iip,ie)*Bu2
     &        +ms_Dmat(ivec,4,3,iip,ie)*Bu3+ms_Dmat(ivec,4,4,iip,ie)*Bu4
     &        +ms_Dmat(ivec,4,5,iip,ie)*Bu5+ms_Dmat(ivec,4,6,iip,ie)*Bu6
        DBu5 = ms_Dmat(ivec,5,1,iip,ie)*Bu1+ms_Dmat(ivec,5,2,iip,ie)*Bu2
     &        +ms_Dmat(ivec,5,3,iip,ie)*Bu3+ms_Dmat(ivec,5,4,iip,ie)*Bu4
     &        +ms_Dmat(ivec,5,5,iip,ie)*Bu5+ms_Dmat(ivec,5,6,iip,ie)*Bu6
        DBu6 = ms_Dmat(ivec,6,1,iip,ie)*Bu1+ms_Dmat(ivec,6,2,iip,ie)*Bu2
     &        +ms_Dmat(ivec,6,3,iip,ie)*Bu3+ms_Dmat(ivec,6,4,iip,ie)*Bu4
     &        +ms_Dmat(ivec,6,5,iip,ie)*Bu5+ms_Dmat(ivec,6,6,iip,ie)*Bu6

        BDBu_iip1 = dNdx11*DBu1 + dNdx21*DBu4 + dNdx31*DBu6
        BDBu_iip2 = dNdx11*DBu4 + dNdx21*DBu2 + dNdx31*DBu5
        BDBu_iip3 = dNdx11*DBu6 + dNdx21*DBu5 + dNdx31*DBu3
        BDBu_iip4 = dNdx12*DBu1 + dNdx22*DBu4 + dNdx32*DBu6
        BDBu_iip5 = dNdx12*DBu4 + dNdx22*DBu2 + dNdx32*DBu5
        BDBu_iip6 = dNdx12*DBu6 + dNdx22*DBu5 + dNdx32*DBu3
        BDBu_iip7 = dNdx13*DBu1 + dNdx23*DBu4 + dNdx33*DBu6
        BDBu_iip8 = dNdx13*DBu4 + dNdx23*DBu2 + dNdx33*DBu5
        BDBu_iip9 = dNdx13*DBu6 + dNdx23*DBu5 + dNdx33*DBu3
        BDBu_iip10 = dNdx14*DBu1 + dNdx24*DBu4 + dNdx34*DBu6
        BDBu_iip11 = dNdx14*DBu4 + dNdx24*DBu2 + dNdx34*DBu5
        BDBu_iip12 = dNdx14*DBu6 + dNdx24*DBu5 + dNdx34*DBu3
        BDBu_iip13 = dNdx15*DBu1 + dNdx25*DBu4 + dNdx35*DBu6
        BDBu_iip14 = dNdx15*DBu4 + dNdx25*DBu2 + dNdx35*DBu5
        BDBu_iip15 = dNdx15*DBu6 + dNdx25*DBu5 + dNdx35*DBu3
        BDBu_iip16 = dNdx16*DBu1 + dNdx26*DBu4 + dNdx36*DBu6
        BDBu_iip17 = dNdx16*DBu4 + dNdx26*DBu2 + dNdx36*DBu5
        BDBu_iip18 = dNdx16*DBu6 + dNdx26*DBu5 + dNdx36*DBu3
        BDBu_iip19 = dNdx17*DBu1 + dNdx27*DBu4 + dNdx37*DBu6
        BDBu_iip20 = dNdx17*DBu4 + dNdx27*DBu2 + dNdx37*DBu5
        BDBu_iip21 = dNdx17*DBu6 + dNdx27*DBu5 + dNdx37*DBu3
        BDBu_iip22 = dNdx18*DBu1 + dNdx28*DBu4 + dNdx38*DBu6
        BDBu_iip23 = dNdx18*DBu4 + dNdx28*DBu2 + dNdx38*DBu5
        BDBu_iip24 = dNdx18*DBu6 + dNdx28*DBu5 + dNdx38*DBu3
        BDBu_iip25 = dNdx19*DBu1 + dNdx29*DBu4 + dNdx39*DBu6
        BDBu_iip26 = dNdx19*DBu4 + dNdx29*DBu2 + dNdx39*DBu5
        BDBu_iip27 = dNdx19*DBu6 + dNdx29*DBu5 + dNdx39*DBu3
        BDBu_iip28 = dNdx1a*DBu1 + dNdx2a*DBu4 + dNdx3a*DBu6
        BDBu_iip29 = dNdx1a*DBu4 + dNdx2a*DBu2 + dNdx3a*DBu5
        BDBu_iip30 = dNdx1a*DBu6 + dNdx2a*DBu5 + dNdx3a*DBu3
        BDBu1 = BDBu1 + weight(iip)*BDBu_iip1
        BDBu2 = BDBu2 + weight(iip)*BDBu_iip2
        BDBu3 = BDBu3 + weight(iip)*BDBu_iip3
        BDBu4 = BDBu4 + weight(iip)*BDBu_iip4
        BDBu5 = BDBu5 + weight(iip)*BDBu_iip5
        BDBu6 = BDBu6 + weight(iip)*BDBu_iip6
        BDBu7 = BDBu7 + weight(iip)*BDBu_iip7
        BDBu8 = BDBu8 + weight(iip)*BDBu_iip8
        BDBu9 = BDBu9 + weight(iip)*BDBu_iip9
        BDBu10 = BDBu10 + weight(iip)*BDBu_iip10
        BDBu11 = BDBu11 + weight(iip)*BDBu_iip11
        BDBu12 = BDBu12 + weight(iip)*BDBu_iip12
        BDBu13 = BDBu13 + weight(iip)*BDBu_iip13
        BDBu14 = BDBu14 + weight(iip)*BDBu_iip14
        BDBu15 = BDBu15 + weight(iip)*BDBu_iip15
        BDBu16 = BDBu16 + weight(iip)*BDBu_iip16
        BDBu17 = BDBu17 + weight(iip)*BDBu_iip17
        BDBu18 = BDBu18 + weight(iip)*BDBu_iip18
        BDBu19 = BDBu19 + weight(iip)*BDBu_iip19
        BDBu20 = BDBu20 + weight(iip)*BDBu_iip20
        BDBu21 = BDBu21 + weight(iip)*BDBu_iip21
        BDBu22 = BDBu22 + weight(iip)*BDBu_iip22
        BDBu23 = BDBu23 + weight(iip)*BDBu_iip23
        BDBu24 = BDBu24 + weight(iip)*BDBu_iip24
        BDBu25 = BDBu25 + weight(iip)*BDBu_iip25
        BDBu26 = BDBu26 + weight(iip)*BDBu_iip26
        BDBu27 = BDBu27 + weight(iip)*BDBu_iip27
        BDBu28 = BDBu28 + weight(iip)*BDBu_iip28
        BDBu29 = BDBu29 + weight(iip)*BDBu_iip29
        BDBu30 = BDBu30 + weight(iip)*BDBu_iip30
      enddo ! iip

      BDBu1 = BDBu1*detj
      BDBu2 = BDBu2*detj
      BDBu3 = BDBu3*detj
      BDBu4 = BDBu4*detj
      BDBu5 = BDBu5*detj
      BDBu6 = BDBu6*detj
      BDBu7 = BDBu7*detj
      BDBu8 = BDBu8*detj
      BDBu9 = BDBu9*detj
      BDBu10 = BDBu10*detj
      BDBu11 = BDBu11*detj
      BDBu12 = BDBu12*detj
      BDBu13 = BDBu13*detj
      BDBu14 = BDBu14*detj
      BDBu15 = BDBu15*detj
      BDBu16 = BDBu16*detj
      BDBu17 = BDBu17*detj
      BDBu18 = BDBu18*detj
      BDBu19 = BDBu19*detj
      BDBu20 = BDBu20*detj
      BDBu21 = BDBu21*detj
      BDBu22 = BDBu22*detj
      BDBu23 = BDBu23*detj
      BDBu24 = BDBu24*detj
      BDBu25 = BDBu25*detj
      BDBu26 = BDBu26*detj
      BDBu27 = BDBu27*detj
      BDBu28 = BDBu28*detj
      BDBu29 = BDBu29*detj
      BDBu30 = BDBu30*detj
      tmp1i=tmp1(ivec,ie)
      tmp2i=tmp2(ivec,ie)
      BDBu1=(
     - m00238*u1+
     - m00039*(u4+u7+u10)+
     - mm0158*(u13+u19+u22)+
     - mm0238*(u16+u25+u28)
     - )*tmp1i+BDBu1*tmp2i
      BDBu4=(
     - m00039*(u1+u7+u10)+
     - m00238*u4+
     - mm0158*(u13+u16+u25)+
     - mm0238*(u19+u22+u28)
     - )*tmp1i+BDBu4*tmp2i
      BDBu7=(
     - m00039*(u1+u4+u10)+
     - m00238*u7+
     - mm0238*(u13+u22+u25)+
     - mm0158*(u16+u19+u28)
     - )*tmp1i+BDBu7*tmp2i
      BDBu10=(
     - m00039*(u1+u4+u7)+
     - m00238*u10+
     - mm0238*(u13+u16+u19)+
     - mm0158*(u22+u25+u28)
     - )*tmp1i+BDBu10*tmp2i
      BDBu13=(
     - mm0158*(u1+u4)+
     - mm0238*(u7+u10)+
     - m01269*u13+
     - m00634*(u16+u19+u22+u25)+
     - m00317*u28
     - )*tmp1i+BDBu13*tmp2i
      BDBu16=(
     - mm0238*(u1+u10)+
     - mm0158*(u4+u7)+
     - m00634*(u13+u19+u25+u28)+
     - m01269*u16+
     - m00317*u22
     - )*tmp1i+BDBu16*tmp2i
      BDBu19=(
     - mm0158*(u1+u7)+
     - mm0238*(u4+u10)+
     - m00634*(u13+u16+u22+u28)+
     - m01269*u19+
     - m00317*u25
     - )*tmp1i+BDBu19*tmp2i
      BDBu22=(
     - mm0158*(u1+u10)+
     - mm0238*(u4+u7)+
     - m00634*(u13+u19+u25+u28)+
     - m00317*u16+
     - m01269*u22
     - )*tmp1i+BDBu22*tmp2i
      BDBu25=(
     - mm0238*(u1+u7)+
     - mm0158*(u4+u10)+
     - m00634*(u13+u16+u22+u28)+
     - m00317*u19+
     - m01269*u25
     - )*tmp1i+BDBu25*tmp2i
      BDBu28=(
     - mm0238*(u1+u4)+
     - mm0158*(u7+u10)+
     - m00317*u13+
     - m00634*(u16+u19+u22+u25)+
     - m01269*u28
     - )*tmp1i+BDBu28*tmp2i
      BDBu2=(
     - m00238*u2+
     - m00039*(u5+u8+u11)+
     - mm0158*(u14+u20+u23)+
     - mm0238*(u17+u26+u29)
     - )*tmp1i+BDBu2*tmp2i
      BDBu5=(
     - m00039*(u2+u8+u11)+
     - m00238*u5+
     - mm0158*(u14+u17+u26)+
     - mm0238*(u20+u23+u29)
     - )*tmp1i+BDBu5*tmp2i
      BDBu8=(
     - m00039*(u2+u5+u11)+
     - m00238*u8+
     - mm0238*(u14+u23+u26)+
     - mm0158*(u17+u20+u29)
     - )*tmp1i+BDBu8*tmp2i
      BDBu11=(
     - m00039*(u2+u5+u8)+
     - m00238*u11+
     - mm0238*(u14+u17+u20)+
     - mm0158*(u23+u26+u29)
     - )*tmp1i+BDBu11*tmp2i
      BDBu14=(
     - mm0158*(u2+u5)+
     - mm0238*(u8+u11)+
     - m01269*u14+
     - m00634*(u17+u20+u23+u26)+
     - m00317*u29
     - )*tmp1i+BDBu14*tmp2i
      BDBu17=(
     - mm0238*(u2+u11)+
     - mm0158*(u5+u8)+
     - m00634*(u14+u20+u26+u29)+
     - m01269*u17+
     - m00317*u23
     - )*tmp1i+BDBu17*tmp2i
      BDBu20=(
     - mm0158*(u2+u8)+
     - mm0238*(u5+u11)+
     - m00634*(u14+u17+u23+u29)+
     - m01269*u20+
     - m00317*u26
     - )*tmp1i+BDBu20*tmp2i
      BDBu23=(
     - mm0158*(u2+u11)+
     - mm0238*(u5+u8)+
     - m00634*(u14+u20+u26+u29)+
     - m00317*u17+
     - m01269*u23
     - )*tmp1i+BDBu23*tmp2i
      BDBu26=(
     - mm0238*(u2+u8)+
     - mm0158*(u5+u11)+
     - m00634*(u14+u17+u23+u29)+
     - m00317*u20+
     - m01269*u26
     - )*tmp1i+BDBu26*tmp2i
      BDBu29=(
     - mm0238*(u2+u5)+
     - mm0158*(u8+u11)+
     - m00317*u14+
     - m00634*(u17+u20+u23+u26)+
     - m01269*u29
     - )*tmp1i+BDBu29*tmp2i
      BDBu3=(
     - m00238*u3+
     - m00039*(u6+u9+u12)+
     - mm0158*(u15+u21+u24)+
     - mm0238*(u18+u27+u30)
     - )*tmp1i+BDBu3*tmp2i
      BDBu6=(
     - m00039*(u3+u9+u12)+
     - m00238*u6+
     - mm0158*(u15+u18+u27)+
     - mm0238*(u21+u24+u30)
     - )*tmp1i+BDBu6*tmp2i
      BDBu9=(
     - m00039*(u3+u6+u12)+
     - m00238*u9+
     - mm0238*(u15+u24+u27)+
     - mm0158*(u18+u21+u30)
     - )*tmp1i+BDBu9*tmp2i
      BDBu12=(
     - m00039*(u3+u6+u9)+
     - m00238*u12+
     - mm0238*(u15+u18+u21)+
     - mm0158*(u24+u27+u30)
     - )*tmp1i+BDBu12*tmp2i
      BDBu15=(
     - mm0158*(u3+u6)+
     - mm0238*(u9+u12)+
     - m01269*u15+
     - m00634*(u18+u21+u24+u27)+
     - m00317*u30
     - )*tmp1i+BDBu15*tmp2i
      BDBu18=(
     - mm0238*(u3+u12)+
     - mm0158*(u6+u9)+
     - m00634*(u15+u21+u27+u30)+
     - m01269*u18+
     - m00317*u24
     - )*tmp1i+BDBu18*tmp2i
      BDBu21=(
     - mm0158*(u3+u9)+
     - mm0238*(u6+u12)+
     - m00634*(u15+u18+u24+u30)+
     - m01269*u21+
     - m00317*u27
     - )*tmp1i+BDBu21*tmp2i
      BDBu24=(
     - mm0158*(u3+u12)+
     - mm0238*(u6+u9)+
     - m00634*(u15+u21+u27+u30)+
     - m00317*u18+
     - m01269*u24
     - )*tmp1i+BDBu24*tmp2i
      BDBu27=(
     - mm0238*(u3+u9)+
     - mm0158*(u6+u12)+
     - m00634*(u15+u18+u24+u30)+
     - m00317*u21+
     - m01269*u27
     - )*tmp1i+BDBu27*tmp2i
      BDBu30=(
     - mm0238*(u3+u6)+
     - mm0158*(u9+u12)+
     - m00317*u15+
     - m00634*(u18+u21+u24+u27)+
     - m01269*u30
     - )*tmp1i+BDBu30*tmp2i
!$acc atomic
      rg(ivec,1,cny(1,ie))=rg(ivec,1,cny(1,ie))+BDBu1
!$acc atomic
      rg(ivec,2,cny(1,ie))=rg(ivec,2,cny(1,ie))+BDBu2
!$acc atomic
      rg(ivec,3,cny(1,ie))=rg(ivec,3,cny(1,ie))+BDBu3
!$acc atomic
      rg(ivec,1,cny(2,ie))=rg(ivec,1,cny(2,ie))+BDBu4
!$acc atomic
      rg(ivec,2,cny(2,ie))=rg(ivec,2,cny(2,ie))+BDBu5
!$acc atomic
      rg(ivec,3,cny(2,ie))=rg(ivec,3,cny(2,ie))+BDBu6
!$acc atomic
      rg(ivec,1,cny(3,ie))=rg(ivec,1,cny(3,ie))+BDBu7
!$acc atomic
      rg(ivec,2,cny(3,ie))=rg(ivec,2,cny(3,ie))+BDBu8
!$acc atomic
      rg(ivec,3,cny(3,ie))=rg(ivec,3,cny(3,ie))+BDBu9
!$acc atomic
      rg(ivec,1,cny(4,ie))=rg(ivec,1,cny(4,ie))+BDBu10
!$acc atomic
      rg(ivec,2,cny(4,ie))=rg(ivec,2,cny(4,ie))+BDBu11
!$acc atomic
      rg(ivec,3,cny(4,ie))=rg(ivec,3,cny(4,ie))+BDBu12
!$acc atomic
      rg(ivec,1,cny(5,ie))=rg(ivec,1,cny(5,ie))+BDBu13
!$acc atomic
      rg(ivec,2,cny(5,ie))=rg(ivec,2,cny(5,ie))+BDBu14
!$acc atomic
      rg(ivec,3,cny(5,ie))=rg(ivec,3,cny(5,ie))+BDBu15
!$acc atomic
      rg(ivec,1,cny(6,ie))=rg(ivec,1,cny(6,ie))+BDBu16
!$acc atomic
      rg(ivec,2,cny(6,ie))=rg(ivec,2,cny(6,ie))+BDBu17
!$acc atomic
      rg(ivec,3,cny(6,ie))=rg(ivec,3,cny(6,ie))+BDBu18
!$acc atomic
      rg(ivec,1,cny(7,ie))=rg(ivec,1,cny(7,ie))+BDBu19
!$acc atomic
      rg(ivec,2,cny(7,ie))=rg(ivec,2,cny(7,ie))+BDBu20
!$acc atomic
      rg(ivec,3,cny(7,ie))=rg(ivec,3,cny(7,ie))+BDBu21
!$acc atomic
      rg(ivec,1,cny(8,ie))=rg(ivec,1,cny(8,ie))+BDBu22
!$acc atomic
      rg(ivec,2,cny(8,ie))=rg(ivec,2,cny(8,ie))+BDBu23
!$acc atomic
      rg(ivec,3,cny(8,ie))=rg(ivec,3,cny(8,ie))+BDBu24
!$acc atomic
      rg(ivec,1,cny(9,ie))=rg(ivec,1,cny(9,ie))+BDBu25
!$acc atomic
      rg(ivec,2,cny(9,ie))=rg(ivec,2,cny(9,ie))+BDBu26
!$acc atomic
      rg(ivec,3,cny(9,ie))=rg(ivec,3,cny(9,ie))+BDBu27
!$acc atomic
      rg(ivec,1,cny(10,ie))=rg(ivec,1,cny(10,ie))+BDBu28
!$acc atomic
      rg(ivec,2,cny(10,ie))=rg(ivec,2,cny(10,ie))+BDBu29
!$acc atomic
      rg(ivec,3,cny(10,ie))=rg(ivec,3,cny(10,ie))+BDBu30
      enddo
      enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
      end
