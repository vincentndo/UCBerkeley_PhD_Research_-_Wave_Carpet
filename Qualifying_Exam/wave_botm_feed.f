      PROGRAM WAVE_BOTM
c
c WAVE_BOTM.F -- program for wave-bottom interaction calculations in 3D. 
c This is the spectral code modified from DKY's resonant wave-wave 
c interactions which is the modification of DGR's moving pressure problem. 
c by YL 12/7/93, modified by ND 2015
c
c for CRAY programs:
c  -> remove programs after "CRAY REMOVE"
c  -> append FFT.F
c  -> comment out the call to OUTMODE
c  -> change the output unit IOUTTT from 6 to 98
c  -> change list directed output to fixed format
c
      IMPLICIT  REAL(A-H,M,O-Z)
C
C #WORDS ON CRAY IS ABOUT NDX*NDY*((18+4*NDP)/8+8+NDP*2)
C
      PARAMETER(NDX=512,NDY=4,NDP=1)
      INTEGER   NUMB(20,2) !,mode(2,ndx)
c      REAL      Y1(NDX,NDY,NDP),Y2(NDX,NDY,NDP),BAMP(NDX,NDY)
      REAL      Y1(NDX,NDY,NDP),Y2(NDX,NDY,NDP)
      REAL      FQ(NDX/2+1,NDY/4+1,6)
      REAL      F(NDX), G(NDX)
      COMPLEX   VP(NDX,NDY),SF(NDX,NDY),W(NDX,NDY)
      COMPLEX   VPB(NDX,NDY),BAMP(NDX,NDY),WB(NDX,NDY)
      COMPLEX   SF_input(10240,NDX,NDY),VP_input(10240,NDX,NDY)
      COMPLEX   VPB_input(10240,NDX,NDY),BAMP_input(10240,NDX,NDY)
      COMPLEX   R1(NDX/2+1,NDY/4+1),R2(NDX/2+1,NDY/4+1)
      COMPLEX   RB1(NDX/2+1,NDY/4+1),RB2(NDX/2+1,NDY/4+1)
      COMPLEX   VPOT(NDX/2+1,NDY/4+1,4),SURF(NDX/2+1,NDY/4+1,4)
      COMPLEX   VPOTB(NDX/2+1,NDY/4+1,4),SURFB(NDX/2+1,NDY/4+1,4)
      COMPLEX   SAVE(NDX/2+1,NDY/4+1,2),SAVEB(NDX/2+1,NDY/4+1,2)
      COMPLEX   P1(NDX/2+1,NDY/4+1,NDP),P2(NDX/2+1,NDY/4+1,NDP)
      COMPLEX   P3(NDX/2+1,NDY/4+1,NDP),P4(NDX/2+1,NDY/4+1,NDP)
      COMPLEX   P24(NDX/2+1,NDY/4+1,NDP)
      COMPLEX   VX(NDX),VY(NDX),WRKX((NDX/2)*5),WRKY((NDY/2)*5)
      COMPLEX   WRK1(NDX,NDY),WRK2(NDX,NDY),WRK3(NDX,NDY),WRK4(NDX,NDY)
      COMMON/ABLOCK/ISTP,DELT,TLEN,TWID,PERIOD,ENERGY0,DEPT,
     1              ZETA,GAMMA,MU
      COMMON/IBLOCK/NUMB 
C      common/inout/in1,infsph,iouttt,iouteng,iout1,iout2,ioutrs
C      data in1/10/,iouttt/6/,iout3/97/
C
C INITIALIZATION 
C
      TIME=0.
C      CALL INIT(NDX,NDY,NDP,SF,VP,BAMP,FQ,ISTP,itout,nmode,mode,
C     $          VX,VY,WRKX,WRKY)
c      CALL INIT(NDX,NDY,NDP,SF,VP,BAMP,FQ,ISTP,VX,VY,WRKX,WRKY)
      CALL INIT(NDX,NDY,NDP,SF,VP,BAMP,VPB,FQ,VX,VY,WRKX,WRKY)
c      WRITE(*,*) (REAL(BAMP(I,1)),I=1,3)
c      CALL BOTM(NDX,NDY,NDP,BAMP,Y2,WRK1,VX,VY,WRKX,WRKY)
      
      call fg(ndx,ndy,ndp,f,g)      
      call treat(ndx,ndy,ndp,sf,vp,bamp,vpb,f)      
      call readInput(ndx,ndy,ndp,sf_input,vp_input,bamp_input,vpb_input)
      
C
C IMPLEMENT THE 4TH ORDER RUNGE-KUTTA ALGORITHM
C
      DO 200 ITIM=1,ISTP
        TPERIOD=TIME/PERIOD
        WRITE(*,*) ' -> time step, time* ',ITIM,TPERIOD 
        DO 100 IRUN=1,4
          IF(IRUN.EQ.2.OR.IRUN.EQ.4) TIME=TIME+.5*DELT
          CALL SUMA(NDX,NDY,IRUN,SURF,VPOT,SAVE,SF,VP)
c          IF(IRUN.EQ.3) THEN
c          WRITE(*,*) (REAL(SURF(I,2,1)),I=1,3)
c          WRITE(*,*) (REAL(SURF(I,2,2)),I=1,3)
c          WRITE(*,*) (REAL(SURF(I,2,3)),I=1,3)
c          WRITE(*,*) (REAL(SURF(I,2,4)),I=1,3)
c          ENDIF
          CALL SUMA(NDX,NDY,IRUN,SURFB,VPOTB,SAVEB,BAMP,VPB)
c          IF(IRUN.EQ.1) THEN
c          WRITE(*,*) (REAL(BAMP(I,1)),I=1,3)          
c          ENDIF
          CALL FREX(NDX,NDY,NDP,SF,Y1,WRK1,VX,VY,WRKX,WRKY)
c          IF(IRUN.EQ.2) THEN
c          WRITE(*,*) (REAL(Y1(I,2,1)),I=1,3)          
c          ENDIF
c          CALL BOTM(NDX,NDY,NDP,BAMP,Y2,WRK1,VX,VY,WRKX,WRKY)
          CALL FREX(NDX,NDY,NDP,BAMP,Y2,WRK1,VX,VY,WRKX,WRKY)
c          WRITE(*,*) (REAL(BAMP(I,1)),I=1,3)
c          CALL  BVP(NDX,NDY,NDP,VP,P1,P2,P3,P4,Y1,Y2,FQ,
          CALL  BVP(NDX,NDY,NDP,VP,VPB,P1,P2,P24,P3,P4,Y1,Y2,FQ,
     1              WRK1,WRK2,WRK3,WRK4,VX,VY,WRKX,WRKY)
c          IF(IRUN.EQ.1) THEN
c          WRITE(*,*) (REAL(P3(I,2,1)),I=1,3)          
c          ENDIF
          CALL VELO(NDX,NDY,NDP,P1,P3,Y1,W,FQ,WRK1,VX,VY,WRKX,WRKY)
c          WRITE(*,*) (REAL(W(I,1)),I=1,3)
c          IF(IRUN.EQ.2) THEN
c          WRITE(*,*) (REAL(W(I,2)),I=1,3)          
c          ENDIF
          CALL VELO(NDX,NDY,NDP,P24,P4,Y2,WB,FQ,WRK1,VX,VY,WRKX,WRKY)
c          WRITE(*,*) (REAL(P24(I,1,1)),I=1,3)
          CALL RIGH(NDX,NDY,NDP,R1,R2,RB1,RB2,SF,VP,W,BAMP,VPB,WB,FQ,
     1              WRK1,WRK2,WRK3,WRK4,VX,VY,WRKX,WRKY)
c          WRITE(*,*) (REAL(W(I,1)),I=1,3)
c          IF(IRUN.EQ.2) THEN
c          WRITE(*,*) (REAL(R1(I,2)),I=1,3)
c          WRITE(*,*) (REAL(R2(I,2)),I=1,3)          
c          ENDIF
          CALL UPFR(NDX,NDY,IRUN,SURF,VPOT,R1,R2)
          CALL UPFR(NDX,NDY,IRUN,SURFB,VPOTB,RB1,RB2)
c          WRITE(*,*) (REAL(RB1(I,1)),I=1,3)
c          WRITE(*,*) (REAL(SURFB(I,1,1)),I=1,3)
c          WRITE(*,*) (REAL(SURFB(I,1,2)),I=1,3)
c          WRITE(*,*) (REAL(SURFB(I,1,3)),I=1,3)
c          WRITE(*,*) (REAL(SURFB(I,1,4)),I=1,3)
c          IF(IRUN.EQ.1) THEN
c          WRITE(*,*) (REAL(R1(I,2)),I=1,3)
c          WRITE(*,*) (REAL(R2(I,2)),I=1,3)
c          WRITE(*,*) (REAL(RB1(I,2)),I=1,3)
c          WRITE(*,*) (REAL(RB2(I,2)),I=1,3)          
c          ENDIF
          IF(IRUN.EQ.1)THEN     ! to do every integral time step
          CALL CHCK(NDX,NDY,ITIM,TPERIOD,SF,VP,R1,WRK1,WRK2,WRK3,
     1              VX,VY,WRKX,WRKY)
          CALL WRIT(NDX,NDY,TIME,SF,VP,BAMP,VX,VY,WRKX,WRKY)
C          CALL output1(iout1,NDX,NDY,tperiod,SF,nmode,mode)
C          if(mod(itim-1,itout).eq.0) then    ! to do every ITOUT steps
C            CALL output2(iout2,NDX,NDY,tperiod,SF,wrk1,
C     $                     vx,vy,wrkx,wrky)
C           call outmode(iout3,NDX,NDY,time,SF)
C          endif
          ENDIF
100     CONTINUE
        CALL SUMA(NDX,NDY,5,SURF,VPOT,SAVE,SF,VP)
        CALL SUMA(NDX,NDY,5,SURFB,VPOTB,SAVEB,BAMP,VPB)
        CALL SMOT1(NDX,NDY,SF,VP,FQ,0.9)
        CALL SMOT1(NDX,NDY,BAMP,VPB,FQ,0.9)
        call treatment(ndx,ndy,ndp,sf,vp,bamp,vpb,sf_input,vp_input,
     1                 bamp_input,vpb_input,f,g,itim)
200   CONTINUE
C      close(iout1)
C      close(iout2)
      close(11)
      close(12)
      STOP
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE output1(iout1,NDX,NDY,TIME,SF,nmode,mode)
C
C output at every time step
C
      IMPLICIT  REAL(A-H,M,O-Z)
      INTEGER   mode(2,nmode)
      INTEGER   m
      COMPLEX   SF(NDX,NDY)
C
      WRITE(iout1,*) time,(cabs(sf(mode(1,m),mode(2,m))),m=1,nmode)
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      subroutine output2(iout2,NDX,NDY,time,SF,wrk1,vx,vy,wrkx,wrky)
C
C output for every ITOUT time steps
c ----> sf in wavenumber and physical spaces
C
      IMPLICIT REAL(A-H,M,O-Z)
      INTEGER  NUMB(20,2)
      complex  sf(ndx,ndy),wrk1(ndx,ndy)
      complex  vx(ndx),vy(ndx),wrkx((ndx/2)*5),wrky((ndy/2)*5)
      COMMON/IBLOCK/NUMB 
C
      write(iout2,*) time
C
C sf in wavenumber space
C
      do 10 j=1,numb(4,2)
        write(iout2,*) (sf(i,j),i=1,numb(4,1)),
     $                 (sf(i,j),i=numb(7,1),numb(8,1))
10    continue
C
C transform fs from wavenumber space to physical space
C
      do 20 j=1,numb(8,2)
      do 20 i=1,numb(8,1)
        wrk1(i,j)=sf(i,j)
20    continue
      call fft(ndx,ndy,0,-1,wrk1,vx,vy,wrkx,wrky)
      do 30 j=1,numb(8,2) 
30    write(iout2,*) (real(wrk1(i,j)),i=1,numb(8,1))
c
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE WRIT(NDX,NDY,TIME,SF,VP,BAMP,VX,VY,WRKX,WRKY)
C
C OUTPUT THE RESULTS
C
      IMPLICIT  REAL(A-H,M,O-Z)
      INTEGER   NUMB(20,2)
      COMPLEX   SF(NDX,NDY),VP(NDX,NDY),BAMP(NDX,NDY)
      COMPLEX   SF_temp(NDX,NDY),VP_temp(NDX,NDY)
      COMPLEX   BAMP_temp(NDX,NDY)
      COMPLEX   VX(NDX),VY(NDX),WRKX((NDX/2)*5),WRKY((NDY/2)*5)
      COMMON/IBLOCK/NUMB 
C
C WRITE INTERMEDIATE RESULTS
C
      DO 100 J=1,NUMB(4,2)
100   WRITE(11,*) TIME,(SF(I,J),I=1,NUMB(4,1)),
     1               (SF(I,J),I=NUMB(7,1),NUMB(8,1))
c 100  WRITE(11,*) TIME,(VP(I,J),I=1,NUMB(4,1)),
c     1               (VP(I,J),I=NUMB(7,1),NUMB(8,1))

      DO 106 J=1,NUMB(8,2)
      DO 106 I=1,NUMB(8,1)
      SF_temp(i,j) = SF(i,j)
      VP_temp(i,j) = VP(i,j)
      BAMP_temp(i,j) = BAMP(i,j)
106   CONTINUE      
      
      CALL FFT(ndx,ndy,0,-1,SF_temp,vx,vy,wrkx,wrky)
      CALL FFT(ndx,ndy,0,-1,BAMP_temp,vx,vy,wrkx,wrky)
      
c      DO 108 J=1,5
c      WRITE(*,*) (REAL(SF_temp(I,J)),I=1,3)
c108   WRITE(*,*) (REAL(BAMP_temp(I,J)),I=1,3)
      
c      WRITE(*,*) "MILESTONE"
      
      DO 107 J=1,NUMB(8,2)
      WRITE(25,*) (REAL(SF_temp(I,J)),I=1,NUMB(8,1))
107   WRITE(25,*) (REAL(BAMP_temp(I,J)),I=1,NUMB(8,1))    
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE CHCK(NDX,NDY,ITIM,TIME,SF,VP,R1,WRK1,WRK2,WRK3,
     1            VX,VY,WRKX,WRKY)
C
C FIND THE POTENTIAL AND KINETIC ENERGIES
C
      IMPLICIT REAL(A-H,M,O-Z)
      REAL     KINE,POTL,ENER,MU
      INTEGER  NUMB(20,2)
      COMPLEX  SF(NDX,NDY),VP(NDX,NDY),R1(NDX/2+1,NDY/4+1)
      COMPLEX  WRK1(NDX,NDY),WRK2(NDX,NDY),WRK3(NDX,NDY)
      complex  vx(ndx),vy(ndx),wrkx((ndx/2)*5),wrky((ndy/2)*5)
      COMMON/ABLOCK/ISTP,DELT,TLEN,TWID,PERIOD,ENERGY0,DEPT,
     1              ZETA,GAMMA,MU
      COMMON/IBLOCK/NUMB 
C      common/inout/in1,infsph,iouttt,iouteng,iout1,iout2,ioutrs
C
C DO VOLUME CONSERVATION
C
      VOLU=TLEN*TWID*REAL(SF(1,1))
C
C DO VOLUME FLUX
C
      FLUX=TLEN*TWID*REAL(R1(1,1))
C
C TRANSFORM FROM WAVENUMBER SPACE TO PHYSICAL SPACE
C
      DO 10 J=1,NUMB(8,2)
      DO 10 I=1,NUMB(8,1)
      WRK1(I,J)=SF(I,J)
10    WRK2(I,J)=VP(I,J)
C
C copy/unpack the compact R1
C
      CALL UN_PACK(NDX,NDY,R1,WRK3)
C
      CALL FFT(NDX,NDY,0,-1,WRK1,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WRK2,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WRK3,VX,VY,WRKX,WRKY)
C
C DO POTENTIAL ENERGY AND KINETIC ENERGIES IN PHYSICAL SPACE
C
      DO 20 J=1,NUMB(8,2)
      DO 20 I=1,NUMB(8,1)
      WRK1(I,J)=WRK1(I,J)*WRK1(I,J)
20    WRK2(I,J)=WRK2(I,J)*WRK3(I,J)
C
C TRANSFORM FROM PHYSICAL SPACE TO WAVENUMBER SPACE
C
      CALL FFT(NDX,NDY,0,1,WRK1,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,WRK2,VX,VY,WRKX,WRKY)
C
C EVALUATE POTENTIAL, KINETIC, AND TOTAL ENERGIES
C
      POTL=0.5*TLEN*TWID*REAL(WRK1(1,1))
      KINE=0.5*TLEN*TWID*REAL(WRK2(1,1))
      ENER=POTL+KINE
      IF(ITIM.EQ.1) energy0=ENER
      WRITE(12,900) TIME,VOLU,FLUX,POTL,KINE,ENER,ener/energy0
C      write(iouttt,910) potl,kine,ener,ener/energy0
900   FORMAT(7F11.16)
910   format(5x,'pe,ke,e,e/e0=',4F12.16)
C
c      IF(ABS(energy0-ENER).GE.0.1*energy0) STOP 'bad ener'
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
C      SUBROUTINE INIT(NDX,NDY,NDP,SF,VP,BAMP,FQ,ISTP,itout,nmode,mode,
C     1                VX,VY,WRKX,WRKY)
c      SUBROUTINE INIT(NDX,NDY,NDP,SF,VP,BAMP,FQ,ISTP,VX,VY,WRKX,WRKY)
      SUBROUTINE INIT(NDX,NDY,NDP,SF,VP,BAMP,VPB,FQ,
     1                 VX,VY,WRKX,WRKY)
C
C INITIALIZE DATA
C
      IMPLICIT  REAL(A-H,M,O-Z)
      INTEGER   NUMB(20,2) !,mode(2,ndx)
c      REAL      BAMP(NDX,NDY),FQ(NDX/2+1,NDY/4+1,6)
      REAL      FQ(NDX/2+1,NDY/4+1,6),MU
c      COMPLEX   SF(NDX,NDY),VP(NDX,NDY)
      COMPLEX   SF(NDX,NDY),VP(NDX,NDY),BAMP(NDX,NDY),VPB(NDX,NDY)
      COMPLEX   SF_temp(NDX,NDY),VP_temp(NDX,NDY)
      COMPLEX   VX(NDX),VY(NDX),WRKX((NDX/2)*5),WRKY((NDY/2)*5)
      COMMON/ABLOCK/ISTP,DELT,TLEN,TWID,PERIOD,ENERGY0,DEPT,
     1              ZETA,GAMMA,MU
      COMMON/IBLOCK/NUMB 
C      common/inout/in1,infsph,iouttt,iouteng,iout1,iout2,ioutrs
      DATA TWOPI/6.283185307179586476925287/
C
C NUMBER OF PERTURBATION TERMS, AND X AND Y aliased FOURIER COEFFICIENTS
C
C      read(in1,*) infsph,iouteng,iout1,iout2,ioutrs
C
      READ(10,*) numb(1,1),Np,Nq                   ! M,p,q
      NUMB(2,1)=Np+1
      NUMB(2,2)=Nq+1
c
c for P==2**p, Q==2**q
c complex alised modes:	  sum_p(-2P,2P) sum_q(-2Q,2Q)= (4P+1)(4Q+1)
c   storage required = [NDX=4P]x[NDY=4Q] complex
c complex unalised modes: sum_p(-P,P) sum_q(-Q,Q)    = (2P+1)(2Q+1)
c
C TIME STEPS PER WAVE PERIOD AND NUMBER OF TIME STEPS
c output frequency: write output when mod(itim-1,itout)=0
C
      READ(10,*) IPER,ISTP !,itout
C
C input modes of interest to output to ioutfs and ioutph
C
C      read(in1,*) nmode,((mode(k,imode),k=1,2),imode=1,nmode)
C      close(in1)
C      if(nmode.gt.ndx) stop 'nmode>nx'
C
C ASSIGN POINTERS FOR X
C
      NUMB(4,1)= 2**(NUMB(2,1)-1)+1       ! P+1
      NUMB(5,1)= 2**(NUMB(2,1)-1)+2       ! P+2
      NUMB(6,1)= 2**(NUMB(2,1)-1)*3       ! 3P
      NUMB(7,1)= 2**(NUMB(2,1)-1)*3+1     ! 3P+1
      NUMB(8,1)= 2**(NUMB(2,1)+1)         ! 4P
      NUMB(9,1)= 2**(NUMB(2,1)+1)+1       ! 4P+1
      NUMB(10,1)=2**(NUMB(2,1))-1         ! 2P-1
C
C ASSIGN POINTERS FOR Y
C
      NUMB(4,2)= 2**(NUMB(2,2)-1)+1       ! Q+1
      NUMB(5,2)= 2**(NUMB(2,2)-1)+2       ! Q+2
      NUMB(6,2)= 2**(NUMB(2,2))+1         ! 2Q+1
      NUMB(8,2)= 2**(NUMB(2,2)+1)         ! 4Q
      NUMB(9,2)= 2**(NUMB(2,2)+1)+1       ! 4Q+1
C
C CHECK FOR DIMENSIONS 
C
      if(numb(1,1).ne.ndp) stop 'bad NDP'
      if(numb(8,1).ne.ndx) stop 'bad NDX'
      if(numb(8,2).ne.ndy) stop 'bad NDY'
      if(numb(8,1).lt.numb(8,2)) stop 'ndx<ndy'
c
c change mode(2,nmode)to pointers to the wavenumber arrays
C
C      mp=2**np
C      mq=2**nq
C      do 70 m=1,nmode
C        if(iabs(mode(1,m)).gt.mp) stop 'bad mode'
C        if(iabs(mode(2,m)).gt.mq) stop 'bad mode'
C choose the complex conjugate if y-mode < 0
C        if(mode(2,m).lt.0) then
C          mode(1,m)=-mode(1,m)
C          mode(2,m)=-mode(2,m)
C        endif          
C        if(mode(1,m).ge.0) then
C          mode(1,m)=mode(1,m)+1
C        else
C          mode(1,m)=numb(9,1)+mode(1,m)
C        endif
C        mode(2,m)=mode(2,m)+1
C70    continue
C
C INPUT DIMENSIONS OF THE COMPUTATIONAL DOMAIN, FUNDAMENTAL WAVE PERIOD, 
C AND WATER DEPTH. 
C
      READ(10,*) TLEN,TWID,PERIOD,DEPT,ZETA,GAMMA,MU
c      WRITE(*,*) TLEN,TWID,PERIOD,DEPT,ZETA,GAMMA,MU
C
C FIND FREQUENCIES(wavenumbers) AND DEPTH DEPENDENT PARAMETERS
C
      WX0= TWOPI/TLEN
      WY0= TWOPI/TWID
      DO 110 J=1,NUMB(4,2)
      WY=REAL(J-1)*WY0
      DO 100 I=1,NUMB(4,1)
      WX=REAL(I-1)*WX0
      WZ=SQRT(WX*WX+WY*WY)
      FQ(I,J,1)=WX
      FQ(I,J,2)=WY
      FQ(I,J,3)=WZ
      FQ(I,J,4)=WZ*DEPT 
      FQ(I,J,5)=COSH(FQ(I,J,4))
100   FQ(I,J,6)=TANH(FQ(I,J,4))
      DO 110 I=NUMB(7,1),NUMB(8,1)
      WX=REAL(I-NUMB(9,1))*WX0
      WZ=SQRT(WX*WX+WY*WY)
      IM=I-NUMB(10,1)
      FQ(IM,J,1)=WX
      FQ(IM,J,2)=WY
      FQ(IM,J,3)=WZ
      FQ(IM,J,4)=WZ*DEPT 
      FQ(IM,J,5)=COSH(FQ(IM,J,4))
110   FQ(IM,J,6)=TANH(FQ(IM,J,4))
C
C input free-surface elevation and potential
C
      do 120 j=1,numb(8,2)
      do 120 i=1,numb(8,1)
        read(10,*) z1,z2 
        sf(i,j)=cmplx(z1,0.D0)
        vp(i,j)=cmplx(z2,0.D0) 
120   continue
      close(10)
      DELT=PERIOD/IPER 
      
c      DO 108 J=1,NUMB(8,2)
c108   WRITE(*,*) (REAL(SF(I,J)),I=1,2)
C
C initialize fft
C
      call fft(ndx,ndy,1,1,sf,vx,vy,wrkx,wrky)
C
C transfrom from physical to wavenumber space
C
      call fft(ndx,ndy,0,1,sf,vx,vy,wrkx,wrky)
      call fft(ndx,ndy,0,1,vp,vx,vy,wrkx,wrky)
C
c re-define MEAN free surface to be zero
C
      sf(1,1)=0.
C
C de-alias the input to remove energy from modes higher than [2**p,2**q]
C
      call alia(ndx,ndy,sf)
      call alia(ndx,ndy,vp)
C
C smooth the input to remove energy near (and less than)  [2**p,2**q]
C
      call smot1(ndx,ndy,sf,vp,fq,0.9)
C
C CHECK TIME STEP
C
      TMAX=SQRT(8./FQ(NUMB(4,1),NUMB(4,2),3))
      IF(DELT.GT.TMAX)           STOP 'bad DELT'
      WRITE(11,*) NUMB(1,1),NUMB(2,1),NUMB(2,2),ISTP,
     1              DELT,TLEN,TWID,PERIOD,DEPT,ZETA,GAMMA,MU
     
      DO 105 J=1,NUMB(4,2)
105   WRITE(13,*) TIME,(SF(I,J),I=1,NUMB(4,1)),
     2               (SF(I,J),I=NUMB(7,1),NUMB(8,1))
     
      DO 106 J=1,NUMB(8,2)
      DO 106 I=1,NUMB(8,1)
      SF_temp(i,j) = SF(i,j)
      VP_temp(i,j) = VP(i,j)
106   CONTINUE

      CALL FFT(ndx,ndy,0,-1,SF_temp,vx,vy,wrkx,wrky)
      DO 107 J=1,NUMB(4,2)
107   WRITE(14,*) TIME,(SF_temp(I,J),I=1,NUMB(4,1)),
     3               (SF_temp(I,J),I=NUMB(7,1),NUMB(8,1))
C
C INPUT BOTTOM CONFIGURATION 
C
c      DO 140 J=1,NUMB(8,2)
c      DO 140 I=1,NUMB(8,1)
c140   READ(9,*) BAMP(I,J)      
c      CLOSE(9)
      do 140 j=1,numb(8,2)
      do 140 i=1,numb(8,1)
        read(9,*) z1b,z2b 
        BAMP(i,j)=cmplx(z1b,0.D0)
        VPB(i,j)=cmplx(z2b,0.D0) 
140   continue
      close(9)
C
C initialize fft
C
      call fft(ndx,ndy,1,1,BAMP,vx,vy,wrkx,wrky)
C
C transfrom from physical to wavenumber space
C
      call fft(ndx,ndy,0,1,BAMP,vx,vy,wrkx,wrky)
      call fft(ndx,ndy,0,1,VPB,vx,vy,wrkx,wrky)
C
c re-define MEAN free surface to be zero
C
      BAMP(1,1)=0.
C
C de-alias the input to remove energy from modes higher than [2**p,2**q]
C
      call alia(ndx,ndy,BAMP)
      call alia(ndx,ndy,VPB)
      
      WRITE(*,*) (REAL(BAMP(I,2)),I=1,3)
C
C smooth the input to remove energy near (and less than)  [2**p,2**q]
C
      call smot1(ndx,ndy,BAMP,VPB,fq,0.9)
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE UPFR(NDX,NDY,K,SURF,VPOT,R1,R2)
C
C MOVE TO THE NEXT TIME STEP
C
      IMPLICIT  REAL(A-H,M,O-Z)
      INTEGER   NUMB(20,2)
      REAL      MU
      COMPLEX   SURF(NDX/2+1,NDY/4+1,4),VPOT(NDX/2+1,NDY/4+1,4)
      COMPLEX   R1(NDX/2+1,NDY/4+1),R2(NDX/2+1,NDY/4+1)
      COMMON/ABLOCK/ISTP,DELT,TLEN,TWID,PERIOD,ENERGY0,DEPT,
     1              ZETA,GAMMA,MU
      COMMON/IBLOCK/NUMB
C
      DO 200 J=1,NUMB(4,2)
      DO 100 I=1,NUMB(4,1)
      SURF(I,J,K)=DELT*R1(I,J)
100   VPOT(I,J,K)=DELT*R2(I,J)
      DO 200 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      SURF(IM,J,K)=DELT*R1(IM,J)
200   VPOT(IM,J,K)=DELT*R2(IM,J)
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE SUMA(NDX,NDY,K,SURF,VPOT,SAVE,SF,VP)
C
C MOVE TO THE NEXT TIME STEP
C
      IMPLICIT  REAL(A-H,M,O-Z)
      INTEGER   NUMB(20,2)
      COMPLEX   SURF(NDX/2+1,NDY/4+1,4),VPOT(NDX/2+1,NDY/4+1,4)
      COMPLEX   SAVE(NDX/2+1,NDY/4+1,2),SF(NDX,NDY),VP(NDX,NDY)
      COMMON/IBLOCK/NUMB
C
C ROUTE THE PROGRAM
C
      GOTO(100,300,300,300,500) K
C
C SAVE THE RESULTS AT THE PRIOR TIME STEP (K=1)
C
100   CONTINUE
      DO 200 J=1,NUMB(4,2)
      DO 150 I=1,NUMB(4,1)
      SAVE(I,J,1)=SF(I,J)
150   SAVE(I,J,2)=VP(I,J)
      DO 200 I=NUMB(7,1),NUMB(8,1)
      SAVE(I-NUMB(10,1),J,1)=SF(I,J)
200   SAVE(I-NUMB(10,1),J,2)=VP(I,J)
      RETURN
C
C SAVE THE STEPS IN THE RUNGE-KUTTA METHOD (K=2,3, OR 4)
C
300   CONTINUE
      FAC=0.5D0
      IF(K.EQ.4) FAC=1.0D0
      DO 400 J=1,NUMB(4,2)
      DO 350 I=1,NUMB(4,1)
      SF(I,J)=SAVE(I,J,1)+FAC*SURF(I,J,K-1)
350   VP(I,J)=SAVE(I,J,2)+FAC*VPOT(I,J,K-1)
      DO 400 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      SF(I,J)=SAVE(IM,J,1)+FAC*SURF(IM,J,K-1)
400   VP(I,J)=SAVE(IM,J,2)+FAC*VPOT(IM,J,K-1)
      RETURN
C
C FIND THE SOLUTION AT THE NEXT TIME STEP (K=5)
C
500   CONTINUE
      DO 600 J=1,NUMB(4,2)
      DO 550 I=1,NUMB(4,1)
      SF(I,J)=SAVE(I,J,1)+(SURF(I,J,1)+2.*(SURF(I,J,2)+SURF(I,J,3))
     1                    +SURF(I,J,4))/6.
550   VP(I,J)=SAVE(I,J,2)+(VPOT(I,J,1)+2.*(VPOT(I,J,2)+VPOT(I,J,3))
     1                    +VPOT(I,J,4))/6.
      DO 600 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      SF(I,J)=SAVE(IM,J,1)+(SURF(IM,J,1)
     1   +2.*(SURF(IM,J,2)+SURF(IM,J,3))
     2       +SURF(IM,J,4))/6.
600   VP(I,J)=SAVE(IM,J,2)+(VPOT(IM,J,1)
     1   +2.*(VPOT(IM,J,2) +VPOT(IM,J,3))
     2       +VPOT(IM,J,4))/6.
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE FREX(NDX,NDY,NDP,SF,Y1,WRK1,VX,VY,WRKX,WRKY)
C
C RAISE THE FREE SURFACE BY POWERS
C
      IMPLICIT  REAL(A-H,M,O-Z)
      INTEGER   NUMB(20,2)
      REAL      Y1(NDX,NDY,NDP)
      COMPLEX   SF(NDX,NDY),WRK1(NDX,NDY)
      COMPLEX   VX(NDX),VY(NDX),WRKX((NDX/2)*5),WRKY((NDY/2)*5)
      COMMON/IBLOCK/NUMB 
C
C INITIALIZE ARRAYS
C
      DO 10 J=1,NUMB(6,2)
      DO 10 I=1,NUMB(8,1)
10    WRK1(I,J)=SF(I,J)
      CALL FFT(NDX,NDY,0,-1,WRK1,VX,VY,WRKX,WRKY)
      DO 20 J=1,NUMB(8,2)
      DO 20 I=1,NUMB(8,1)
20    Y1(I,J,1)=REAL(WRK1(I,J))
C
C RAISE FREE SURFACE BY INTEGER POWERS
C
      DO 50 K=2,NUMB(1,1)-1
      DO 30 J=1,NUMB(8,2)
      DO 30 I=1,NUMB(8,1)
30    WRK1(I,J)=CMPLX(Y1(I,J,K-1)*Y1(I,J,1),0.)
C
C ELIMINATE ALIASING
C
      CALL FFT(NDX,NDY,0,1,WRK1,VX,VY,WRKX,WRKY)
      CALL ALIA(NDX,NDY,WRK1)
      CALL FFT(NDX,NDY,0,-1,WRK1,VX,VY,WRKX,WRKY)
C
C STORE POWERS OF FREE SURFACE
C
      DO 40 J=1,NUMB(8,2)
      DO 40 I=1,NUMB(8,1)
40    Y1(I,J,K)=REAL(WRK1(I,J))
50    CONTINUE
C
      RETURN
      END
C
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE BOTM(NDX,NDY,NDP,BAMP,Y2,WRK1,VX,VY,WRKX,WRKY)
C
C RAISE THE BOTTOM SURFACE BY POWERS
C
      IMPLICIT  REAL(A-H,M,O-Z)
      INTEGER   NUMB(20,2)
      REAL      Y2(NDX,NDY,NDP),BAMP(NDX,NDY)
      COMPLEX   VX(NDX),VY(NDX),WRKX((NDX/2)*5),WRKY((NDY/2)*5)
      COMPLEX   WRK1(NDX,NDY)
      COMMON/IBLOCK/NUMB 
C
C INITIALIZE ARRAYS
C
      DO 10 J=1,NUMB(8,2)
      DO 10 I=1,NUMB(8,1)
10    Y2(I,J,1)=BAMP(I,J) 
C
C RAISE BOTTOM SURFACE BY INTEGER POWERS
C
      DO 50 K=2,NUMB(1,1)-1
      DO 30 J=1,NUMB(8,2)
      DO 30 I=1,NUMB(8,1)
30    WRK1(I,J)=CMPLX(Y2(I,J,K-1)*Y2(I,J,1),0.D0)
C
C ELIMINATE ALIASING
C
      CALL FFT(NDX,NDY,0,1,WRK1,VX,VY,WRKX,WRKY)
      CALL ALIA(NDX,NDY,WRK1)
      CALL FFT(NDX,NDY,0,-1,WRK1,VX,VY,WRKX,WRKY)
C
C STORE POWERS OF FREE SURFACE
C
      DO 40 J=1,NUMB(8,2)
      DO 40 I=1,NUMB(8,1)
40    Y2(I,J,K)=REAL(WRK1(I,J))
50    CONTINUE
C
      RETURN
      END
C
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
c      SUBROUTINE BVP(NDX,NDY,NDP,VP,P1,P2,P3,P4,Y1,Y2,FQ,
c     1               WRK1,WRK2,WRK3,WRK4,VX,VY,WRKX,WRKY)
      SUBROUTINE BVP(NDX,NDY,NDP,VP,VPB,P1,P2,P24,P3,P4,Y1,Y2,FQ,
     1               WRK1,WRK2,WRK3,WRK4,VX,VY,WRKX,WRKY)
C
C SOLVE LAPLACE'S EQUATION
C
      IMPLICIT  REAL(A-H,M,O-Z)
      INTEGER   NUMB(20,2)
      REAL      Y1(NDX,NDY,NDP),Y2(NDX,NDY,NDP)
      REAL      FQ(NDX/2+1,NDY/4+1,6),MU
      COMPLEX   P1(NDX/2+1,NDY/4+1,NDP),P2(NDX/2+1,NDY/4+1,NDP)
      COMPLEX   P3(NDX/2+1,NDY/4+1,NDP),P4(NDX/2+1,NDY/4+1,NDP)
      COMPLEX   P24(NDX/2+1,NDY/4+1,NDP)
      COMPLEX   VX(NDX),VY(NDX),WRKX((NDX/2)*5),WRKY((NDY/2)*5)
      COMPLEX   VP(NDX,NDY),WRK1(NDX,NDY),WRK2(NDX,NDY)
      COMPLEX   VPB(NDX,NDY)
      COMPLEX   WRK3(NDX,NDY),WRK4(NDX,NDY),ZERO
      COMMON/ABLOCK/ISTP,DELT,TLEN,TWID,PERIOD,ENERGY0,DEPT,
     1              ZETA,GAMMA,MU
      COMMON/IBLOCK/NUMB 
      DATA ZERO/(0.,0.)/,TWOPI/6.283185307179586476925287/
cC
cC STORE FIRST TERM IN PERTURBATON SERIES
cC
c      DO 10 J=1,NUMB(4,2)
c      DO 5  I=1,NUMB(4,1)
c      P1(I,J,1)=VP(I,J)
c      P2(I,J,1)=ZERO 
c      P3(I,J,1)=P1(I,J,1)*FQ(I,J,3)*FQ(I,J,6) 
c5     P4(I,J,1)=P1(I,J,1)/FQ(I,J,5)
c      DO 10 I=NUMB(7,1),NUMB(8,1)
c      IM=I-NUMB(10,1)
c      P1(IM,J,1)=VP(I,J)     
c      P2(IM,J,1)=ZERO 
c      P3(IM,J,1)=P1(IM,J,1)*FQ(IM,J,3)*FQ(IM,J,6) 
c10    P4(IM,J,1)=P1(IM,J,1)/FQ(IM,J,5)
C
C STORE FIRST TERM IN PERTURBATON SERIES
C
      DO 10 J=1,NUMB(4,2)
      DO 5  I=1,NUMB(4,1)
      P1(I,J,1)=VP(I,J)
      IF (I.EQ.1.AND.J.EQ.1) THEN
      P2(I,J,1)=ZERO
      ELSE
      P2(I,J,1)=(P1(I,J,1)/FQ(I,J,5)-VPB(I,J))/FQ(I,J,6)
      ENDIF
      P24(I,J,1)=VPB(I,J)
      P3(I,J,1)=P1(I,J,1)*FQ(I,J,3)*FQ(I,J,6)
     1                      +P2(I,J,1)*FQ(I,J,3)/FQ(I,J,5)
5     P4(I,J,1)=P2(I,J,1)*FQ(I,J,3)
      DO 10 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      P1(IM,J,1)=VP(I,J)
      IF (I.EQ.1.AND.J.EQ.1) THEN
      P2(I,J,1)=ZERO
      ELSE
      P2(IM,J,1)=(P1(IM,J,1)/FQ(IM,J,5)-VPB(I,J))/FQ(IM,J,6)
      ENDIF
      P24(IM,J,1)=VPB(I,J)
      P3(IM,J,1)=P1(IM,J,1)*FQ(IM,J,3)*FQ(IM,J,6)
     1                      +P2(IM,J,1)*FQ(IM,J,3)/FQ(IM,J,5)
10    P4(IM,J,1)=P2(IM,J,1)*FQ(IM,J,3)

c      WRITE(*,*) (REAL(P4(I,1,1)),I=1,3)
C
C SOLVE THE BVP BY PERTURBATION
C
      DO 330 K=2,NUMB(1,1)
C
C NON-ZERO DIRICHLET BOUNDARY CONDITION ON THE MEAN FREE SURFACE 
C ZERO MATRIX THAT STORES SUM OF PERTURBATION SERIES
C
      DO 20 J=1,NUMB(8,2)
      DO 20 I=1,NUMB(8,1)
20    WRK1(I,J)=ZERO
C
C SUM PERTURBATION SERIES
C
      DERI=1.
      DO 80 K1=1,K-1
      K2=K-K1
      DERI=DERI/REAL(K1)
C
C PLACE INTO WORK MATRIX
C
      CALL  ALIA(NDX,NDY,WRK2)
      IF((K1/2)*2.EQ.K1) THEN 
        SIGN=(-1.)**(K1/2)
        DO 40 J=1,NUMB(4,2)
        DO 30 I=1,NUMB(4,1)
30      WRK2(I,J)=SIGN*P1(I,J,K2)*CMPLX(0.,FQ(I,J,3))**K1 
        DO 40 I=NUMB(7,1),NUMB(8,1)
        IM=I-NUMB(10,1)
40      WRK2(I,J)=SIGN*P1(IM,J,K2)*CMPLX(0.,FQ(IM,J,3))**K1 
      ELSE 
        IF(K1.EQ.1) THEN 
          DO 46 J=1,NUMB(4,2)
          DO 45 I=1,NUMB(4,1)
45        WRK2(I,J)=P3(I,J,K2)
          DO 46 I=NUMB(7,1),NUMB(8,1)
          IM=I-NUMB(10,1)
46        WRK2(I,J)=P3(IM,J,K2)
        ELSE 
          SIGN=(-1.D0)**((K1-1)/2) 
          DO 60 J=1,NUMB(4,2)
          DO 50 I=1,NUMB(4,1)
50        WRK2(I,J)=SIGN*P3(I,J,K2)*CMPLX(0.,FQ(I,J,3))**(K1-1) 
          DO 60 I=NUMB(7,1),NUMB(8,1)
          IM=I-NUMB(10,1)
60        WRK2(I,J)=SIGN*P3(IM,J,K2)*CMPLX(0.,FQ(IM,J,3))**(K1-1) 
        ENDIF 
      ENDIF 
C
C TRANSFORM FROM WAVE NUMBER DOMAIN TO PHYSICAL DOMAIN
C
      CALL FFT(NDX,NDY,0,-1,WRK2,VX,VY,WRKX,WRKY)
C
C STORE SUM
C
      DO 70 J=1,NUMB(8,2)
      DO 70 I=1,NUMB(8,1)
70    WRK1(I,J)=WRK1(I,J)-Y1(I,J,K1)*DERI*REAL(WRK2(I,J))
80    CONTINUE
C
C TRANSFORM FROM PHYSICAL DOMAIN TO WAVENUMBER DOMAIN
C
      CALL FFT(NDX,NDY,0,1,WRK1,VX,VY,WRKX,WRKY)
C
C ELIMINATE ALIASING AND ASSEMBLE THE RESULTS
C
      DO 95 J=1,NUMB(4,2)
      DO 90 I=1,NUMB(4,1)
90    P1(I,J,K)=WRK1(I,J)
      DO 95 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
95    P1(IM,J,K)=WRK1(I,J)
cC
cC NON-ZERO NEUMANN BOUNDARY CONDITION ON THE BOTTOM 
cC ZERO MATRIX THAT STORES SUM OF PERTURBATION SERIES
cC
c      DO 110 J=1,NUMB(8,2)
c      DO 110 I=1,NUMB(8,1)
c      WRK1(I,J)=ZERO
c110   WRK3(I,J)=ZERO 
cC
cC SUM PERTURBATION SERIES
cC
c      DERI=1.
c      DO 180 K1=1,K-1
c      K2=K-K1
c      DERI=DERI/REAL(K1)
cC
cC PLACE INTO WORK MATRIX
cC
c      CALL  ALIA(NDX,NDY,WRK2)
c      CALL  ALIA(NDX,NDY,WRK4)
c      IF(K1.EQ.1) THEN 
c        DO 120 J=1,NUMB(4,2)
c        DO 115 I=1,NUMB(4,1)
c        WRK2(I,J)=P4(I,J,K2)*CMPLX(0.,FQ(I,J,1))
c115     WRK4(I,J)=P4(I,J,K2)*CMPLX(0.,FQ(I,J,2))
c        DO 120 I=NUMB(7,1),NUMB(8,1)
c        IM=I-NUMB(10,1)
c        WRK2(I,J)=P4(IM,J,K2)*CMPLX(0.,FQ(IM,J,1)) 
c120     WRK4(I,J)=P4(IM,J,K2)*CMPLX(0.,FQ(IM,J,2))
c      ELSE 
c        IF(((K1-1)/2)*2.EQ.(K1-1)) THEN 
c          SIGN=(-1.)**((K1-1)/2)
c          DO 140 J=1,NUMB(4,2)
c          DO 130 I=1,NUMB(4,1)
c          WRK2(I,J)=SIGN*P4(I,J,K2)*CMPLX(0.,FQ(I,J,1))*
c     1          CMPLX(0.D0,FQ(I,J,3))**(K1-1) 
c130       WRK4(I,J)=SIGN*P4(I,J,K2)*CMPLX(0.,FQ(I,J,2))*
c     1          CMPLX(0.D0,FQ(I,J,3))**(K1-1) 
c          DO 140 I=NUMB(7,1),NUMB(8,1)
c          IM=I-NUMB(10,1)
c          WRK2(I,J)=SIGN*P4(IM,J,K2)*CMPLX(0.,FQ(IM,J,1))*
c     1          CMPLX(0.,FQ(IM,J,3))**(K1-1) 
c140       WRK4(I,J)=SIGN*P4(IM,J,K2)*CMPLX(0.,FQ(IM,J,2))*
c     1          CMPLX(0.,FQ(IM,J,3))**(K1-1) 
c        ELSE 
c          DO 145 J=1,NUMB(4,2)
c          WRK2(1,J)=P2(1,J,K2)*CMPLX(0.,FQ(I,J,1)) 
c145       WRK4(1,J)=P2(1,J,K2)*CMPLX(0.,FQ(I,J,2)) 
c          SIGN=(-1.)**((K1-2)/2) 
c          DO 160 J=1,NUMB(4,2)
c          DO 150 I=2,NUMB(4,1)
c          WRK2(I,J)=SIGN*P2(I,J,K2)*CMPLX(0.,FQ(I,J,1))*
c     1          CMPLX(0.,FQ(I,J,3))**(K1-2) 
c150       WRK4(I,J)=SIGN*P2(I,J,K2)*CMPLX(0.,FQ(I,J,2))*
c     1          CMPLX(0.,FQ(I,J,3))**(K1-2) 
c          DO 160 I=NUMB(7,1),NUMB(8,1)
c          IM=I-NUMB(10,1)
c          WRK2(I,J)=SIGN*P2(IM,J,K2)*CMPLX(0.,FQ(IM,J,1))*
c     1          CMPLX(0.,FQ(IM,J,3))**(K1-2) 
c160       WRK4(I,J)=SIGN*P2(IM,J,K2)*CMPLX(0.,FQ(IM,J,2))*
c     1          CMPLX(0.,FQ(IM,J,3))**(K1-2) 
c        ENDIF 
c      ENDIF 
cC
cC TRANSFORM FROM WAVE NUMBER DOMAIN TO PHYSICAL DOMAIN
cC
c      CALL FFT(NDX,NDY,0,-1,WRK2,VX,VY,WRKX,WRKY)
c      CALL FFT(NDX,NDY,0,-1,WRK4,VX,VY,WRKX,WRKY)
cC
cC STORE SUM
cC
c      DO 170 J=1,NUMB(8,2)
c      DO 170 I=1,NUMB(8,1)
c      WRK1(I,J)=WRK1(I,J)+Y2(I,J,K1)*DERI*REAL(WRK2(I,J))
c170   WRK3(I,J)=WRK3(I,J)+Y2(I,J,K1)*DERI*REAL(WRK4(I,J))
c180   CONTINUE
cC
cC TRANSFORM FROM PHYSICAL DOMAIN TO WAVENUMBER DOMAIN
cC
c      CALL FFT(NDX,NDY,0,1,WRK1,VX,VY,WRKX,WRKY)
c      CALL FFT(NDX,NDY,0,1,WRK3,VX,VY,WRKX,WRKY)
cC
cC ELIMINATE ALIASING AND ASSEMBLE THE RESULTS
cC
c      DO 195 J=1,NUMB(4,2)
c      DO 190 I=1,NUMB(4,1)
c190   P2(I,J,K)=WRK1(I,J)*CMPLX(0.,FQ(I,J,1))
c     1    +WRK3(I,J)*CMPLX(0.,FQ(I,J,2)) 
c      DO 195 I=NUMB(7,1),NUMB(8,1)
c      IM=I-NUMB(10,1)
c195   P2(IM,J,K)=WRK1(I,J)*CMPLX(0.,FQ(IM,J,1))
c     1    +WRK3(I,J)*CMPLX(0.,FQ(IM,J,2))
C
C NON-ZERO DIRICHLET BOUNDARY CONDITION ON THE BOTTOM 
C ZERO MATRIX THAT STORES SUM OF PERTURBATION SERIES
C
      DO 120 J=1,NUMB(8,2)
      DO 120 I=1,NUMB(8,1)
120    WRK1(I,J)=ZERO
C
C SUM PERTURBATION SERIES
C
      DERI=1.
      DO 180 K1=1,K-1
      K2=K-K1
      DERI=DERI/REAL(K1)
C
C PLACE INTO WORK MATRIX
C
      CALL  ALIA(NDX,NDY,WRK2)
      IF((K1/2)*2.EQ.K1) THEN 
        SIGN=(-1.)**(K1/2)
        DO 140 J=1,NUMB(4,2)
        DO 130 I=1,NUMB(4,1)
130      WRK2(I,J)=SIGN*P24(I,J,K2)*CMPLX(0.,FQ(I,J,3))**K1 
        DO 140 I=NUMB(7,1),NUMB(8,1)
        IM=I-NUMB(10,1)
140      WRK2(I,J)=SIGN*P24(IM,J,K2)*CMPLX(0.,FQ(IM,J,3))**K1 
      ELSE 
        IF(K1.EQ.1) THEN 
          DO 146 J=1,NUMB(4,2)
          DO 145 I=1,NUMB(4,1)
145        WRK2(I,J)=P4(I,J,K2)
          DO 146 I=NUMB(7,1),NUMB(8,1)
          IM=I-NUMB(10,1)
146        WRK2(I,J)=P4(IM,J,K2)
        ELSE 
          SIGN=(-1.D0)**((K1-1)/2) 
          DO 160 J=1,NUMB(4,2)
          DO 150 I=1,NUMB(4,1)
150        WRK2(I,J)=SIGN*P4(I,J,K2)*CMPLX(0.,FQ(I,J,3))**(K1-1) 
          DO 160 I=NUMB(7,1),NUMB(8,1)
          IM=I-NUMB(10,1)
160        WRK2(I,J)=SIGN*P4(IM,J,K2)*CMPLX(0.,FQ(IM,J,3))**(K1-1) 
        ENDIF 
      ENDIF 
C
C TRANSFORM FROM WAVE NUMBER DOMAIN TO PHYSICAL DOMAIN
C
      CALL FFT(NDX,NDY,0,-1,WRK2,VX,VY,WRKX,WRKY)
C
C STORE SUM
C
      DO 170 J=1,NUMB(8,2)
      DO 170 I=1,NUMB(8,1)
170    WRK1(I,J)=WRK1(I,J)-Y2(I,J,K1)*DERI*REAL(WRK2(I,J))
180    CONTINUE
C
C TRANSFORM FROM PHYSICAL DOMAIN TO WAVENUMBER DOMAIN
C
      CALL FFT(NDX,NDY,0,1,WRK1,VX,VY,WRKX,WRKY)
cC
cC ELIMINATE ALIASING AND ASSEMBLE THE RESULTS
cC
c      DO 95 J=1,NUMB(4,2)
c      DO 90 I=1,NUMB(4,1)
c90    P1(I,J,K)=WRK1(I,J)
c      DO 95 I=NUMB(7,1),NUMB(8,1)
c      IM=I-NUMB(10,1)
c95    P1(IM,J,K)=WRK1(I,J)
C
C ELIMINATE ALIASING AND ASSEMBLE THE RESULTS
C
      DO 195 J=1,NUMB(4,2)
      DO 190 I=1,NUMB(4,1)
      IF (I.EQ.1.AND.J.EQ.1) THEN
      P2(I,J,K)=ZERO
      ELSE
      P2(I,J,K)=(-WRK1(I,J)+P1(I,J,K)/FQ(I,J,5))/FQ(I,J,6)
      ENDIF
190    P24(I,J,K)=WRK1(I,J)
      DO 195 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      IF (I.EQ.1.AND.J.EQ.1) THEN
      P2(IM,J,K)=ZERO
      ELSE
      P2(IM,J,K)=(-WRK1(I,J)+P1(I,J,K)/FQ(I,J,5))/FQ(I,J,6)
      ENDIF
195    P24(IM,J,K)=WRK1(I,J)
cC
cC FIND P3(I,J,K) AND P4(I,J,K)
cC
c      DO 220 J=1,NUMB(4,2)
c      DO 230 I=1,NUMB(4,1)
c230   P3(I,J,K)=P1(I,J,K)*FQ(I,J,3)*FQ(I,J,6)+P2(I,J,K)/FQ(I,J,5)
c      DO 220 I=NUMB(7,1),NUMB(8,1)
c      IM=I-NUMB(10,1)
c220   P3(IM,J,K)=P1(IM,J,K)*FQ(IM,J,3)*FQ(IM,J,6)
c     1    +P2(IM,J,K)/FQ(IM,J,5)
c      P4(1,1,K)=P1(1,1,K)-P2(1,1,K)*DEPT
c      DO 250 J=2,NUMB(4,2)
c250   P4(1,J,K)=P1(1,J,K)/FQ(1,J,5)-P2(1,J,K)*FQ(1,J,6)/FQ(1,J,3)
c      DO 260 J=1,NUMB(4,2)
c      DO 270 I=2,NUMB(4,1)
c270   P4(I,J,K)=P1(I,J,K)/FQ(I,J,5)-P2(I,J,K)*FQ(I,J,6)/FQ(I,J,3)
c      DO 260 I=NUMB(7,1),NUMB(8,1)
c      IM=I-NUMB(10,1)
c260   P4(IM,J,K)=P1(IM,J,K)/FQ(IM,J,5)
c     1    -P2(IM,J,K)*FQ(IM,J,6)/FQ(IM,J,3)
C
C FIND P3(I,J,K) AND P4(I,J,K)
C
      DO 220 J=1,NUMB(4,2)
      DO 230 I=1,NUMB(4,1)
230   P3(I,J,K)=P1(I,J,K)*FQ(I,J,3)*FQ(I,J,6)
     1                      +P2(I,J,K)*FQ(I,J,3)/FQ(I,J,5)
      DO 220 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
220   P3(IM,J,K)=P1(IM,J,K)*FQ(IM,J,3)*FQ(IM,J,6)
     1    +P2(IM,J,K)*FQ(IM,J,3)/FQ(IM,J,5)
            
      DO 260 J=1,NUMB(4,2)
      DO 270 I=1,NUMB(4,1)
270   P4(I,J,K)=P2(I,J,K)*FQ(I,J,3)
      DO 260 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
260   P4(IM,J,K)=P2(IM,J,K)*FQ(IM,J,3)
C
330   CONTINUE
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE VELO(NDX,NDY,NDP,P1,P3,Y1,W,FQ,WRK1,
     1                   VX,VY,WRKX,WRKY)
C
C FIND WATER-PARTICLE VELOCITIES ON THE FREE SURFACE 
C
      IMPLICIT  REAL(A-H,M,O-Z)
      INTEGER   NUMB(20,2)
      REAL      Y1(NDX,NDY,NDP),FQ(NDX/2+1,NDY/4+1,6),MU
      COMPLEX   P1(NDX/2+1,NDY/4+1,NDP),P3(NDX/2+1,NDY/4+1,NDP)
      COMPLEX   VX(NDX),VY(NDX),WRKX((NDX/2)*5),WRKY((NDY/2)*5)
      COMPLEX   W(NDX,NDY),WRK1(NDX,NDY),ZERO
      COMMON/ABLOCK/ISTP,DELT,TLEN,TWID,PERIOD,ENERGY0,DEPT,
     1              ZETA,GAMMA,MU
      COMMON/IBLOCK/NUMB
      DATA ZERO/(0.,0.)/,TWOPI/6.283185307179586476925287/
C
C ZERO MATRICES WHICH SUM PERTURBATION SERIES
C
      DO 10 J=1,NUMB(8,2)
      DO 10 I=1,NUMB(8,1)
10    W(I,J)=ZERO
C
C ITERATE OVER ALL TERMS IN PERTURBATION SERIES
C
      DO 170 K=1,NUMB(1,1)
C
C FIND LEADING TERMS IN TAYLOR SERIES 
C
      CALL  ALIA(NDX,NDY,WRK1)
      DO 120 J=1,NUMB(4,2)
      DO 110 I=1,NUMB(4,1)
110   WRK1(I,J)=P3(I,J,K)
      DO 120 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
120   WRK1(I,J)=P3(IM,J,K)
C
C TRANSFORM FORM WAVENUMBER TO PHYSICAL DOMAIN 
C
      CALL FFT(NDX,NDY,0,-1,WRK1,VX,VY,WRKX,WRKY)
C
C STORE THE SUM 
C
      DO 130 J=1,NUMB(8,2)
      DO 130 I=1,NUMB(8,1)
130   W(I,J)=W(I,J)+REAL(WRK1(I,J))
C      
C DO A TAYLOR EXPANSION ABOUT MEAN WATERLINE
C
      DERI=1.
      DO 70 K1=1,NUMB(1,1)-K
      KPOW=K1+1
      DERI=DERI/REAL(K1)
C
C PLACE INTO WORK MATRIX
C
      CALL ALIA(NDX,NDY,WRK1)
      IF((KPOW/2)*2.EQ.KPOW) THEN 
      SIGN=(-1.D0)**(KPOW/2) 
      DO 30 J=1,NUMB(4,2)
      DO 20 I=1,NUMB(4,1)
20    WRK1(I,J)=SIGN*P1(I,J,K)*CMPLX(0.,FQ(I,J,3))**KPOW
      DO 30 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
30    WRK1(I,J)=SIGN*P1(IM,J,K)*CMPLX(0.,FQ(IM,J,3))**KPOW
      ELSE 
      SIGN=(-1.D0)*(K1/2) 
      DO 50 J=1,NUMB(4,2)
      DO 40 I=1,NUMB(4,1)
40    WRK1(I,J)=SIGN*P3(I,J,K)*CMPLX(0.,FQ(I,J,3))**K1 
      DO 50 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
50    WRK1(I,J)=SIGN*P3(IM,J,K)*CMPLX(0.,FQ(IM,J,3))**K1 
      ENDIF 
C
C TRANSFORM FROM WAVE NUMBER DOMAIN TO PHYSICAL DOMAIN
C
      CALL FFT(NDX,NDY,0,-1,WRK1,VX,VY,WRKX,WRKY)
C
C STORE SUM
C
      DO 60 J=1,NUMB(8,2)
      DO 60 I=1,NUMB(8,1)
60    W(I,J)=W(I,J)+Y1(I,J,K1)*DERI*REAL(WRK1(I,J))
70    CONTINUE
170   CONTINUE
      
c      WRITE(*,*) (REAL(W(I,1)),I=1,3)
C
C ELIMINATE ALIASING
C
      CALL   FFT(NDX,NDY,0,1,W,VX,VY,WRKX,WRKY)
      CALL  ALIA(NDX,NDY,W)
      
c      WRITE(*,*) (REAL(P3(I,1,1)),I=1,3)
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE RIGH(NDX,NDY,NDP,R1,R2,RB1,RB2,SF,VP,W,BAMP,VPB,WB,FQ,
     1                WRK1,WRK2,WRK3,WRK4,VX,VY,WRKX,WRKY)
C
C ASSEMBLE THE RIGHT HAND SIDE
C
      IMPLICIT  REAL(A-H,M,O-Z)
      INTEGER   NUMB(20,2)
      REAL      FQ(NDX/2+1,NDY/4+1,6),MU
      COMPLEX   R1(NDX/2+1,NDY/4+1),R2(NDX/2+1,NDY/4+1)
      COMPLEX   RB1(NDX/2+1,NDY/4+1),RB2(NDX/2+1,NDY/4+1)
      COMPLEX   RB1_temp(NDX,NDY),RB2_temp(NDX,NDY)
      COMPLEX   SF(NDX,NDY),VP(NDX,NDY),W(NDX,NDY),ONE
      COMPLEX   BAMP(NDX,NDY),VPB(NDX,NDY),WB(NDX,NDY)
      COMPLEX   VX(NDX),VY(NDX),WRKX((NDX/2)*5),WRKY((NDY/2)*5)
      COMPLEX   WRK1(NDX,NDY),WRK2(NDX,NDY)
      COMPLEX   WRK3(NDX,NDY),WRK4(NDX,NDY)
      COMMON/ABLOCK/ISTP,DELT,TLEN,TWID,PERIOD,ENERGY0,DEPT,
     1              ZETA,GAMMA,MU
      COMMON/IBLOCK/NUMB
      DATA ONE/(0.,1.)/
C
C FIND THE FREE SURFACE AND POTENTIAL'S SLOPES
C
      CALL  ALIA(NDX,NDY,WRK1)
      CALL  ALIA(NDX,NDY,WRK2)
      CALL  ALIA(NDX,NDY,WRK3)
      CALL  ALIA(NDX,NDY,WRK4)
      DO 20 J=1,NUMB(4,2)
      DO 10 I=1,NUMB(4,1)
      WRK1(I,J)=SF(I,J)*FQ(I,J,1)*ONE
      WRK2(I,J)=SF(I,J)*FQ(I,J,2)*ONE
      WRK3(I,J)=VP(I,J)*FQ(I,J,1)*ONE
10    WRK4(I,J)=VP(I,J)*FQ(I,J,2)*ONE
      DO 20 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      WRK1(I,J)=SF(I,J)*FQ(IM,J,1)*ONE
      WRK2(I,J)=SF(I,J)*FQ(IM,J,2)*ONE
      WRK3(I,J)=VP(I,J)*FQ(IM,J,1)*ONE
20    WRK4(I,J)=VP(I,J)*FQ(IM,J,2)*ONE
C
C EVALUATE LINEAR TERMS
C
      DO 60 J=1,NUMB(4,2)
      DO 50 I=1,NUMB(4,1)
      R1(I,J)= W(I,J)
50    R2(I,J)=-SF(I,J)
      DO 60 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      R1(IM,J)= W(I,J)
60    R2(IM,J)=-SF(I,J)
C
C EVALUATE LINEAR TERMS - BOTTOM
C
c      DO 260 J=1,NUMB(4,2)
c      DO 250 I=1,NUMB(4,1)
c      RB1(I,J)= WB(I,J)
c250    RB2(I,J)=(1.0/GAMMA-1)*BAMP(I,J)+ZETA*RB1(I,J)*SQRT(DEPT)
c      DO 260 I=NUMB(7,1),NUMB(8,1)
c      IM=I-NUMB(10,1)
c      RB1(IM,J)= WB(I,J)
c260    RB2(IM,J)=(1.0/GAMMA-1)*BAMP(I,J)+ZETA*RB1(IM,J)*SQRT(DEPT)

c---------- BEGIN - Ninh 07/28/2016 ----------
      
C
C TRANSFORM FROM WAVE NUMBER DOMAIN TO PHYSICAL DOMAIN
C
      CALL FFT(NDX,NDY,0,-1,BAMP,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WB,VX,VY,WRKX,WRKY)
      
      DO 250 J=1,NUMB(8,2)
        DO 250 I=1,NUMB(8,1)
            RB1_temp(I,J)= WB(I,J)
            if (i.le.2*ndx/6.or.i.ge.4*ndx/6) then
                RB2_temp(I,J)=(1.0/GAMMA-1)*BAMP(I,J)
     1                          +ZETA*RB1_temp(I,J)*SQRT(DEPT)
            else
                RB2_temp(I,J)=(1.0/0.6-1)*BAMP(I,J)
     1                          +0.8*RB1_temp(I,J)*SQRT(DEPT)
            end if
250     continue
      
C
C TRANSFORM FROM PHYSICAL DOMAIN TO WAVE NUMBER DOMAIN
C
      CALL FFT(NDX,NDY,0,1,RB1_temp,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,RB2_temp,VX,VY,WRKX,WRKY)

      DO 260 J=1,NUMB(4,2)
        DO 255 I=1,NUMB(4,1)
            RB1(I,J)= RB1_temp(I,J)
            RB2(I,J)= RB2_temp(I,J)
255     continue
        DO 260  I=NUMB(7,1),NUMB(8,1)
            IM=I-NUMB(10,1)
            RB1(IM,J)= RB1_temp(I,J)
            RB2(IM,J)= RB2_temp(I,J)
260     continue

C
C TRANSFORM FROM PHYSICAL DOMAIN TO WAVE NUMBER DOMAIN
C
      CALL FFT(NDX,NDY,0,1,BAMP,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,WB,VX,VY,WRKX,WRKY)
      
c---------- END - Ninh 07/28/2016 ----------
            
c      WRITE(*,*) "GAMMA ",GAMMA
C      WRITE(*,*) (REAL(W(I,1)),I=1,3)

      IF(NUMB(1,1).EQ.1)RETURN
C
C TRANSFORM FROM WAVE NUMBER DOMAIN TO PHYSICAL DOMAIN
C
      CALL FFT(NDX,NDY,0,-1,WRK1,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WRK2,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WRK3,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WRK4,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,   W,VX,VY,WRKX,WRKY)
C
C CALCULATE NONLINEAR TERMS
C
      DO 90 J=1,NUMB(8,2)
      DO 90 I=1,NUMB(8,1)
      SX=REAL(WRK1(I,J))
      SY=REAL(WRK2(I,J))
      PX=REAL(WRK3(I,J))
      PY=REAL(WRK4(I,J))
      PZ=REAL(W(I,J))
      WRK1(I,J)= PX*SX+PY*SY
      WRK2(I,J)=(PX*PX+PY*PY)*0.5
      WRK3(I,J)= PZ*PZ
90    WRK4(I,J)= SX*SX+SY*SY
C
C ELIMINATE ALIASING
C
      CALL FFT(NDX,NDY,0,1,WRK1,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,WRK2,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,WRK3,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,WRK4,VX,VY,WRKX,WRKY)
      CALL  ALIA(NDX,NDY,WRK3)
      CALL  ALIA(NDX,NDY,WRK4)
      CALL FFT(NDX,NDY,0,-1,WRK3,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WRK4,VX,VY,WRKX,WRKY)
C
C CALCULATE NONLINEAR TERMS
C
      DO 100 J=1,NUMB(8,2)
      DO 100 I=1,NUMB(8,1)
      PZ1= REAL(W(I,J))
      PZ2= REAL(WRK3(I,J))
      SGSG=REAL(WRK4(I,J))
      WRK3(I,J)= PZ1*SGSG
100   WRK4(I,J)=(PZ2*SGSG+PZ2)*0.5
C
C TRANSFORM FROM PHYSICAL DOMAIN TO WAVENUMBER DOMAIN
C
      CALL FFT(NDX,NDY,0,1,WRK3,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,WRK4,VX,VY,WRKX,WRKY)
C
C COMPLETE EVALUATION OF KINEMATIC AND DYNAMIC BOUNDARY CONDITIONS
C
      DO 120 J=1,NUMB(4,2)
      DO 110 I=1,NUMB(4,1)
      R1(I,J)=R1(I,J)-WRK1(I,J)+WRK3(I,J)
110   R2(I,J)=R2(I,J)-WRK2(I,J)+WRK4(I,J)
      DO 120 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      R1(IM,J)=R1(IM,J)-WRK1(I,J)+WRK3(I,J)
120   R2(IM,J)=R2(IM,J)-WRK2(I,J)+WRK4(I,J)
C
c ---------BOTTOM---------
C
      DO 220 J=1,NUMB(4,2)
      DO 210 I=1,NUMB(4,1)
      WRK1(I,J)=BAMP(I,J)*FQ(I,J,1)*ONE
      WRK2(I,J)=BAMP(I,J)*FQ(I,J,2)*ONE
      WRK3(I,J)=VPB(I,J)*FQ(I,J,1)*ONE
210    WRK4(I,J)=VPB(I,J)*FQ(I,J,2)*ONE
      DO 220 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      WRK1(I,J)=BAMP(I,J)*FQ(IM,J,1)*ONE
      WRK2(I,J)=BAMP(I,J)*FQ(IM,J,2)*ONE
      WRK3(I,J)=VPB(I,J)*FQ(IM,J,1)*ONE
220    WRK4(I,J)=VPB(I,J)*FQ(IM,J,2)*ONE
C
C TRANSFORM FROM WAVE NUMBER DOMAIN TO PHYSICAL DOMAIN
C
      CALL FFT(NDX,NDY,0,-1,WRK1,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WRK2,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WRK3,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WRK4,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,   WB,VX,VY,WRKX,WRKY)
C
C CALCULATE NONLINEAR TERMS
C
      DO 290 J=1,NUMB(8,2)
      DO 290 I=1,NUMB(8,1)
      SBX=REAL(WRK1(I,J))
      SBY=REAL(WRK2(I,J))
      PBX=REAL(WRK3(I,J))
      PBY=REAL(WRK4(I,J))
      PBZ=REAL(WB(I,J))
      WRK1(I,J)= PBX*SBX+PBY*SBY
      WRK2(I,J)=(PBX*PBX+PBY*PBY)*0.5
      WRK3(I,J)= PBZ*PBZ
290    WRK4(I,J)= SBX*SBX+SBY*SBY
C
C ELIMINATE ALIASING
C
      CALL FFT(NDX,NDY,0,1,WRK1,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,WRK2,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,WRK3,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,WRK4,VX,VY,WRKX,WRKY)
      CALL  ALIA(NDX,NDY,WRK3)
      CALL  ALIA(NDX,NDY,WRK4)
      CALL FFT(NDX,NDY,0,-1,WRK3,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,-1,WRK4,VX,VY,WRKX,WRKY)
C
C CALCULATE NONLINEAR TERMS
C
      DO 300 J=1,NUMB(8,2)
      DO 300 I=1,NUMB(8,1)
      PBZ1= REAL(WB(I,J))
      PBZ2= REAL(WRK3(I,J))
      SBGSBG=REAL(WRK4(I,J))
      WRK3(I,J)= PBZ1*SBGSBG
300   WRK4(I,J)=(PBZ2*SBGSBG+PBZ2)*0.5
C
C TRANSFORM FROM PHYSICAL DOMAIN TO WAVENUMBER DOMAIN
C
      CALL FFT(NDX,NDY,0,1,WRK3,VX,VY,WRKX,WRKY)
      CALL FFT(NDX,NDY,0,1,WRK4,VX,VY,WRKX,WRKY)
C
C COMPLETE EVALUATION OF KINEMATIC AND DYNAMIC BOUNDARY CONDITIONS
C
      DO 320 J=1,NUMB(4,2)
      DO 310 I=1,NUMB(4,1)
      RB1(I,J)=RB1(I,J)-WRK1(I,J)+WRK3(I,J)
310   RB2(I,J)=RB2(I,J)-WRK2(I,J)+WRK4(I,J)+ZETA*RB1(I,J)*SQRT(DEPT)
     1                  -ZETA*WB(I,J)*SQRT(DEPT)
      DO 320 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      RB1(IM,J)=RB1(IM,J)-WRK1(I,J)+WRK3(I,J)
320   RB2(IM,J)=RB2(IM,J)-WRK2(I,J)+WRK4(I,J)+ZETA*RB1(IM,J)*SQRT(DEPT)
     1                  -ZETA*WB(I,J)*SQRT(DEPT)
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE ALIA(NDX,NDY,VECT)
C
C ELIMINATE ALIASING
C
      IMPLICIT REAL(A-H,M,O-Z)
      INTEGER  NUMB(20,2)
      COMPLEX  VECT(NDX,NDY),ZERO
      COMMON/IBLOCK/NUMB
      DATA ZERO/(0.,0.)/
C
      DO 10 J=        1,NUMB(6,2)
      DO 10 I=NUMB(5,1),NUMB(6,1) 
10    VECT(I,J)=ZERO
      DO 30 J=NUMB(5,2),NUMB(6,2)
      DO 20 I=        1,NUMB(4,1)
20    VECT(I,J)=ZERO
      DO 30 I=NUMB(7,1),NUMB(8,1)
30    VECT(I,J)=ZERO
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE SMOT1(NDX,NDY,SF,VP,FQ,FAC)
C
C SMOOTH THE DATA
C
      IMPLICIT REAL(A-H,M,O-Z)
      INTEGER  NUMB(20,2)
      REAL     FQ(NDX/2+1,NDY/4+1,6)
      COMPLEX  SF(NDX,NDY),VP(NDX,NDY)
      COMMON/IBLOCK/NUMB
      DATA ZERO/(0.,0.)/
C
      WXMAX=FAC*FQ(NUMB(4,1),NUMB(4,2),1)
      WYMAX=FAC*FQ(NUMB(4,1),NUMB(4,2),2)
C
      DO 20 J=1,NUMB(4,2)
      DO 10 I=1,NUMB(4,1)
      IF(FQ(I,J,1).GE.WXMAX.OR.
     1   FQ(I,J,2).GE.WYMAX)THEN
      VP(I,J)=ZERO
      SF(I,J)=ZERO
      ENDIF
10    CONTINUE
      DO 20 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      IF(ABS(FQ(IM,J,1)).GE.WXMAX.OR.ABS(FQ(IM,J,2)).GE.WYMAX)THEN
      VP(I,J)=ZERO
      SF(I,J)=ZERO
      ENDIF
20    CONTINUE
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      SUBROUTINE SMOT2(NDX,NDY,WK,FQ,FAC)
C
C SMOOTH THE DATA
C
      IMPLICIT REAL(A-H,M,O-Z)
      INTEGER  NUMB(20,2)
      REAL     FQ(NDX/2+1,NDY/4+1,6)
      COMPLEX  WK(NDX,NDY)
      COMMON/IBLOCK/NUMB
      DATA ZERO/(0.,0.)/
C
      WXMAX=FAC*FQ(NUMB(4,1),NUMB(4,2),1)
      WYMAX=FAC*FQ(NUMB(4,1),NUMB(4,2),2)
C
      DO 20 J=1,NUMB(4,2)
      DO 10 I=1,NUMB(4,1)
      IF(FQ(I,J,1).GE.WXMAX.OR.
     1   FQ(I,J,2).GE.WYMAX) WK(I,J)=ZERO
10    CONTINUE
      DO 20 I=NUMB(7,1),NUMB(8,1)
      IM=I-NUMB(10,1)
      IF(ABS(FQ(IM,J,1)).GE.WXMAX.OR.ABS(FQ(IM,J,2)).GE.WYMAX)THEN
      WK(I,J)=ZERO
      ENDIF
20    CONTINUE
C
      RETURN
      END
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      subroutine un_pack(ndx,ndy,ap,a)
c
c copies a packed array AP (such as r1,r2,vpot,surf,save) into
c a full (unpacked) array with 0's in the "expanded area" 
C
      implicit real(a-h,m,o-z)
      integer  numb(20,2)
      complex  a(ndx,ndy),ap(ndx/2+1,ndy/4+1)
      COMMON/IBLOCK/NUMB 
C
      call alia(ndx,ndy,a)	! "zero" a first
      do 200 j=1,numb(4,2)
        do 100 i=1,numb(4,1)
          a(i,j)=ap(i,j)
100     continue
        do 110 i=numb(7,1),numb(8,1)
          a(i,j)=ap(i-numb(10,1),j)
110     continue
200   continue
c
      return
      end
C****************************
C CRAY REMOVE BELOW
C****************************
      subroutine outmode(iout2,NDX,NDY,time,SF)
C
C output modal amplitudes in a compact form
C
      IMPLICIT REAL(A-H,M,O-Z)
      INTEGER  NUMB(20,2)
      character*1 amp,asign
      complex  sf(ndx,ndy)
      COMMON/IBLOCK/NUMB 
C
      n41=min0(numb(4,1),24)
      write(iout2,900) time,(i,i=0,n41-1)
      nQ=numb(4,2)-1
      do 10 j=-nQ,-1
        jj=-j+1
        write(iout2,910) j,asign(sf(1,jj)),amp(sf(1,jj)),
     $           (asign(sf(numb(9,1)-i,jj)),amp(sf(numb(9,1)-i,jj)),
     $                i=1,n41-1)
10    continue
      do 20 j=0,nQ
        jj=j+1
        write(iout2,910) j,(asign(sf(i,jj)),amp(sf(i,jj)),i=1,n41)
20    continue
      return
900   format(1x,'time=',f11.7/2x,'j\i',24i3)
910   format(1x,i3,1x,24(1x,2a1))
      end
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      character*1 function asign(c)
      character*1  plus,minus
      complex c
      data plus/' '/,minus/'-'/
      ac=cabs(c)
      if(ac.ge.1.D0) then        ! ac >= 1
        asign=plus
      else
        asign=minus
      endif
      return
      end
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      character*1 function amp(c)
      character*1  blank,star
      complex c
      data blank/' '/,star/'*'/
      data eps/1.e-7/
c amp='a' for |c|=xEa  where 1 <= x < 10
c amp=' ' for |c| <= xE-9
c amp='*' for |c| >= xE+9
      ac=cabs(c)
      if(ac.eq.0.) then
        amp=blank
      else                       ! ac .ne. 0
        alc=alog10(ac)
        if(ac.ge.1.) then        ! ac >= 1
          if(alc.ge.10.) then    ! ac >= 1**10
            amp=star
          else                   ! 1 <= ac < 1**10 
            amp=char(48+int(alc))
          endif
        else                     ! ac < 1
          if(alc.le.-9.) then   ! ac <= 1**(-10)
            amp=blank
          else
c eps for the case where x=1 "exactly"
            amp=char(49+iabs(int(alc+eps)))
          endif
        endif
      endif
      return
      end

cC****************************
cC BEGIN THE NEXT SUBROUTINE *
cC****************************
c      SUBROUTINE FFT(NDX,NDY,INIT,IFFT,A,X,Y,WX,WY)
cC
cC 2-D FFT (NDX.GE.NDY)
cC
c      COMPLEX A(NDX,NDY),ZERO
c      COMPLEX X(NDX),Y(NDX),WX((NDX/2)*5),WY((NDY/2)*5)
c      DATA    ZERO/(0.,0.)/
cC
cC INITIALIZE
cC
c      IF(INIT.NE.0)THEN
cCMIC$ PROCESS
c      CALL CFFT2(1,1,NDX,X,WX,Y)
c      CALL CFFT2(1,1,NDY,X,WY,Y)
cCMIC$ ENDPROCESS
c      ENDIF
cC
cC PHYSICAL SPACE TO WAVENUMBER SPACE
cC
c      IF(INIT.EQ.0.AND.IFFT.EQ.1)THEN
c      FACT=1./REAL(NDX*NDY)
cCMIC$ DOGLOBAL
c      DO 20 J=1,NDY
c      DO 10 I=1,NDX
c10    X(I)=A(I,J)
c      CALL CFFT2(0,-1,NDX,X,WX,Y)
c      DO 20 I=1,NDX
c20    A(I,J)=Y(I)
cC
cCMIC$ DOGLOBAL
c      DO 40 I=1,NDX
c      DO 30 J=1,NDY
c30    X(J)=A(I,J)
c      CALL CFFT2(0,-1,NDY,X,WY,Y)
c      DO 40 J=1,NDY
c40    A(I,J)=Y(J)*FACT
c      ENDIF
cC
cC WAVENUMBER SPACE TO PHYSICAL SPACE
cC
c      IF(INIT.EQ.0.AND.IFFT.EQ.-1)THEN
c      NDX2=NDX+2
c      NDY2=NDY+2
c      DO 50 J=2,NDY/2
c      A(     1,NDY2-J)=CONJG(A(1,J))
c      DO 50 I=2,NDX
c50    A(NDX2-I,NDY2-J)=CONJG(A(I,J))
cC
cCMIC$ DOGLOBAL
c      DO 70 J=1,(NDY/4)+1
c      DO 60 I=1,NDX
c60    X(I)=A(I,J)
c      CALL CFFT2(0,1,NDX,X,WX,Y)
c      DO 70 I=1,NDX
c70    A(I,J)=Y(I)
cC
cCMIC$ DOGLOBAL
c      DO 80 J=(NDY/4)*3+1,NDY
c      DO 75 I=1,NDX
c75    X(I)=A(I,J)
c      CALL CFFT2(0,1,NDX,X,WX,Y)
c      DO 80 I=1,NDX
c80    A(I,J)=Y(I)
cC
cCMIC$ DOGLOBAL
c      DO 100 I=1,NDX
c      DO  90 J=1,NDY
c90    X(J)=A(I,J)
c      CALL CFFT2(0,1,NDY,X,WY,Y)
c      DO 100 J=1,NDY
c100   A(I,J)=Y(J)
c      ENDIF
cC
c      RETURN
c      END

cC****************************
cC BEGIN THE NEXT SUBROUTINE *
cC****************************

c      SUBROUTINE FFT(NDX,NDY,INIT,IFFT,A,X,Y,WX,WY)
cC
cC 2-D FFT (NDX.GE.NDY)
cC  "MIMIC" VERSION CALLING FAST INSTEAD OF CFFT2
cC
c      COMPLEX A(NDX,NDY)
c      COMPLEX X(NDX),Y(NDX),WX((NDX/2)*5),WY((NDY/2)*5)
cC
c      IF(INIT.NE.0) RETURN
cC
cC PHYSICAL SPACE TO WAVENUMBER SPACE
cC
c      IF(IFFT.EQ.1)THEN
c        FACT=1./REAL(NDX*NDY)
c        DO 20 J=1,NDY
c          CALL FAST(A(1,J),NDX,-1.)
c20      CONTINUE
cC
c        DO 40 I=1,NDX
c          DO 30 J=1,NDY
c            X(J)=A(I,J)
c30        CONTINUE
c          CALL FAST(X,NDY,-1.)
c          DO 35 J=1,NDY
c            A(I,J)=X(J)*FACT
c35        CONTINUE
c40      CONTINUE
c      ENDIF
cC
cC WAVENUMBER SPACE TO PHYSICAL SPACE
cC
c      IF(IFFT.EQ.-1)THEN
c        NDX2=NDX+2
c        NDY2=NDY+2
c        DO 50 J=2,NDY/2
c          A(1,NDY2-J)=CONJG(A(1,J))
c          DO 55 I=2,NDX
c            A(NDX2-I,NDY2-J)=CONJG(A(I,J))
c55        CONTINUE
c50      CONTINUE
cC
c        DO 70 J=1,NDY
c          CALL FAST(A(1,J),NDX,1.)
c70      CONTINUE
cC
c        DO 100 I=1,NDX
c          DO 90 J=1,NDY
c            X(J)=A(I,J)
c90        CONTINUE
c          CALL FAST(X,NDY,1.)
c          DO 95 J=1,NDY
c            A(I,J)=X(J)
c95        CONTINUE
c100     CONTINUE
c      ENDIF
cC
c      RETURN
c      END
cC
cC****************************
cC BEGIN THE NEXT SUBROUTINE *
cC****************************
c      SUBROUTINE FAST(X,NDIM,ONE)
cC
cC IMPLEMENT THE FAST FOURIER TRANSFORM
cC
c      IMPLICIT REAL(A-H,M,O-Z),INTEGER(I-L,N)
c      COMPLEX  T,U,W,X(NDIM)
c      DATA PI/3.1415926535897932384/
c      NPOW = 8
c      N1=2**NPOW
c      N2=N1/2
c      N3=N1-1
c      J= 1
c      DO 300 I=1,N3
c      IF(I.GE.J) GO TO 100
c      T=   X(J)
c      X(J)=X(I)
c      X(I)=T
c100   K=N2
c200   IF(K.GE.J) GO TO 300
c      J=J-K
c      K=K/2
c      GO TO 200
c300   J=J+K
c      DO 500 L=1,NPOW
c      LE1=2**L
c      LE2=LE1/2
c      ANG=PI/LE2
c      U=CMPLX(1.,0.)
c      W=CMPLX(COS(ANG),ONE*SIN(ANG))
c      DO 500 J=1,LE2
c      DO 400 I=J,N1,LE1
c      IP=I+LE2
c      T=    X(IP)*U
c      X(IP)=X(I)-T
c400   X(I)= X(I)+T
c500   U=U*W
c      IF(ONE.EQ.1.)RETURN
c      SCL= 1./REAL(N1)
c      DO 600 I=1,N1
c600   X(I)=X(I)*SCL
c      RETURN
c      END

C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************

      SUBROUTINE FFT(NDX,NDY,INIT,IFFT,A,X,Y,WX,WY)
C
C 2-D FFT (NDX.GE.NDY)
C  "MIMIC" VERSION CALLING FAST INSTEAD OF CFFT2
C
      COMPLEX A(NDX,NDY)
      COMPLEX X(NDX),Y(NDX),WX((NDX/2)*5),WY((NDY/2)*5)
C
      IF(INIT.NE.0) RETURN
C
C PHYSICAL SPACE TO WAVENUMBER SPACE
C
      IF(IFFT.EQ.1)THEN
        FACT=1./REAL(NDX*NDY)
        DO 20 J=1,NDY
          CALL FAST(A(1,J),NDX,-1)
20      CONTINUE
C
        DO 40 I=1,NDX
          DO 30 J=1,NDY
            X(J)=A(I,J)
30        CONTINUE
          CALL FAST(X,NDY,-1)
          DO 35 J=1,NDY
            A(I,J)=X(J)*FACT
35        CONTINUE
40      CONTINUE
      ENDIF
C
C WAVENUMBER SPACE TO PHYSICAL SPACE
C
      IF(IFFT.EQ.-1)THEN
        NDX2=NDX+2
        NDY2=NDY+2
        DO 50 J=2,NDY/2
          A(1,NDY2-J)=CONJG(A(1,J))
          DO 55 I=2,NDX
            A(NDX2-I,NDY2-J)=CONJG(A(I,J))
55        CONTINUE
50      CONTINUE
C
        DO 70 J=1,NDY
          CALL FAST(A(1,J),NDX,1)
70      CONTINUE
C
        DO 100 I=1,NDX
          DO 90 J=1,NDY
            X(J)=A(I,J)
90        CONTINUE
          CALL FAST(X,NDY,1)
          DO 95 J=1,NDY
            A(I,J)=X(J)
95        CONTINUE
100     CONTINUE
      ENDIF
C
      RETURN
      END
C********************************
      SUBROUTINE FAST(X,N1,IONE)
C********************************
C
C IMPLEMENT THE FAST FOURIER TRANSFORM
C
      IMPLICIT REAL (A-H,M,O-Z),INTEGER(I-L,N)
      COMPLEX  T,U,W,X(N1)
      DATA PI/3.1415926535897932/
      N2=N1/2
      N3=N1-1
C
      J= 1
      DO 300 I=1,N3
      IF(I.GE.J) GO TO 100
      T=   X(J)
      X(J)=X(I)
      X(I)=T
100   K=N2
200   IF(K.GE.J) GO TO 300
      J=J-K
      K=K/2
      GO TO 200
300   J=J+K
C
      LE2=1
      ANG=PI
C LOOP:
C DO XXX L=1,NPOW
10      LE1=2*LE2
        U=CMPLX(1.,0.)
        W=CMPLX(COS(ANG),IONE*SIN(ANG))
        DO 500 J=1,LE2
          DO 400 I=J,N1,LE1
            IP=I+LE2
            T=    X(IP)*U
            X(IP)=X(I)-T
            X(I)= X(I)+T
400       CONTINUE
          U=U*W
500     CONTINUE
        LE2=LE1
        ANG=ANG/2
C XXX   CONTINUE
      IF(LE1.LT.N1) GOTO 10
99    RETURN
      END
C
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      subroutine fg(ndx,ndy,ndp,f,g)
C
C Get the filter function f and the damping function g
C      
      integer   i
      real      f(ndx), g(ndx)
      
      do 100 i = 1,ndx
        read(23,*) f(i),g(i)
100   continue

      return
      end
C
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      subroutine treat(ndx,ndy,ndp,sf,vp,bamp,vpb,f)
C
C Initially treat the initial condition to mimic a wave maker
C      
      implicit  real(a-h,m,o-z),integer(i-l,n)
      integer   numb(20,2)
      real      f(ndx), g(ndx)
      COMPLEX   SF(NDX,NDY),VP(NDX,NDY)
      COMPLEX   BAMP(NDX,NDY),VPB(NDX,NDY)
      COMPLEX   VX(NDX),VY(NDX),WRKX((NDX/2)*5),WRKY((NDY/2)*5)
      COMMON/ABLOCK/ISTP,DELT,TLEN,TWID,PERIOD,ENERGY0,DEPT,
     1              ZETA,GAMMA,MU
      COMMON/IBLOCK/NUMB
      
      do 10 j = 1,ndy
        do 10 i=1,ndx
c            READ(20,*) YPOS,POTL
            READ(20,*) YPOS
            SF(I,J)=CMPLX(YPOS,0.0)
c            VP(I,J)=CMPLX(POTL,0.0)
          VP(I,J)=CMPLX(0,0)
c          BAMP(I,J)=CMPLX(BPOS,0)
          BAMP(I,J)=CMPLX(0,0)
c          VPB(I,J)=CMPLX(POTB,0)
          VPB(I,J)=CMPLX(0,0)
10    continue
C
C TRANSFORM FROM PHYSICAL DOMAIN TO WAVE NUMBER DOMAIN
C
      call fft(ndx,ndy,0,1,sf,vx,vy,wrkx,wrky)
      call fft(ndx,ndy,0,1,vp,vx,vy,wrkx,wrky)
      call fft(ndx,ndy,0,1,bamp,vx,vy,wrkx,wrky)
      call fft(ndx,ndy,0,1,vpb,vx,vy,wrkx,wrky)
            
      return
      end
C
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      subroutine readInput(ndx,ndy,ndp,sf_input,vp_input,
     1                      bamp_input,vpb_input)
C
C Read the data that serves as input for each time step
C      
      implicit  real(a-h,m,o-z),integer(i-l,n)
      integer   numb(20,2)            
      COMPLEX   SF_input(istp,ndx,ndy),VP_input(istp,ndx,ndy)
      COMPLEX   BAMP_input(istp,ndx,ndy),VPB_input(istp,ndx,ndy)
      COMMON/ABLOCK/ISTP,DELT,TLEN,TWID,PERIOD,ENERGY0,DEPT,
     1              ZETA,GAMMA,MU
      COMMON/IBLOCK/NUMB
C
C Read from file 20
C
      do 10 k=1,istp
        do 10 j=1,ndy
            do 10 i = 1,ndx
c                READ(20,*) YPOS,POTL
                READ(20,*) YPOS
                SF_input(k,i,j)=CMPLX(YPOS,0.0)
c                VP_input(k,i,j)=CMPLX(POTL,0.0)
                VP_input(k,i,j)=CMPLX(0,0.0)
c                BAMP_input(k,i,j)=CMPLX(BPOS,0)
                BAMP_input(k,i,j)=CMPLX(0,0)
c                VPB_input(k,i,j)=CMPLX(POTB,0)
                VPB_input(k,i,j)=CMPLX(0,0)
10    continue

c      write(*,*) 'readInput ', SF_input(1,1), SF_input(2,1)
      return
      end
C
C****************************
C BEGIN THE NEXT SUBROUTINE *
C****************************
      subroutine treatment(ndx,ndy,ndp,sf,vp,bamp,vpb,sf_input,vp_input,
     1                        bamp_input,vpb_input,f,g,itim)
C
C Treat the output after each time step to form new initial condition
C      
      implicit  real(a-h,m,o-z),integer(i-l,n)
      integer   numb(20,2)
      real      f(ndx),g(ndx)
      COMPLEX   SF(NDX,NDY),VP(NDX,NDY),BAMP(NDX,NDY),VPB(NDX,NDY)
      COMPLEX   SF_input(istp,NDX,NDY),VP_input(istp,NDX,NDY)
      COMPLEX   BAMP_input(istp,NDX,NDY),VPB_input(istp,NDX,NDY)
      COMPLEX   VX(NDX),VY(NDX),WRKX((NDX/2)*5),WRKY((NDY/2)*5)
      COMMON/ABLOCK/ISTP,DELT,TLEN,TWID,PERIOD,ENERGY0,DEPT,
     1              ZETA,GAMMA,MU
      COMMON/IBLOCK/NUMB
C
C TRANSFORM FROM WAVE NUMBER DOMAIN TO PHYSICAL DOMAIN
C
      call fft(ndx,ndy,0,-1,sf,vx,vy,wrkx,wrky)
      call fft(ndx,ndy,0,-1,vp,vx,vy,wrkx,wrky)
cC
cC Form new boundary condition
cC
c      k = mod(itim,iper)
c      if (k.eq.0) k=iper
      
c      do 20 i=1,ndx
c          sf(i) = sf(i)*g(i) + sf_input(k,i)
c          vp(i) = vp(i)*g(i) + vp_input(k,i)
c20    continue
C
C Form new boundary condition
C      
      do 20 j=1,ndy
        do 20 i=1,ndx
          sf(i,j) = sf(i,j)*g(i) + sf_input(itim,i,j)
c          vp(i,j) = vp(i,j)*g(i) + vp_input(itim,i,j)
20    continue
C
C TRANSFORM FROM PHYSICAL DOMAIN TO WAVE NUMBER DOMAIN
C
      call fft(ndx,ndy,0,1,sf,vx,vy,wrkx,wrky)
      call fft(ndx,ndy,0,1,vp,vx,vy,wrkx,wrky)
C
C ELIMINATE ALIASING
C
      CALL ALIA(NDX,NDY,sf)
      CALL ALIA(NDX,NDY,vp)

      return
      end
C
