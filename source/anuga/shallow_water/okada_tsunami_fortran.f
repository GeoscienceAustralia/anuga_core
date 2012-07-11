      subroutine fault_disp(alp,elon,elat,edep,strk,dip,length,wdt,
     &     disl1,disl2,rlon,rlat,dstmx,edsp,ndsp,zdsp,m,n)
c***********************************************************************
c   
c***********************************************************************
      Implicit REAL*8 (A-H,L,O-Z)
      integer m,n

      real*8,intent(in):: alp,elat(n),elon(n),edep(n),strk(n),dip(n)
      real*8,intent(in):: dstmx,length(n),wdt(n),disl1(n),disl2(n)
      real*8,intent(in):: rlon(m),rlat(m)
      real*8,intent(out):: zdsp(m),edsp(m),ndsp(m)
      
      logical:: DEBUG
c     Convenience flag for debugging -- will write various outputs to file
      DEBUG=.FALSE. 
   
c      print*, 'IN FORTRAN'

      IF(DEBUG.eqv..TRUE.) THEN  

          open(32, file='mylog.log')
          write(32,*) ' MYLOGGINGZ'

          write(32,*) 'WRITING TO FILE'
          write(32,*) 'Input pars:'
          write(32,*) 'alp=', alp
          write(32,*)  'elon=', elon
          write(32,*)  'elat=',elat
          write(32,*)  'edep=',edep
          write(32,*)  'strk=',strk
          write(32,*)  'dip=', dip
          write(32,*)  'length=', length
          write(32,*)  'wdt=',wdt
          write(32,*)  'disl1, disl2, disl3=', disl1, disl2, disl3
          write(32,*)  'rlon, rlat = ', rlon, rlat
          write(32,*)  'm , n = ', m, n
 
      END IF

      pi  = 4.*atan2(1.,1.)
      dtr = pi/180.d0 ! degrees to radians
      lt2k = 6371.d0*dtr ! latitude to kilometres

c Initialize displacements to zero
      do 50 j=1,m
         ndsp(j) = 0.
         edsp(j) = 0.
 50      zdsp(j) = 0.

c Loop over faults
      do 200 i=1,n
c Translate the origin from fault center to Okada's lower left corner
         tmp2 = cos(dip(i)*dtr)
         tmp1 = 0.5*sqrt(length(i)*length(i) + tmp2*wdt(i)*tmp2*wdt(i)) ! Half 'diagonal length' of surface projection of the slip area
         tmp2 = atan2(tmp2*wdt(i),length(i)) ! Angle of triangle connecting corners of surface projection of slip area
         odep = edep(i) + 0.5*wdt(i)*sin(dip(i)*dtr) ! Origin depth in Okada's reference frame

c       Compute origin for Okada's reference frame
c       GD CHANGES HERE
         oy = elat(i) - 1000.*tmp1*cos(tmp2-strk(i)*dtr)
         ox = elon(i) + 1000.*tmp1*sin(tmp2-strk(i)*dtr)
c         ox = elon(i) - 1000.*tmp1*cos(tmp2-strk(i)*dtr)
c         oy = elat(i) - 1000.*tmp1*sin(tmp2-strk(i)*dtr)

         IF (DEBUG.eqv..TRUE.) THEN
             write(32,*) 'Okada origin of subfault ',i,': ',ox,oy,odep
             write(32,*) 'Subfault centroid ',elon(i),elat(i)
             write(32,*) 'Strike, dip',strk(i), dip(i)
             write(32,*) 'More pars:', tmp1,tmp2, n, m, length(i),wdt(i)
             write(32,*) 'END'
         END IF
c Loop over displacement points
         do 100 j=1,m
c          Translate to new origin
           x = (rlon(j) - ox)*.001
           y = (rlat(j) - oy)*.001
           d = sqrt(x**2+y**2) ! Find distance from origin
           az = 90.-atan2(y,x)/dtr

           IF(DEBUG.eqv..TRUE.) THEN
               write(32,*) x,y,rlon(j),rlat(j),d,az,atan2(y,x)/dtr
               write(32,*) 'END'
           END IF

c          Rotate into Okada's reference frame, which has the strike
c          direction = positive x axis, with the origin at the deep end
c          of the slip with the most negative value along the x axis.
           y = d*sin(dtr*(strk(i)-az))
           x = d*cos(dtr*(strk(i)-az))
C           x = -d*sin(dtr*(strk(i)-az))
c           y = d*cos(dtr*(strk(i)-az))
           IF(DEBUG.eqv..TRUE.) THEN
               write(32,*) x,y
               write(32,*) 'END'
           END IF
c            write(*,*) x,y,d,az
c Skip this contribution of distance exceeds the threshhold            
            if (dstmx.gt.0. .and. sqrt(x*x+y*y).gt.dstmx) goto 100
c      SUBROUTINE  SRECTF(ALP,X,Y,DEP,AL1,AL2,AW1,AW2,                   01840000
c     *                   SD,CD,DISL1,DISL2,DISL3,                       01850000
c     *                   U1,U2,U3,U11,U12,U21,U22,U31,U32)              01860000
            IF(DEBUG.eqv..TRUE.) THEN
                write(32,*) 'SRECTF INPUTS: ', alp, x, y, odep, 0.d0
                write(32,*) length(i), 0.d0, wdt(i), sin(dip(i)*dtr)
                write(32,*) cos(dip(i)*dtr), disl1(i), disl2(i), 0.d0
            END IF

            call srectf(alp,x,y,odep,0.d0,length(i),0.d0,wdt(i),
     &           sin(dip(i)*dtr),cos(dip(i)*dtr),disl1(i),
     &           disl2(i),0.d0,u1,u2,u3,u11,u12,u21,u22,u31,u32)
            IF(DEBUG.eqv..TRUE.) THEN
                write(32,*) 'U3 is: ', u3
            END IF
            edsp(j) = edsp(j) - u2*cos(strk(i)*dtr)+u1*sin(strk(i)*dtr)
            ndsp(j) = ndsp(j) + u2*sin(strk(i)*dtr)+u1*cos(strk(i)*dtr)
            zdsp(j) = zdsp(j) + u3
 100     continue
 200  continue

      IF(DEBUG.eqv..TRUE.) THEN
          close(32)
      END IF
      RETURN
      END

C********************************************************************   00570000
      SUBROUTINE  SPOINT(ALP,X,Y,DEP,SD,CD,POT1,POT2,POT3,              00530000
     *                   U1,U2,U3,U11,U12,U21,U22,U31,U32)              00540000
      IMPLICIT REAL*8 (A-H,O-Z)                                         00550000
C                                                                       00560000
C********************************************************************   00570000
C*****                                                          *****   00580000
C*****    SURFACE DISPLACEMENT,STRAIN,TILT                      *****   00590000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   00600000
C*****                         CODED BY  Y.OKADA ... JAN 1985   *****   00610000
C*****                                                          *****   00620000
C********************************************************************   00630000
C                                                                       00640000
C***** INPUT                                                            00650000
C*****   ALP   : MEDIUM CONSTANT  MYU/(LAMDA+MYU)=1./((VP/VS)**2-1)     00660000
C*****   X,Y   : COORDINATE OF STATION                                  00670000
C*****   DEP     : SOURCE DEPTH                                         00680000
C*****   SD,CD : SIN,COS OF DIP-ANGLE                                   00690000
C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)00700000
C*****   POT1,POT2,POT3 : STRIKE-, DIP- AND TENSILE-POTENCY             00710000
C*****       POTENCY=(  MOMENT OF DOUBLE-COUPLE  )/MYU    FOR POT1,2    00720000
C*****       POTENCY=(INTENSITY OF ISOTROPIC PART)/LAMDA  FOR POT3      00730000
C                                                                       00740000
C***** OUTPUT                                                           00750000
C*****   U1, U2, U3      : DISPLACEMENT ( UNIT=(UNIT OF POTENCY) /      00760000
C*****                   :                     (UNIT OF X,Y,D)**2  )    00770000
C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF POTENCY) /      00780000
C*****   U31,U32         : TILT                (UNIT OF X,Y,D)**3  )    00790000
C                                                                       00800000
      DATA  F0,F1,F2,F3,F4,F5,F8,F9                                     00810000
     *      /0.D0, 1.D0, 2.D0, 3.D0, 4.D0, 5.D0, 8.D0, 9.D0/            00820000
      PI2=6.283185307179586D0                                           00830000
C-----                                                                  00840000
      D =DEP                                                            00850000
      P =Y*CD + D*SD                                                    00860000
      Q =Y*SD - D*CD                                                    00870000
      S =P*SD + Q*CD                                                    00880000
      X2=X*X                                                            00890000
      Y2=Y*Y                                                            00900000
      XY=X*Y                                                            00910000
      D2=D*D                                                            00920000
      R2=X2 + Y2 + D2                                                   00930000
      R =SQRT(R2)                                                       00940000
      R3=R *R2                                                          00950000
      R5=R3*R2                                                          00960000
      QR=F3*Q/R5                                                        00970000
      XR =F5*X2/R2                                                      00980000
      YR =F5*Y2/R2                                                      00990000
      XYR=F5*XY/R2                                                      01000000
      DR =F5*D /R2                                                      01010000
      RD =R + D                                                         01020000
      R12=F1/(R*RD*RD)                                                  01030000
      R32=R12*(F2*R + D)/ R2                                            01040000
      R33=R12*(F3*R + D)/(R2*RD)                                        01050000
      R53=R12*(F8*R2 + F9*R*D + F3*D2)/(R2*R2*RD)                       01060000
      R54=R12*(F5*R2 + F4*R*D +    D2)/R3*R12                           01070000
C-----                                                                  01080000
      A1= ALP*Y*(R12-X2*R33)                                            01090000
      A2= ALP*X*(R12-Y2*R33)                                            01100000
      A3= ALP*X/R3 - A2                                                 01110000
      A4=-ALP*XY*R32                                                    01120000
      A5= ALP*( F1/(R*RD) - X2*R32 )                                    01130000
      B1= ALP*(-F3*XY*R33      + F3*X2*XY*R54)                          01140000
      B2= ALP*( F1/R3 - F3*R12 + F3*X2*Y2*R54)                          01150000
      B3= ALP*( F1/R3 - F3*X2/R5) - B2                                  01160000
      B4=-ALP*F3*XY/R5 - B1                                             01170000
      C1=-ALP*Y*(R32 - X2*R53)                                          01180000
      C2=-ALP*X*(R32 - Y2*R53)                                          01190000
      C3=-ALP*F3*X*D/R5 - C2                                            01200000
C-----                                                                  01210000
      U1 =F0                                                            01220000
      U2 =F0                                                            01230000
      U3 =F0                                                            01240000
      U11=F0                                                            01250000
      U12=F0                                                            01260000
      U21=F0                                                            01270000
      U22=F0                                                            01280000
      U31=F0                                                            01290000
      U32=F0                                                            01300000
C======================================                                 01310000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 01320000
C======================================                                 01330000
      IF(POT1.EQ.F0)   GO TO 200                                        01340000
      UN=POT1/PI2                                                       01350000
      QRX=QR*X                                                          01360000
      FX=F3*X/R5*SD                                                     01370000
      U1 =U1 - UN*( QRX*X + A1*SD )                                     01380000
      U2 =U2 - UN*( QRX*Y + A2*SD )                                     01390000
      U3 =U3 - UN*( QRX*D + A4*SD )                                     01400000
      U11=U11- UN*( QRX* (F2-XR)        + B1*SD )                       01410000
      U12=U12- UN*(-QRX*XYR      + FX*X + B2*SD )                       01420000
      U21=U21- UN*( QR*Y*(F1-XR)        + B2*SD )                       01430000
      U22=U22- UN*( QRX *(F1-YR) + FX*Y + B4*SD )                       01440000
      U31=U31- UN*( QR*D*(F1-XR)        + C1*SD )                       01450000
      U32=U32- UN*(-QRX*DR*Y     + FX*D + C2*SD )                       01460000
C======================================                                 01470000
C=====    DIP-SLIP CONTRIBUTION   =====                                 01480000
C======================================                                 01490000
  200 IF(POT2.EQ.F0)   GO TO 300                                        01500000
      UN=POT2/PI2                                                       01510000
      SDCD=SD*CD                                                        01520000
      QRP=QR*P                                                          01530000
      FS=F3*S/R5                                                        01540000
      U1 =U1 - UN*( QRP*X - A3*SDCD )                                   01550000
      U2 =U2 - UN*( QRP*Y - A1*SDCD )                                   01560000
      U3 =U3 - UN*( QRP*D - A5*SDCD )                                   01570000
      U11=U11- UN*( QRP*(F1-XR)        - B3*SDCD )                      01580000
      U12=U12- UN*(-QRP*XYR     + FS*X - B1*SDCD )                      01590000
      U21=U21- UN*(-QRP*XYR            - B1*SDCD )                      01600000
      U22=U22- UN*( QRP*(F1-YR) + FS*Y - B2*SDCD )                      01610000
      U31=U31- UN*(-QRP*DR*X           - C3*SDCD )                      01620000
      U32=U32- UN*(-QRP*DR*Y    + FS*D - C1*SDCD )                      01630000
C========================================                               01640000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               01650000
C========================================                               01660000
  300 IF(POT3.EQ.F0)   GO TO 900                                        01670000
      UN=POT3/PI2                                                       01680000
      SDSD=SD*SD                                                        01690000
      QRQ=QR*Q                                                          01700000
      FQ=F2*QR*SD                                                       01710000
      U1 =U1 + UN*( QRQ*X - A3*SDSD )                                   01720000
      U2 =U2 + UN*( QRQ*Y - A1*SDSD )                                   01730000
      U3 =U3 + UN*( QRQ*D - A5*SDSD )                                   01740000
      U11=U11+ UN*( QRQ*(F1-XR)        - B3*SDSD )                      01750000
      U12=U12+ UN*(-QRQ*XYR     + FQ*X - B1*SDSD )                      01760000
      U21=U21+ UN*(-QRQ*XYR            - B1*SDSD )                      01770000
      U22=U22+ UN*( QRQ*(F1-YR) + FQ*Y - B2*SDSD )                      01780000
      U31=U31+ UN*(-QRQ*DR*X           - C3*SDSD )                      01790000
      U32=U32+ UN*(-QRQ*DR*Y    + FQ*D - C1*SDSD )                      01800000
C-----                                                                  01810000
  900 RETURN                                                            01820000
      END                                                               01830000
      SUBROUTINE  SRECTF(ALP,X,Y,DEP,AL1,AL2,AW1,AW2,                   01840000
     *                   SD,CD,DISL1,DISL2,DISL3,                       01850000
     *                   U1,U2,U3,U11,U12,U21,U22,U31,U32)              01860000
      IMPLICIT REAL*8 (A-H,O-Z)                                         01870000
C                                                                       01880000
C*********************************************************              01890000
C*****                                               *****              01900000
C*****    SURFACE DISPLACEMENT,STRAIN,TILT           *****              01910000
C*****    DUE TO RECTANGULAR FAULT IN A HALF-SPACE   *****              01920000
C*****              CODED BY  Y.OKADA ... JAN 1985   *****              01930000
C*****                                               *****              01940000
C*********************************************************              01950000
C                                                                       01960000
C***** INPUT                                                            01970000
C*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)=1./((VP/VS)**2-1)   01980000
C*****   X,Y     : COORDINATE OF STATION                                01990000
C*****   DEP     : SOURCE DEPTH                                         02000000
C*****   AL1,AL2 : FAULT LENGTH RANGE                                   02010000
C*****   AW1,AW2 : FAULT WIDTH RANGE                                    02020000
C*****   SD,CD   : SIN,COS OF DIP-ANGLE                                 02030000
C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)02040000
C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION      02050000
C                                                                       02060000
C***** OUTPUT                                                           02070000
C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL     )      02080000
C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /          02090000
C*****   U31,U32         : TILT                 UNIT OF X,Y,,,AW )      02100000
C                                                                       02110000
C***** SUBROUTINE USED...SRECTG                                         02120000
C                                                                       02130000
      DIMENSION  U(9),DU(9)                                             02140000
      DATA  F0, F1 / 0.D0, 1.D0 /                                       02150000
C-----                                                                  02160000
      P = Y*CD + DEP*SD                                                 02170000
      Q = Y*SD - DEP*CD                                                 02180000
      DO 1111  I=1,9                                                    02190000
 1111 U(I)=F0                                                           02200000
C-----                                                                  02210000
      DO 5555  K=1,2                                                    02220000
       IF(K.EQ.1)  ET=P-AW1                                             02230000
       IF(K.EQ.2)  ET=P-AW2                                             02240000
       DO 4444  J=1,2                                                   02250000
        IF(J.EQ.1)  XI=X-AL1                                            02260000
        IF(J.EQ.2)  XI=X-AL2                                            02270000
        JK=J+K                                                          02280000
        IF(JK.NE.3)  SIGN= F1                                           02290000
        IF(JK.EQ.3)  SIGN=-F1                                           02300000
        CALL SRECTG(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3,                02310000
     *           DU(1),DU(2),DU(3),DU(4),DU(5),DU(6),DU(7),DU(8),DU(9)) 02320000
        DO 3333  I=1,9                                                  02330000
         U(I)=U(I)+SIGN*DU(I)                                           02340000
 3333   CONTINUE                                                        02350000
 4444  CONTINUE                                                         02360000
 5555 CONTINUE                                                          02370000
      U1 =U(1)                                                          02380000
      U2 =U(2)                                                          02390000
      U3 =U(3)                                                          02400000
      U11=U(4)                                                          02410000
      U12=U(5)                                                          02420000
      U21=U(6)                                                          02430000
      U22=U(7)                                                          02440000
      U31=U(8)                                                          02450000
      U32=U(9)                                                          02460000
      RETURN                                                            02470000
      END                                                               02480000
      SUBROUTINE  SRECTG(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3,           02490000
     *                   U1,U2,U3,U11,U12,U21,U22,U31,U32)              02500000
      IMPLICIT REAL*8 (A-H,O-Z)                                         02510000
C                                                                       02520000
C*********************************************************************  02530000
C*****                                                           *****  02540000
C*****  INDEFINITE INTEGRAL OF SURFACE DISPLACEMENT,STRAIN,TILT  *****  02550000
C*****  DUE TO RECTANGULAR FAULT IN A HALF-SPACE                 *****  02560000
C*****                          CODED BY  Y.OKADA ... JAN 1985   *****  02570000
C*****                                                           *****  02580000
C*********************************************************************  02590000
C                                                                       02600000
C***** INPUT                                                            02610000
C*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)=1./((VP/VS)**2-1)   02620000
C*****   XI,ET,Q : FAULT COORDINATE                                     02630000
C*****   SD,CD   : SIN,COS OF DIP-ANGLE                                 02640000
C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)02650000
C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION      02660000
C                                                                       02670000
C***** OUTPUT                                                           02680000
C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL    )       02690000
C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /          02700000
C*****   U31,U32         : TILT                 UNIT OF XI,ET,Q )       02710000
C                                                                       02720000
      DATA  F0,F1,F2/ 0.D0, 1.D0, 2.D0 /                                02730000
      PI2=6.283185307179586D0                                           02740000
C-----                                                                  02750000
      XI2=XI*XI                                                         02760000
      ET2=ET*ET                                                         02770000
      Q2=Q*Q                                                            02780000
      R2=XI2+ET2+Q2                                                     02790000
      R =DSQRT(R2)                                                      02800000
      R3=R*R2                                                           02810000
      D =ET*SD-Q*CD                                                     02820000
      Y =ET*CD+Q*SD                                                     02830000
      RET=R+ET                                                          02840000
      IF(RET.LT.F0)  RET=F0                                             02850000
      RD =R+D                                                           02860000
      RRD=F1/(R*RD)                                                     02870000
C-----                                                                  02880000
      IF( Q .NE.F0)  TT = DATAN( XI*ET/(Q*R) )                          02890000
      IF( Q .EQ.F0)  TT = F0                                            02900000
      IF(RET.NE.F0)  RE = F1/RET                                        02910000
      IF(RET.EQ.F0)  RE = F0                                            02920000
      IF(RET.NE.F0)  DLE= DLOG(RET)                                     02930000
      IF(RET.EQ.F0)  DLE=-DLOG(R-ET)                                    02940000
      RRX=F1/(R*(R+XI))                                                 02950000
      RRE=RE/R                                                          02960000
      AXI=(F2*R+XI)*RRX*RRX/R                                           02970000
      AET=(F2*R+ET)*RRE*RRE/R                                           02980000
      IF(CD.EQ.F0)  GO TO 20                                            02990000
C==============================                                         03000000
C=====   INCLINED FAULT   =====                                         03010000
C==============================                                         03020000
      TD=SD/CD                                                          03030000
      X =DSQRT(XI2+Q2)                                                  03040000
      IF(XI.EQ.F0)  A5=F0                                               03050000
      IF(XI.NE.F0)                                                      03060000
     *A5= ALP*F2/CD*DATAN( (ET*(X+Q*CD)+X*(R+X)*SD) / (XI*(R+X)*CD) )   03070000
      A4= ALP/CD*( DLOG(RD) - SD*DLE )                                  03080000
      A3= ALP*(Y/RD/CD - DLE) + TD*A4                                   03090000
      A1=-ALP/CD*XI/RD        - TD*A5                                   03100000
      C1= ALP/CD*XI*(RRD - SD*RRE)                                      03110000
      C3= ALP/CD*(Q*RRE - Y*RRD)                                        03120000
      B1= ALP/CD*(XI2*RRD - F1)/RD - TD*C3                              03130000
      B2= ALP/CD*XI*Y*RRD/RD       - TD*C1                              03140000
      GO TO 30                                                          03150000
C==============================                                         03160000
C=====   VERTICAL FAULT   =====                                         03170000
C==============================                                         03180000
   20 RD2=RD*RD                                                         03190000
      A1=-ALP/F2*XI*Q/RD2                                               03200000
      A3= ALP/F2*( ET/RD + Y*Q/RD2 - DLE )                              03210000
      A4=-ALP*Q/RD                                                      03220000
      A5=-ALP*XI*SD/RD                                                  03230000
      B1= ALP/F2*  Q  /RD2*(F2*XI2*RRD - F1)                            03240000
      B2= ALP/F2*XI*SD/RD2*(F2*Q2 *RRD - F1)                            03250000
      C1= ALP*XI*Q*RRD/RD                                               03260000
      C3= ALP*SD/RD*(XI2*RRD - F1)                                      03270000
C-----                                                                  03280000
   30 A2=-ALP*DLE - A3                                                  03290000
      B3=-ALP*XI*RRE - B2                                               03300000
      B4=-ALP*( CD/R + Q*SD*RRE ) - B1                                  03310000
      C2= ALP*(-SD/R + Q*CD*RRE ) - C3                                  03320000
C-----                                                                  03330000
      U1 =F0                                                            03340000
      U2 =F0                                                            03350000
      U3 =F0                                                            03360000
      U11=F0                                                            03370000
      U12=F0                                                            03380000
      U21=F0                                                            03390000
      U22=F0                                                            03400000
      U31=F0                                                            03410000
      U32=F0                                                            03420000
C======================================                                 03430000
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 03440000
C======================================                                 03450000
      IF(DISL1.EQ.F0)  GO TO 200                                        03460000
      UN=DISL1/PI2                                                      03470000
      REQ=RRE*Q                                                         03480000
      U1 =U1 - UN*( REQ*XI +   TT    + A1*SD )                          03490000
      U2 =U2 - UN*( REQ*Y  + Q*CD*RE + A2*SD )                          03500000
      U3 =U3 - UN*( REQ*D  + Q*SD*RE + A4*SD )                          03510000
      U11=U11+ UN*( XI2*Q*AET - B1*SD )                                 03520000
      U12=U12+ UN*( XI2*XI*( D/(ET2+Q2)/R3 - AET*SD ) - B2*SD )         03530000
      U21=U21+ UN*( XI*Q/R3*CD + (XI*Q2*AET - B2)*SD )                  03540000
      U22=U22+ UN*( Y *Q/R3*CD + (Q*SD*(Q2*AET-F2*RRE)                  03550000
     *                            -(XI2+ET2)/R3*CD - B4)*SD )           03560000
      U31=U31+ UN*(-XI*Q2*AET*CD + (XI*Q/R3 - C1)*SD )                  03570000
      U32=U32+ UN*( D*Q/R3*CD + (XI2*Q*AET*CD - SD/R + Y*Q/R3 - C2)*SD )03580000
C===================================                                    03590000
C=====  DIP-SLIP CONTRIBUTION  =====                                    03600000
C===================================                                    03610000
  200 IF(DISL2.EQ.F0)  GO TO 300                                        03620000
      UN=DISL2/PI2                                                      03630000
      SDCD=SD*CD                                                        03640000
      U1 =U1 - UN*( Q/R             - A3*SDCD )                         03650000
      U2 =U2 - UN*( Y*Q*RRX + CD*TT - A1*SDCD )                         03660000
      U3 =U3 - UN*( D*Q*RRX + SD*TT - A5*SDCD )                         03670000
      U11=U11+ UN*( XI*Q/R3            + B3*SDCD )                      03680000
      U12=U12+ UN*( Y *Q/R3 - SD/R     + B1*SDCD )                      03690000
      U21=U21+ UN*( Y *Q/R3 + Q*CD*RRE + B1*SDCD )                      03700000
      U22=U22+ UN*( Y*Y*Q*AXI - (F2*Y*RRX + XI*CD*RRE)*SD + B2*SDCD )   03710000
      U31=U31+ UN*( D *Q/R3 + Q*SD*RRE + C3*SDCD )                      03720000
      U32=U32+ UN*( Y*D*Q*AXI - (F2*D*RRX + XI*SD*RRE)*SD + C1*SDCD )   03730000
C========================================                               03740000
C=====  TENSILE-FAULT CONTRIBUTION  =====                               03750000
C========================================                               03760000
  300 IF(DISL3.EQ.F0)  GO TO 900                                        03770000
      UN=DISL3/PI2                                                      03780000
      SDSD=SD*SD                                                        03790000
      U1 =U1 + UN*( Q2*RRE                       - A3*SDSD )            03800000
      U2 =U2 + UN*(-D*Q*RRX - SD*(XI*Q*RRE - TT) - A1*SDSD )            03810000
      U3 =U3 + UN*( Y*Q*RRX + CD*(XI*Q*RRE - TT) - A5*SDSD )            03820000
      U11=U11- UN*( XI*Q2*AET             + B3*SDSD )                   03830000
      U12=U12- UN*(-D*Q/R3 - XI2*Q*AET*SD + B1*SDSD )                   03840000
      U21=U21- UN*( Q2*(CD/R3 + Q*AET*SD) + B1*SDSD )                   03850000
      U22=U22- UN*((Y*CD-D*SD)*Q2*AXI - F2*Q*SD*CD*RRX                  03860000
     *                      - (XI*Q2*AET - B2)*SDSD )                   03870000
      U31=U31- UN*( Q2*(SD/R3 - Q*AET*CD) + C3*SDSD )                   03880000
      U32=U32- UN*((Y*SD+D*CD)*Q2*AXI + XI*Q2*AET*SD*CD                 03890000
     *                       - (F2*Q*RRX - C1)*SDSD )                   03900000
C-----                                                                  03910000
  900 RETURN                                                            03920000
      END                                                               03930000
