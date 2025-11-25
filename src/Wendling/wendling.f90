!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   wendling : Wendling Neural Mass Model for a cortical neuron.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION sigm(v,vr)
       DOUBLE PRECISION  e0,r
       DOUBLE PRECISION :: v
       DOUBLE PRECISION :: vr
       r=0.56
       e0 = 2.5
       sigm = 2.0 * e0 / (1.0 + EXP(r * (vr - v)))
      END FUNCTION sigm
      DOUBLE PRECISION FUNCTION dsigm(v,vr)
       DOUBLE PRECISION  e0,r
       DOUBLE PRECISION :: v
       DOUBLE PRECISION :: vr
       r=0.56
       e0 = 2.5
       dsigm =2*e0*r*exp(r*(vr-v))/((1+exp(r*(vr-v)))**2)
      END FUNCTION dsigm


      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)


      DOUBLE PRECISION :: C(7) = (/135.0, 108.0, 33.75, 33.75, 40.5, 13.5, 108.0/)
      DOUBLE PRECISION :: A(5) = (/3.25, 3.25, 22.0, 10.0, 22.0/)
      DOUBLE PRECISION :: b(5) = (/100.0, 100.0, 50.0, 500.0, 50.0/) * 1.0
      DOUBLE PRECISION :: incoming(5) = (/0.0, 0.0, 0.0, 0.0, 0.0/)
      DOUBLE PRECISION :: w(5) = (/6.0, 6.0, 6.0, 6.0, 6.0/)
      DOUBLE PRECISION y0, y1, y2, y3, y4, y5, y6, y7, p1, p2
      DOUBLE PRECISION sigm,dsigm,s
      INTEGER i,j



       y0=U(1)
       y1=U(2)
       y2=U(3)
       y3=U(4)

       y4=U(5)
       y5=U(6)
       y6=U(7)
       y7=U(8)


       p1 =PAR(1) 


       F(1)= y4
       F(2)= y5
       F(3)= y6
       F(4)= y7

       
       F(5)= A(1) * b(1) * sigm(y1 - y2  - y3, w(1)) - 2 * b(1) * y4 - b(1)**2 * y0
       F(6)= A(2) * b(2) * (p1 - C(2)*sigm(C(1) * y0, w(2))) - 2 * b(2) * y5 - b(2)**2 * y1
       F(7)= A(3) * b(3) * C(4)*sigm(C(3) * y0, w(3)) - 2 * b(3) * y6 - b(3)**2 * y2
       F(8)= A(4) * b(4) * C(7)*sigm(C(5) * y3 -  C(6) * y2, w(4)) - 2 * b(4) * y7 - b(4)**2 * y3


        IF(IJAC.EQ.0)RETURN
        do i = 1, 8
                do j = 1, 8
                        DFDU(i,j)=0;
                end do
        end do


        DFDU(1,5)=1;
        DFDU(2,6)=1;
        DFDU(3,7)=1;
        DFDU(4,8)=1;

        DFDU(5,5)=-2*b(1);
        DFDU(6,6)=-2*b(2);
        DFDU(7,7)=-2*b(3);
        DFDU(8,8)=-2*b(4);

        DFDU(5,1)=-b(1)**2;
        DFDU(5,2)=A(1)*b(1)*dsigm(y1 -y2 - y3,w(1));
        DFDU(5,3)=-A(1)*b(1)*dsigm(y1 - y2 -y3 ,w(1));
        DFDU(5,4)=-A(1)*b(1)*dsigm(y1 - y2 - y3,w(1));

        DFDU(6,1)=-A(2)*b(2)*C(2)*C(1)*dsigm(C(1) * y0,w(2));
        DFDU(6,2)=-b(2)**2;

        DFDU(7,1)=A(3)*b(3)*C(4)*C(3)*dsigm(C(3) * y0,w(3));
        DFDU(7,3)=-b(3)**2;

        DFDU(8,3) = -A(4)*b(4)*C(7)*C(6)*dsigm(C(5) * y3 -  C(6) * y2 ,w(4)); 
        DFDU(8,4) = A(4)*b(4)*C(7)*C(5)*dsigm(C(5) * y3 -  C(6) * y2 ,w(4)) - b(4)**2;


        IF(IJAC.EQ.1)RETURN

        do i = 1, 8
                do j = 1, 1
                        DFDP(i,j)=0;
                end do
        end do

        DFDP(6,1) = A(2)*b(2)

      END SUBROUTINE FUNC
!-----------------------------------------------------------
!-----------------------------------------------------------
      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=00.0 ! P

      U(1)=0.01
      U(2)=0.01
      U(3)=0.01
      U(4)=0.01
      U(5)=0.01
      U(6)=0.01
      U(7)=0.01
      U(8)=0.01

      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      DOUBLE PRECISION, EXTERNAL :: GETP,GETUMN1, GETUMX1
       PAR(2)= GETP('MIN',1,U)!GETUMN1(U,NDX,NTST,NCOL) 
       PAR(3)= GETP('MAX',1,U)!GETUMX1(U,NDX,NTST,NCOL)
      END SUBROUTINE PVLS
