!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION :: tau(2) = (/7.5,2.0/)
      DOUBLE PRECISION r2,v2,s2,z2, Delta2, C, etta, PI
      INTEGER i,j
      PI=3.141592653589793238462643383279502884197 

       v2=U(1)
       r2=U(2)
       s2=U(3)
       z2=U(4)

       etta =PAR(1)
       C =PAR(2)
       Delta2 = 1.0/(PI * tau(1))

       F(1)= (1/tau(1)) * (etta - (PI*r2*tau(1))**2 + v2**2) - C*s2
       F(2)= (1/tau(1)) * (Delta2 + 2 * r2 * v2)
       F(3)= z2/tau(2)
       F(4)= (r2 - 2*z2 - s2)/tau(2)

        IF(IJAC.EQ.0)RETURN

        do i = 1, 4
                do j = 1, 4
                        DFDU(i,j)=0;
                end do
        end do


        DFDU(1,1)=2*v2/tau(1);
        DFDU(1,2)=-(2*r2)*(PI**2)*tau(1);
        DFDU(1,4)= -C;

        DFDU(2,2)=2*r2/tau(1);
        DFDU(2,1)=2*v2/tau(1);

        DFDU(3,4)=1/tau(2);
   
        DFDU(4,2) = 1/tau(2);
        DFDU(4,3) = -1/tau(2);
        DFDU(4,4) = -2/tau(2);

        IF(IJAC.EQ.1)RETURN

        do i = 1, 4
                do j = 1, 4
                        DFDP(i,j)=0;
                end do
        end do

        DFDP(1,1)= 1/tau(1)
        DFDP(1,2)= -s2

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=0.0 ! P
      PAR(2)=0.0 ! EPSILON

      U(1)=0.0
      U(2)=0.0
      U(3)=0.0
      U(4)=0.0


      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT


!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUY_MAX(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: AUX_Y(0:NTST*NCOL)

      AUX_Y=U(2,:)

      GETUY_MAX=MAXVAL(AUX_Y)

      END FUNCTION GETUY_MAX
!-----------------------------------------------------------
!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUY_MIN(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: AUX_Y(0:NTST*NCOL)

      AUX_Y=U(2,:)

      GETUY_MIN=MINVAL(AUX_Y)

      END FUNCTION GETUY_MIN
!-----------------------------------------------------------

      SUBROUTINE PVLS(NDIM,U,PAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      DOUBLE PRECISION, EXTERNAL :: GETP,GETUY_MAX,GETUY_MIN
      INTEGER NDX,NCOL,NTST

      NDX=NINT(GETP('NDX',0,U))
      NTST=NINT(GETP('NTST',0,U))
      NCOL=NINT(GETP('NCOL',0,U))

     PAR(3)=GETUY_MIN(U,NDX,NTST,NCOL)
     PAR(4)=GETUY_MAX(U,NDX,NTST,NCOL)

      END SUBROUTINE PVLS
