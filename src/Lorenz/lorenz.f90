!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Here we study the bifurcation diagram of the classic Lorenz model.    
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)


      DOUBLE PRECISION x,y,z,b,sigma,r
      INTEGER i,j
       b=8/3
       sigma=10
       x=U(1)
       y=U(2)
       z=U(3)
       r=PAR(1) 

       F(1)= sigma*(y - x)
       F(2)= r*x - y - x*z
       F(3)= x*y - b*z


      IF(IJAC.EQ.0)RETURN

        do i = 1, 3
                do j = 1, 3
                        DFDU(i,j)=0;
                end do
        end do

        DFDU(1,1)=-sigma;
        DFDU(1,2)=sigma;
        DFDU(2,1)=r - z;
        DFDU(2,2)=-1;
        DFDU(2,3)=-x;
        DFDU(3,1)=y;
        DFDU(3,2)=x;
        DFDU(3,3)=-b;


      IF(IJAC.EQ.1)RETURN

        do i = 1, 3
                do j = 1, 1
                        DFDP(i,j)=0;
                end do
        end do

       DFDP(2,1)=x

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=0.0 ! mu

      U(1)=7.7E-02
      U(2)=-6.05E-02
      U(3)=-6.05E-02


      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUY_MAX(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: AUX_Y(0:NTST*NCOL)

      AUX_Y=U(1,:)

      GETUY_MAX=MAXVAL(AUX_Y)

      END FUNCTION GETUY_MAX
!-----------------------------------------------------------
!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUY_MIN(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: AUX_Y(0:NTST*NCOL)

      AUX_Y=U(1,:)

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
      PAR(2)=GETUY_MIN(U,NDX,NTST,NCOL)
      PAR(3)=GETUY_MAX(U,NDX,NTST,NCOL)
      END SUBROUTINE PVLS
