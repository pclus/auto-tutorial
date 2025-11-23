!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Here we study the bifurcation diagram of the Wilson-Cowan model.    
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION x,y,w_xx,w_yy,w_xy,w_yx,I_x,I_y,tau_y,tau_x	
      DOUBLE PRECISION e0,r,r_x,r_y,v0_x,v0_y,c,sigm,s,a_x,a_y, dsigm
      INTEGER i,j


       sigm(s)=1 / (1+exp(r*(v0_x - s)))
       dsigm(s)=r*exp(r*(v0_x-s))/((1+exp(r*(v0_x-s)))**2)

!       w_xx = 12
       w_xy = 12   
       w_yx = 12
       w_yy = 2
       tau_x = 1   
       tau_y = 1
       r = 1
       v0_x = 0
       v0_y = 0

       x=U(1)
       y=U(2)
       I_x =PAR(1) 
       I_y =PAR(2) 
       w_xx =PAR(7)

       F(1)= (-x + sigm(w_xx*x - w_xy*y + I_x) )/tau_x
       F(2)= (-y + sigm(w_yx*x - w_yy*y + I_y) )/tau_y


      IF(IJAC.EQ.0)RETURN

        do i = 1, 2
                do j = 1, 2
                        DFDU(i,j)=0;
                end do
        end do

        DFDU(1,1)= (-1 + w_xx*dsigm(w_xx*x - w_xy*y + I_x) )/tau_x ;

        DFDU(1,2)= (-w_xy*dsigm(w_xx*x - w_xy*y + I_x) )/tau_x;

        DFDU(2,1)= (w_yx*dsigm(w_yx*x - w_yy*y + I_y))/tau_y;

        DFDU(2,2)= (-1 - w_yy*dsigm(w_yx*x - w_yy*y + I_y))/tau_y ;


      IF(IJAC.EQ.1)RETURN

        do i = 1, 2
                do j = 1, 3
                        DFDP(i,j)=0;
                end do
        end do

       DFDP(1,1)=(dsigm(w_xx*x - w_xy*y + I_x) )/tau_x
       DFDP(2,2)=(dsigm(w_yx*x - w_yy*y + I_y) )/tau_y
       DFDP(1,3)=(dsigm(w_xx*x - w_xy*y + I_x)*x)/tau_x


      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=0.0 ! mu
      PAR(2)=0.0 ! mu
      PAR(7)=0.0
      U(1)=7.7E-02
      U(2)=-6.05E-02

      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT


!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUX_MAX(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: AUX_Y(0:NTST*NCOL)

      AUX_Y=U(1,:)

      GETUX_MAX=MAXVAL(AUX_Y)

      END FUNCTION GETUX_MAX
!-----------------------------------------------------------
!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUX_MIN(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: AUX_Y(0:NTST*NCOL)

      AUX_Y=U(1,:)

      GETUX_MIN=MINVAL(AUX_Y)

      END FUNCTION GETUX_MIN
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
      DOUBLE PRECISION, EXTERNAL :: GETP,GETUX_MAX,GETUX_MIN,GETUY_MAX,GETUY_MIN
      INTEGER NDX,NCOL,NTST
      DOUBLE PRECISION U1, ONE,r,v0

      NDX=NINT(GETP('NDX',0,U))
      NTST=NINT(GETP('NTST',0,U))
      NCOL=NINT(GETP('NCOL',0,U))
      PAR(3)=GETUX_MIN(U,NDX,NTST,NCOL)
      PAR(4)=GETUX_MAX(U,NDX,NTST,NCOL)
      PAR(5)=GETUY_MIN(U,NDX,NTST,NCOL)
      PAR(6)=GETUY_MAX(U,NDX,NTST,NCOL)

      END SUBROUTINE PVLS
