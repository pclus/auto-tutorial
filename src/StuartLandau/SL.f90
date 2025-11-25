!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Here we study the bifurcation diagram of the Stuart-Landau model.    
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION x,y,a,w
      INTEGER i,j
   

       x=U(1)
       y=U(2)
       a =PAR(1) 
       w =PAR(2) 


       F(1)= (a - (x**2 + y**2))*x - w*y 
       F(2)= (a - (x**2 + y**2))*y + w*x 


      IF(IJAC.EQ.0)RETURN

        do i = 1, 2
                do j = 1, 2
                        DFDU(i,j)=0;
                end do
        end do

        DFDU(1,1)= a - 3*x**2 ;
        DFDU(1,2)= -y**2 - w;
        DFDU(2,1)= -x**2 + w;
        DFDU(2,2)= a - 3*y**2 ;


      IF(IJAC.EQ.1)RETURN

        do i = 1, 2
                do j = 1, 1
                        DFDP(i,j)=0;
                end do
        end do

       DFDP(2,1)=(1-x**2)*y

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=0.0 ! mu
      PAR(2)=0.0 ! mu

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
!-----------------------------------------------------------

      SUBROUTINE PVLS(NDIM,U,PAR)
      END SUBROUTINE PVLS
