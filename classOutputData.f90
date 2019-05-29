!> @file classOutputData.f90
!! @brief 
!! @detail F08 format
!! @date 2017.2.18 
!! @date  Last Update 2017.7.7
!! @author 
MODULE classOutputData
    IMPLICIT NONE
    PRIVATE
    TYPE, PUBLIC :: outputData
        INTEGER, PRIVATE :: nserial
    CONTAINS
        PROCEDURE, PUBLIC :: initialize
        PROCEDURE, PUBLIC :: outputAll
    END TYPE
CONTAINS

!> @brief initialize
! @detail 
SUBROUTINE initialize(self)
    CLASS(outputData) :: self
    self%nserial = 0
END SUBROUTINE

!> @brief output
! @detail 
SUBROUTINE outputALL(self, c, t)
    USE classCell
    CLASS(outputData) :: self
    TYPE(cell), POINTER :: c(:,:)
    DOUBLE PRECISION, INTENT(IN) :: t
    INTEGER :: i, j
    CHARACTER(5) :: a
    CHARACTER(20) :: fn

    WRITE(a, '(I5.5)') self%nserial

    WRITE(fn, '(A12)') trim('out'//a//'.csv')
    OPEN(2001, FILE = fn, STATUS = 'REPLACE')
!header
    WRITE(2001, '( A2, "," ,E18.7 )' )'t=', t
    WRITE(2001, '( *(A, ",") )' ) 'I','J','xcoord','ycoord','pressure','uxplus','uxminus','uzplus', 'uzminus', 'alpha'

    DO j = LBOUND(c, 2), UBOUND(c, 2)
    DO i = LBOUND(c, 1), UBOUND(c, 1)
        WRITE(2001, '( I5, "," , I5 , "," , *(E20.10, ",") )') &
        i, j, c(i,j)%coord(1), c(i,j)%coord(2), c(i,j)%p, c(i,j)%Vxp, c(i,j)%Vxm, c(i,j)%Vzp, c(i,j)%Vzm, c(i,j)%alpha
    END DO
    END DO
    
    CLOSE(2001)

    self%nserial = self%nserial + 1
    nullify(c)
END SUBROUTINE

END MODULE
