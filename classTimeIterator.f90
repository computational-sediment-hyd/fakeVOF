!> @file classTimeIterator.f90
!! @brief time Iterator
!! @detail F03 format
!! @date 2017.2.18
!! @date  Last Update 
!! @author 
MODULE classTimeIterator
    IMPLICIT NONE
    PRIVATE

    TYPE, PUBLIC :: timeIterator
        DOUBLE PRECISION, PRIVATE :: tend
        DOUBLE PRECISION, PRIVATE :: now
        DOUBLE PRECISION, PRIVATE :: tout
        DOUBLE PRECISION, PRIVATE :: toutstep
        LOGICAL, PRIVATE :: isIterate
        LOGICAL, PRIVATE :: isOut
    CONTAINS
        PROCEDURE, PUBLIC :: initialize
        PROCEDURE, PUBLIC :: next
        PROCEDURE, PUBLIC :: getTnow
        PROCEDURE, PUBLIC :: getIsiterate
        PROCEDURE, PUBLIC :: getIsOut
    END TYPE
CONTAINS

!> @brief initialize
! @detail 
SUBROUTINE initialize( self, now, tend, toutstep )
    CLASS(timeIterator) :: self
    DOUBLE PRECISION, INTENT(IN) :: now, tend, toutstep
    self%now = now
    self%tend = tend
    self%toutstep = toutstep
    self%tout = now + toutstep
    self%isIterate = .TRUE.
    self%isOut = .FALSE.
END SUBROUTINE

!> @brief time march
! @detail 
SUBROUTINE next(self, dt )
    CLASS(timeIterator) :: self
    DOUBLE PRECISION, INTENT(IN) :: dt
    self%now = self%now + dt
    IF( self%tend <= self%now ) self%isIterate = .FALSE.
    IF( self%tout <= self%now ) THEN
        self%tout = self%tout + self%toutstep
        self%isOut = .TRUE.
    ELSE
        self%isOut = .FALSE.
    END IF
END SUBROUTINE

!> @brief getter
! @detail 
FUNCTION getTnow(self) RESULT(res)
    CLASS(timeIterator) :: self
    DOUBLE PRECISION :: res
    res = self%now
END FUNCTION

!> @brief getter
! @detail 
FUNCTION getIsiterate(self) RESULT(res)
    CLASS(timeIterator) :: self
    logical :: res
    res = self%isIterate
END FUNCTION

!> @brief getter
! @detail 
FUNCTION getIsOut(self) RESULT(res)
    CLASS(timeIterator) :: self
    logical :: res
    res = self%isOut
END FUNCTION

END MODULE

