!> @file classSpatialInfo.f90
!! @brief 
!! @detail F08 format
!! @date 2017.7.7
!! @date Last Update
!! @author 
MODULE classSpatialInfo
    USE classCell
    IMPLICIT NONE
    PRIVATE

    TYPE, PUBLIC :: spatialInfo
        TYPE(cell), PRIVATE, POINTER :: c(:,:)
    CONTAINS
        PROCEDURE, PUBLIC :: inputSpatialInfo
        PROCEDURE, PUBLIC :: getSpatialInfo
    END TYPE
CONTAINS

!> @brief initialize
! @detail 
subroutine inputSpatialInfo(self,nx,nz,dx,dz) 
    CLASS(spatialInfo) :: self
    TYPE(cell), POINTER :: c(:,:)
    TYPE(cell), POINTER :: c_
    INTEGER, intent(in) :: nx, nz
    DOUBLE PRECISION ,intent(in) :: dx, dz
    INTEGER :: i, k
    DOUBLE PRECISION :: xcoord, zcoord

    ! nx = 150
    ! nz = 100
    ALLOCATE( c(1:nx, 1:nz) )
    ! dx = 0.05
    ! dz = 0.05

    DO k = 1, nz
    DO i = 1, nx
        xcoord = dx * (DBLE(i) - 0.5) 
        zcoord = dz * (DBLE(k) - 0.5)
        CALL c(i, k)%initialize( dx, dz, xcoord, zcoord )
    END DO
    END DO

    DO k = 1, nz
    DO i = 1, nx
        IF(i/=nx)THEN
            c_ => c(i+1,k)
            CALL c(i,k)%setNode( cellxp = c_ , BCxp = 0)
        ELSE
            CALL c(i,k)%setNode( BCxp = 1)
        END IF

        IF(i/=1)THEN
            c_ => c(i-1,k)
            CALL c(i,k)%setNode( cellxm = c_ , BCxm = 0)
        ELSE
            CALL c(i,k)%setNode( BCxm = 1)
        END IF

        IF(k/=nz)THEN 
            c_ => c(i,k+1)
            CALL c(i,k)%setNode( cellzp =  c_ , BCzp = 0)
        ELSE
            ! CALL c(i,k)%setNode( BCzp = 1) 
            CALL c(i,k)%setNode( BCzp = 2) 
        END IF

        IF(k/=1)THEN
            c_ => c(i,k-1)
            CALL c(i,k)%setNode( cellzm = c_ , BCzm = 0)
        ELSE
            CALL c(i,k)%setNode( BCzm = 1)
        END IF

        CALL c(i,k)%setInterfaceFlux()
    END DO
    END DO

    self%c => c
end subroutine


!> @brief initialize
! @detail 
function getSpatialInfo(self) result(c)
    CLASS(spatialInfo) :: self
    TYPE(cell), POINTER :: c(:,:)

    c => self%c
end function


END MODULE
