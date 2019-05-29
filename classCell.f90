!> @file classCell.f90
!! @brief cellのクラス
!! @detail F08 format
!! @date 2017.07.03
!! @date  Last Update 
!! @author 
MODULE classCell
    IMPLICIT NONE
    PRIVATE
    INTEGER, PRIVATE, PARAMETER :: ndim = 2  ! number of dimension

    TYPE, PUBLIC :: cell
        TYPE(cell), POINTER :: cellxp => null() ! x-plus direction
        TYPE(cell), POINTER :: cellxm => null() ! x-minus direction
        ! TYPE(cell), POINTER :: yp => null() ! y-plus direction
        ! TYPE(cell), POINTER :: ym => null() ! y-minus direction
        TYPE(cell), POINTER :: cellzp => null() ! z-plus direction
        TYPE(cell), POINTER :: cellzm => null() ! z-minus direction

        INTEGER, ALLOCATABLE :: BCp(:), BCm(:) ! 0:通常境界、1:壁境界、2:開境界
        ! LOGICAL :: dsbound, usbound
        DOUBLE PRECISION, ALLOCATABLE :: dx(:) !cell size
        DOUBLE PRECISION, ALLOCATABLE :: coord(:) !cell center coordinate

        DOUBLE PRECISION, ALLOCATABLE :: p 
        DOUBLE PRECISION, ALLOCATABLE :: alpha
        DOUBLE PRECISION, POINTER :: Vxp => null()
        DOUBLE PRECISION, POINTER :: Vxm => null() 
        ! DOUBLE PRECISION, POINTER :: uyp => null() 
        ! DOUBLE PRECISION, POINTER :: uym => null() 
        DOUBLE PRECISION, POINTER :: Vzp => null()
        DOUBLE PRECISION, POINTER :: Vzm => null() 

!pointer for momentum Eq
        TYPE(cell), POINTER :: cellMainp => null()
        TYPE(cell), POINTER :: cellMainm => null() 
        TYPE(cell), POINTER :: cellSub1p => null()
        TYPE(cell), POINTER :: cellSub1m => null() 
        ! TYPE(cell), POINTER :: cellSub2p => null()
        ! TYPE(cell), POINTER :: cellSub2m => null()

        DOUBLE PRECISION, POINTER :: VMainp => null()
        DOUBLE PRECISION, POINTER :: VMainm => null() 
        DOUBLE PRECISION, POINTER :: VSub1p => null()
        DOUBLE PRECISION, POINTER :: VSub1m => null() 
        ! DOUBLE PRECISION, POINTER :: VSub2p => null()
        ! DOUBLE PRECISION, POINTER :: VSub2m => null()

        ! ! DOUBLE PRECISION, ALLOCATABLE :: qbed(:)
        ! DOUBLE PRECISION, ALLOCATABLE :: cManning
        ! LOGICAL :: isFluid

     CONTAINS
        PROCEDURE, PUBLIC :: initialize
        PROCEDURE, PUBLIC :: setInterfaceFlux
! ! !        PROCEDURE, PUBLIC  :: finalize
        PROCEDURE, PUBLIC :: setNode
        PROCEDURE, PUBLIC :: setCalculationDirection
!         PROCEDURE, PUBLIC :: setInitialCondition
    END TYPE
CONTAINS

! > @brief initialize
! @detail 
SUBROUTINE initialize(s, dx, dz, xcoor, zcoor)
    CLASS(cell) :: s
    DOUBLE PRECISION, INTENT(IN) :: dx, dz !cell delta
    DOUBLE PRECISION, INTENT(IN) :: xcoor, zcoor !cell center coordinate

    ALLOCATE( s%dx(ndim), s%coord(ndim) )
    ALLOCATE( s%BCp(ndim), s%BCm(ndim) )

    s%dx(1) = dx
    s%dx(2) = dz
    s%coord(1) = xcoor
    s%coord(2) = zcoor

    ALLOCATE(s%p, s%Vxp, s%Vxm, s%Vzp, s%Vzm, s%alpha)
    s%Vxp = 0.0
    s%Vzp = 0.0
    
END SUBROUTINE

!> @brief initialize
! @detail 
SUBROUTINE setNode(s, cellxp, cellxm, cellzp, cellzm, BCxp, BCxm, BCzp, BCzm)
    CLASS(cell) :: s
    TYPE(cell), POINTER, OPTIONAL :: cellxp,cellxm,cellzp,cellzm
    INTEGER, INTENT(IN), OPTIONAL :: BCxp, BCxm, BCzp, BCzm

    IF( PRESENT(cellxp) ) s%cellxp => cellxp
    IF( PRESENT(cellxm) ) s%cellxm => cellxm
    IF( PRESENT(cellzp) ) s%cellzp => cellzp
    IF( PRESENT(cellzm) ) s%cellzm => cellzm
    IF( PRESENT(BCxp) ) s%BCp(1) = BCxp
    IF( PRESENT(BCxm) ) s%BCm(1) = BCxm
    IF( PRESENT(BCzp) ) s%BCp(2) = BCzp
    IF( PRESENT(BCzm) ) s%BCm(2) = BCzm
END SUBROUTINE

! !> @brief initialize
! ! @detail セル境界の流速をポインタに割り当てる。
SUBROUTINE setInterfaceFlux(s)
    CLASS(cell) :: s
    
    IF( s%BCm(1) == 0 ) s%Vxm => s%cellxm%Vxp
    IF( s%BCm(2) == 0 ) s%Vzm => s%cellzm%Vzp

    IF( s%BCp(2) == 2 ) s%Vzp => s%Vzm

END SUBROUTINE

! !> @brief initialize
! ! @detail セル境界の流速をポインタに割り当てる。
SUBROUTINE setCalculationDirection(s, axis)
    CLASS(cell) :: s
    INTEGER, INTENT(IN) :: axis

    SELECT CASE(axis)
    CASE(1)
        s%cellMainp => s%cellxp
        s%cellMainm => s%cellxm
        s%cellSub1p => s%cellzp
        s%cellSub1m => s%cellzm

        s%VMainp => s%Vxp
        s%VMainm => s%Vxm
        s%VSub1p => s%Vzp
        s%VSub1m => s%Vzm
    
    CASE(2)
        s%cellMainp => s%cellzp
        s%cellMainm => s%cellzm
        s%cellSub1p => s%cellxp
        s%cellSub1m => s%cellxm

        s%VMainp => s%Vzp
        s%VMainm => s%Vzm
        s%VSub1p => s%Vxp
        s%VSub1m => s%Vxm

    END SELECT

END SUBROUTINE


END MODULE