!> @mainpage @file DA2DF.f90
!! @brief  VoF wide averaged 2 dimensional model
!! @detail version 1.0.0
!! @detail F08 format
!! @date  2017.07.03
!! @date  Last Update 
!! @author 
program main
    USE classSpatialInfo
    USE classFluidModel
    USE classTimeIterator
    USE classOutputData
    IMPLICIT NONE
    CLASS(SpatialInfo), ALLOCATABLE :: spatial
    CLASS(fluidModel), ALLOCATABLE :: fluid
    CLASS(timeIterator), ALLOCATABLE :: time
    CLASS(outputData), ALLOCATABLE :: out

    DOUBLE PRECISION :: now, tend, toutstep, dt
    INTEGER :: nx, nz
    DOUBLE PRECISION :: dx, dz
    INTEGER :: is,ie,ke

    ! OPEN(119, FILE = 'error.log', STATUS = 'REPLACE')

    allocate(spatial)
    dx = 0.1 ; dz = 0.1 ; nx = int(7.5/dx) ; nz = int(5.0/dz)
    call spatial%inputSpatialInfo(nx,nz,dx,dz)
    
    ALLOCATE(fluid)
    CALL fluid%initialize( spatial%getSpatialInfo() )

!tmp 水柱の範囲id
    is = 1 ; ie = int(1.5/dx) ; ke = int(3.0/dz)
    CALL fluid%setInitialCondition(is,ie,ke) !set hydro static pressure

    dt = 0.0001
    now = 0.0
    tend = 10.0
    toutstep = 0.05

    ALLOCATE(time)
    CALL time%initialize(now, tend, toutstep)

    ALLOCATE(out)
    CALL out%initialize( )
    CALL out%outputALL( fluid%getSpatialInfo(), time%getTnow() )
    
    DO WHILE( time%getIsiterate() )

        CALL fluid%updateByHsmac( dt )
        CALL time%next( dt )

        IF( time%getIsOut() ) THEN
            CALL out%outputALL( fluid%getSpatialInfo(), time%getTnow() )
            print *, 'out', time%getTnow(), 'dt=', dt
        END IF

    END DO

    DEALLOCATE( out )
    DEALLOCATE( fluid )
    DEALLOCATE( time )

END PROGRAM