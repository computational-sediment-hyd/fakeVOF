!> @file classFluidModel.f90
!! @brief cellのクラス
!! @detail F08 format
!! @date 2017.2.15
!! @date  Last Update 
!! @author 
MODULE classFluidModel
    USE classCell
    IMPLICIT NONE
    PRIVATE
    DOUBLE PRECISION, PRIVATE, PARAMETER :: rhoLiquid = 1000.0 !重力加速度
    DOUBLE PRECISION, PRIVATE, PARAMETER :: rhoGas = 1.0 !重力加速度
    DOUBLE PRECISION, PRIVATE, PARAMETER :: nuLiquid = 0.000001 !動粘性係数：水
    DOUBLE PRECISION, PRIVATE, PARAMETER :: nuGas = 0.00001 !動粘性係数：気体
    DOUBLE PRECISION, PRIVATE, PARAMETER :: gravity = 9.8 !重力加速度
    INTEGER, PRIVATE, PARAMETER :: ndim = 2  ! number of dimension

    type :: celltmp
        DOUBLE PRECISION :: Vxp,Vxm,Vzp,Vzm
        DOUBLE PRECISION :: alpha
    end type

    TYPE, PUBLIC :: fluidModel
        TYPE(cell), PRIVATE, POINTER :: cells(:,:)
    CONTAINS
        PROCEDURE, PUBLIC  :: initialize
        PROCEDURE, PUBLIC  :: setInitialCondition
        PROCEDURE, PUBLIC  :: poissonEq
        PROCEDURE, PUBLIC  :: momentumEq
        PROCEDURE, PUBLIC  :: updateByHsmac
        PROCEDURE, PUBLIC  :: getSpatialInfo
        PROCEDURE, PUBLIC  :: vofFunction

!        PROCEDURE, PUBLIC  :: finalize
    END TYPE

CONTAINS

!> @brief initialize
! @detail 
SUBROUTINE initialize(self, c)
    CLASS(fluidModel) :: self
    TYPE(cell), POINTER :: c(:,:)

    self%cells => c
END SUBROUTINE

!> @brief set Initial Condition
! @detail hydro static pressure 
SUBROUTINE setInitialCondition(self,is,ie,ke)
    CLASS(fluidModel) :: self
    TYPE(cell), POINTER :: cp
    INTEGER,INTENT(in) :: is,ie,ke
    INTEGER :: i, k
    DOUBLE PRECISION :: depth

    depth = 0.0
    ! DO k = 1, UBOUND(self%cells, 2)
    !     depth = depth + self%cells(1,k)%dx(2)
    ! END DO

    ! DO k = 1, UBOUND(self%cells, 2)
    !     DO i = 1, UBOUND(self%cells, 1)
    !         cp => self%cells(i,k)
    !         cp%alpha = 1.0
    !         cp%p = 0.0 !rhoLiquid * gravity * ( depth - cp%coord(2) )
    !     END DO
    ! END DO

    DO k = 1, UBOUND(self%cells, 2)
    DO i = 1, UBOUND(self%cells, 1)
        cp => self%cells(i,k)
        cp%alpha = 0.0
        cp%p = 0.0 !rhoLiquid * gravity * ( depth - cp%coord(2) )
    END DO
    END DO

    depth = self%cells(1,1)%dx(2)*ke
    DO k = 1, ke
    DO i = is, ie
        cp => self%cells(i,k)
        cp%alpha = 1.0
        cp%p = rhoLiquid * gravity * ( depth - cp%coord(2) )
    END DO
    END DO


    DO k = 1, UBOUND(self%cells, 2)
        DO i = 1, UBOUND(self%cells, 1)
            cp => self%cells(i, k)
            cp%Vxp = 0.0
            cp%Vxm = 0.0
            cp%Vzp = 0.0
            cp%Vzm = 0.0
        END DO
    END DO

END SUBROUTINE

!> @brief initialize
! @detail 
function getSpatialInfo(self) result(c)
    CLASS(fluidModel) :: self
    TYPE(cell), POINTER :: c(:,:)

    c => self%cells
end function


!> @brief hsmac Equation
!! @detail 
SUBROUTINE updateByHsmac(self, dt)
    CLASS(fluidModel) :: self
    DOUBLE PRECISION, INTENT(IN) :: dt
    INTEGER :: i, k, axis
    DOUBLE PRECISION :: ks, vis
    TYPE(celltmp), allocatable :: cTmp(:,:)

!make tmp array
    ALLOCATE( cTmp( UBOUND(self%cells, 1), UBOUND(self%cells, 2) ) )
    DO k = 1, UBOUND(self%cells, 2)
        DO i = 1, UBOUND(self%cells, 1)
            cTmp(i,k)%Vxp = self%cells(i,k)%Vxp
            cTmp(i,k)%Vzp = self%cells(i,k)%Vzp
            if(i == 1) cTmp(i,k)%Vxm = self%cells(i, k)%Vxm
            if(k == 1) cTmp(i,k)%Vzm = self%cells(i, k)%Vzm
            cTmp(i,k)%alpha = self%cells(i,k)%alpha
        END DO
    END DO

! update alfa
    call self%vofFunction(cTmp, dt)

    ks = 0.000000001
    vis = 0.00000001
    axis = 1
    CALL self%momentumEq( axis, cTmp, dt, ks, vis )
    axis = 2
    CALL self%momentumEq( axis, cTmp, dt, ks, vis )


! update variable
    DO k = 1, UBOUND(self%cells, 2)
        DO i = 1, UBOUND(self%cells, 1)
            self%cells(i, k)%Vxp = cTmp(i, k)%Vxp
            self%cells(i, k)%Vzp = cTmp(i, k)%Vzp
            IF(i == 1) self%cells(i, k)%Vxm = cTmp(i, k)%Vxm
            IF(k == 1) self%cells(i, k)%Vzm = cTmp(i, k)%Vzm
            self%cells(i,k)%alpha = cTmp(i,k)%alpha
        END DO
    END DO

   CALL self%poissonEq( dt )

END SUBROUTINE

!> @brief momentum Equation
!! @detail 正方向のみを計算する。
SUBROUTINE vofFunction(self, cTmp, dt)
    CLASS(fluidModel) :: self
    INTEGER :: axis !1:x,2:z
    INTEGER :: i, k
    DOUBLE PRECISION, INTENT(IN) :: dt
    DOUBLE PRECISION :: qp(ndim), qm(ndim)
    TYPE(cell), POINTER :: c => null()
    TYPE(celltmp), allocatable :: cTmp(:,:)


    DO k = 1, UBOUND(self%cells, 2)
    DO i = 1, UBOUND(self%cells, 1)
        c => self%cells(i,k)

        do axis = 1, ndim
            call c%setCalculationDirection(axis)

            if( c%BCp(axis) == 1 )THEN
                qp(axis) = 0.0
            else if( c%BCp(axis) == 2 )THEN
                qp(axis) = c%Vmainp * c%alpha
            else IF( c%Vmainp > 0.0 )THEN
                qp(axis) = c%Vmainp * c%alpha
            ELSE
                qp(axis) = c%Vmainp * c%cellMainp%alpha
            END IF

            if( c%BCm(axis) == 1 )THEN
                qm(axis) = 0.0
            else if( c%BCm(axis) == 2 )THEN
                qm(axis) = c%Vmainm * c%alpha
            else IF( c%Vmainm > 0.0 )THEN
                qm(axis) = c%Vmainm * c%cellMainm%alpha
            ELSE
                qm(axis) = c%Vmainm * c%alpha
            END IF

        end do

        cTmp(i,k)%alpha = c%alpha - dt * SUM( ( qp(:)-qm(:) ) / c%dx(:) )
        if( cTmp(i,k)%alpha < 10.0**(-20.0)) cTmp(i,k)%alpha = 0.0

    end do
    end do
END SUBROUTINE

!> @brief Poisson Equation
! @detail 
SUBROUTINE poissonEq(self, dt)
    CLASS(fluidModel) :: self
    INTEGER :: i, k, n, nd
    DOUBLE PRECISION, INTENT(IN) :: dt 
    DOUBLE PRECISION :: beta, dxp(ndim), dxm(ndim), dp, dudx, deltamax
    TYPE(cell), POINTER :: cp
    DOUBLE PRECISION, POINTER :: ufluxp, ufluxm
   
    beta = 1.7
    deltamax = 1.0
    n = 1

    DO WHILE ( (n < 100000) .AND. (deltamax > 10.0**(-8.0)) )
        DO k = 1, UBOUND(self%cells, 2)
        DO i = 1, UBOUND(self%cells, 1)
            cp => self%cells(i,k)

            DO nd = 1, ndim
                IF( cp%BCp(nd) /= 0 )THEN
                    dxp(nd) = cp%dx(nd)
                ELSE
                    SELECT CASE(nd)
                        CASE(1) !x
                            dxp(nd) = 0.5*( cp%dx(nd)+cp%cellxp%dx(nd) )
                        CASE(2) !z
                            dxp(nd) = 0.5*( cp%dx(nd)+cp%cellzp%dx(nd) )
                    END SELECT
                END IF

                IF( cp%BCm(nd) /= 0 )THEN
                    dxm(nd) = cp%dx(nd)
                ELSE
                     SELECT CASE(nd)
                        CASE(1) !x
                            dxm(nd) = 0.5*(cp%dx(nd)+cp%cellxm%dx(nd))
                        CASE(2) !z
                            dxm(nd) = 0.5*(cp%dx(nd)+cp%cellzm%dx(nd))
                    END SELECT
                END IF
            END DO

            dudx = ( cp%Vxp-cp%Vxm ) / cp%dx(1) &
                    + ( cp%Vzp-cp%Vzm ) / cp%dx(2)

            dp = -beta*dudx / dt / ( &
            + (1.0/dxp(1)+1.0/dxm(1))/cp%dx(1) & 
            + (1.0/dxp(2)+1.0/dxm(2))/cp%dx(2) &
            )

            cp%p = cp%p + dp

            DO nd = 1, ndim
                SELECT CASE(nd)
                    CASE(1) !x
                        ufluxp => cp%Vxp
                        ufluxm => cp%Vxm
                    CASE(2) !y
                        ufluxp => cp%Vzp
                        ufluxm => cp%Vzm
                END SELECT
                ufluxp = MERGE( ufluxp + dp*dt/dxp(nd), ufluxp, cp%BCp(nd) == 0 )
                ufluxm = MERGE( ufluxm - dp*dt/dxm(nd), ufluxm, cp%BCm(nd) == 0 )
            END DO
        END DO
        END DO

        deltamax = 0.0
        DO k = 1, UBOUND(self%cells, 2)
        DO i = 1, UBOUND(self%cells, 1)
            cp => self%cells(i,k)
            dudx = ABS( &
                    ( cp%Vxp-cp%Vxm ) / cp%dx(1) &
                  + ( cp%Vzp-cp%Vzm ) / cp%dx(2) &
            )
            IF(dudx*dt > deltamax) deltamax = dudx*dt
        END DO
        END DO

        n = n + 1
    END DO
    IF(n >= 99999) print *, 'poisson error', n

END SUBROUTINE

!> @brief Umomentum Equation
!! @detail 正方向のみを計算する。
SUBROUTINE momentumEq(self, axis, cTmp, dt, ks, vis)
    CLASS(fluidModel) :: self
    INTEGER, INTENT(IN) :: axis !1:x,2:z
    INTEGER :: i, k, n, axSub1
    DOUBLE PRECISION, INTENT(IN) :: dt, ks, vis
    DOUBLE PRECISION :: up(ndim), um(ndim)
    DOUBLE PRECISION :: advection(ndim), qp, qm, dx(ndim)
    DOUBLE PRECISION :: pressure, shstp(ndim), shstm(ndim), shst, cf, exforce, rho_, nu_
    TYPE(cell), POINTER :: c => null()
    TYPE(celltmp), allocatable :: cTmp(:,:)

    axSub1 = 1
    exforce=0.0
    rho_ = 1000.0
!計算する方向をaxisとして、その他をaxSubと設定
    SELECT CASE(axis)
    CASE(1)
        axSub1 = 2
    CASE(2)
        axSub1 = 1
    END SELECT

!重要なポイント：計算する軸方向に合わせて、変数名を変更
! axis=1はx方向の計算となり、変数名のMain：x方向、Sub1：z方向
! axis=2はz方向の計算となり、変数名のMain：z方向、Sub1：x方向

    DO k = 1, UBOUND(self%cells, 2)
    DO i = 1, UBOUND(self%cells, 1)
        c => self%cells(i,k)
        call c%setCalculationDirection(axis)
    END DO
    END DO

    DO k = 1, UBOUND(self%cells, 2)
    DO i = 1, UBOUND(self%cells, 1)
        c => self%cells(i,k)
! Wall Boundary condition
        IF( c%BCp(axis) == 1 )THEN
            SELECT CASE(axis)
            CASE(1)
                cTmp(i,k)%Vxp = 0.0
            CASE(2)
                cTmp(i,k)%Vzp = 0.0
            END SELECT
            CYCLE
        END IF

! free boundary condition
!ポインタにより0勾配を実装
        IF( c%BCp(axis) == 2 ) cycle

        up(axis) = 0.5*( c%cellMainp%VMainp + c%cellMainp%VMainm )
        um(axis) = 0.5*( c%VMainp + c%VMainm )
        up(axSub1) = 0.5*( c%VSub1p + c%cellMainp%VSub1p )
        um(axSub1) = 0.5*( c%VSub1m + c%cellMainp%VSub1m )
    
        dx(axis) = 0.5*( c%dx(axis) + c%cellMainp%dx(axis) )
        dx(axSub1) = c%dx( axSub1 )

!** 移流項の計算 ***-------------------------------------------------
        n = axis
        qp = merge( rho(c%cellMainp%alpha)*c%cellMainp%VMainm &
                   ,rho(c%cellMainp%alpha)*c%cellMainp%VMainp, up(n) > 0.0) 
        qm = merge( rho(c%alpha)*c%VMainm, rho(c%alpha)*c%VMainp, um(n) > 0.0)
        advection(n) = ( up(n)*qp - um(n)*qm ) / dx(n)

        n = axSub1
!境界セルのいづれかが壁面の場合、移流は0とする。
        IF( (c%BCp(n)==1) .OR. (c%cellMainp%BCp(n)==1) )THEN
            qp = 0.0
        else IF( (c%BCp(n)==2) .OR. (c%cellMainp%BCp(n)==2) )THEN
! 開境界の場合：外側と内側で同じプロファイルとして界面流速の正負に関わらず同じ値を与える。
            rho_ = 0.5*( rho(c%cellMainp%alpha) + rho(c%alpha) )
            qp = rho_ * c%VMainp
        else IF(up(n) > 0.0) THEN
            rho_ = 0.5*( rho(c%cellMainp%alpha) + rho(c%alpha) )
            qp = rho_ * c%VMainp
        else
            rho_ = 0.5*( rho( c%cellSub1p%alpha ) + rho( c%cellSub1p%cellMainp%alpha ))
            qp = rho_ * c%cellSub1p%VMainp
        END IF

!境界セルのいづれかが壁面の場合、移流は0とする。
        IF( (c%BCm(n)==1) .OR. (c%cellMainp%BCm(n)==1) )THEN
            qm = 0.0
        else IF( (c%BCm(n)==2) .OR. (c%cellMainp%BCm(n)==2) )THEN
! 開境界の場合：外側と内側で同じプロファイルとして界面流速の正負に関わらず同じ値を与える。
            rho_ = 0.5*( rho( c%cellSub1m%alpha ) + rho( c%cellSub1m%cellMainp%alpha ))
            qm = rho_ * c%cellSub1m%VMainp
        else IF(um(n) > 0.0) THEN
            rho_ = 0.5*( rho( c%cellSub1m%alpha ) + rho( c%cellSub1m%cellMainp%alpha ))
            qm = rho_ * c%cellSub1m%VMainp
        ELSE
            rho_ = 0.5*( rho(c%cellMainp%alpha) + rho(c%alpha) )
            qm = rho_ * c%VMainp
        END IF

        advection(n) = ( up(n)*qp - um(n)*qm ) / dx(n)

!** 圧力項(分離解法)の計算 ***------------------------------------------------
        pressure = ( c%cellMainp%p - c%p )/ dx(axis)  !/ rho

!** 重力項の計算 ***----------------------------------------------------------
        SELECT CASE(axis)
        CASE(1)
            exforce = 0.0
        CASE(2)
            rho_ = 0.5*( rho(c%cellMainp%alpha) + rho(c%alpha) )
            exforce = - rho_ * gravity
        END SELECT

!** 粘性項の計算 ***----------------------------------------------------------
        n = axis
        shstp(n) = rho(c%cellMainp%alpha)*nu(c%cellMainp%alpha) &
        *(c%cellMainp%VMainp - c%cellMainp%VMainm) / c%cellMainp%dx(n)
        shstm(n) = rho(c%alpha)*nu(c%alpha)*(c%VMainp - c%VMainm) / c%dx(n)

        n = axSub1
        cf  = ( log( 30.0*0.5*dx(n)/ks) / 0.4 )**(-2.0)

        IF( (c%BCp(n)==1) .OR. (c%cellMainp%BCp(n)==1) )THEN
            rho_ = 0.5*( rho(c%cellMainp%alpha) + rho(c%alpha) )
            shstp(n) = rho_ * cf * c%VMainp**2.0 
        else IF( (c%BCp(n)==2) .OR. (c%cellMainp%BCp(n)==2) )THEN
            shstp(n) = 0.0
        ELSE
            rho_ = 0.25*( rho(c%alpha) + rho(c%cellMainp%alpha) &
            + rho(c%cellSub1p%alpha) + rho(c%cellSub1p%cellMainp%alpha) )
            nu_ = 0.25*( nu(c%alpha) + nu(c%cellMainp%alpha) &
            + nu(c%cellSub1p%alpha) + nu(c%cellSub1p%cellMainp%alpha) )

            shstp(n) = rho_ * nu_*(c%cellSub1p%VMainp - c%VMainp) / 0.5 / (c%dx(n)+c%cellSub1p%dx(n))
        END IF


        IF( (c%BCm(n)==1) .OR. (c%cellMainp%BCm(n)==1) )THEN
            rho_ = 0.5*( rho(c%cellMainp%alpha) + rho(c%alpha) )
            shstm(n) = rho_ *  cf * c%VMainp**2.0
        else IF( (c%BCm(n)==2) .OR. (c%cellMainp%BCm(n)==2) )THEN
            shstm(n) = 0.0
        ELSE
            rho_ = 0.25*( rho(c%alpha) + rho(c%cellMainp%alpha) &
            + rho(c%cellSub1m%alpha) + rho(c%cellSub1m%cellMainp%alpha) )
            nu_ = 0.25*( nu(c%alpha) + nu(c%cellMainp%alpha) &
            + nu(c%cellSub1m%alpha) + nu(c%cellSub1m%cellMainp%alpha) )
            
            shstm(n) = rho_ * nu_*(c%VMainp - c%cellSub1m%VMainp) / 0.5 / (c%dx(n)+c%cellSub1m%dx(n))
        END IF

!tmporary:watersurface 2.0m
        ! IF(axis == 1)THEN
        !     rho_ = 0.5 * ( rho( c%cellMainp%alpha ) + rho( c%alpha ) )
        !     IF(k == UBOUND(self%cells, 2)) &
        !     shstp(n) = rho_ * vis*(2.0 - c%VMainp)/c%dx(n)
        ! END IF
!------------------------------------------------------------------

        shst = SUM( (shstp(:) - shstm(:)) / dx(:) )

        SELECT CASE(axis)
        CASE(1)
            cTmp(i,k)%Vxp = ( &
            rho(c%alpha) * c%VMainp &
            + dt*( -SUM(advection(:)) + exforce - pressure + shst ) &
            ) / rho( cTmp(i,k)%alpha )
        CASE(2)
            cTmp(i,k)%Vzp = ( &
            rho(c%alpha)*c%VMainp &
            + dt*( -SUM(advection(:)) + exforce - pressure + shst ) &
            ) / rho( cTmp(i,k)%alpha )
        END SELECT

    END DO
    END DO

!boundary 1st cell minus direciton Velocity
    SELECT CASE(axis)
    CASE(1)
        DO k = 1, UBOUND(self%cells, 2)
            IF( self%cells(1,k)%BCm(axis) == 1 ) cTmp(1,k)%Vxm = 0.0
        END DO
    CASE(2)
        DO i = 1, UBOUND(self%cells, 1)
            IF( self%cells(i,1)%BCm(axis) == 1 ) cTmp(i,1)%Vzm = 0.0
        END DO
    END SELECT

CONTAINS

    function rho(alpha) result(res)
        DOUBLE PRECISION, intent(in) :: alpha
        DOUBLE PRECISION :: res
        res = alpha * rhoLiquid + (1.0 - alpha)*rhoGas
    end function

    function nu(alpha) result(res)
        DOUBLE PRECISION, intent(in) :: alpha
        DOUBLE PRECISION :: res, mu
        mu = alpha * rhoLiquid*nuLiquid + (1.0 - alpha)*rhoGas*nuGas
        res = mu / rho(alpha)
    end function

END SUBROUTINE

END MODULE
