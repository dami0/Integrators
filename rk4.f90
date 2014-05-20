! Released under the CC0 1.0 Universal (CC0 1.0) 
! A copy of it has been included in this repository as the file "LICENSE"

! Damian A. Wloch. University of Surrey, FEPS. Department of Physics.
! Runge-Kutta 4th order integrator with distributed potential and reference frame forces

PROGRAM integrator
    IMPLICIT NONE

    INTEGER :: j, m, l ,i
    DOUBLE PRECISION :: irg, nrg, energy
    DOUBLE PRECISION ::  q, dt, w
    DOUBLE PRECISION, DIMENSION(1:3) :: pos = 0.0d+0, vel = 0.0d+0

!   ____Initial Conditions_____________________________________________________
    i  = 5                      !dt is 2 raised to -i
    l  = 2**(i - 1)             !data output at half-t
    dt = 2**(-REAL(i))          !dt assigned
    m  = 100000                 !number of time units to simulate for
    m  = m*INT(dt**(-1))        !divide no. of t-u by dt to get run time in int
    q  = 0.58904862d+0          !virial core radius
    w  = 0.04921478110900686d+0 !angular speed (centrifugal force)

    READ(5,*) pos, vel !defined by stdin ( <values> | ./this > output.file )

    nrg = energy(pos, vel, w, q)
    WRITE(6,'(5e14.6)') pos, SQRT(SUM(vel*vel)), nrg

!   ____Iterate rk4/energy_____________________________________________________
    DO j=1,m
        IF (MODULO(j, l) == 0) THEN     !cut down I/O operations
            nrg = energy(pos, vel, w, q)
            WRITE(6,'(5e14.6)') pos, SQRT(SUM(vel*vel)), nrg
        END IF

        CALL rk4(pos, vel, dt, w, q)
    END DO

END PROGRAM

SUBROUTINE rk4(x, v, h, w, q)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: h, w, q
    DOUBLE PRECISION :: ht, hs, hh !, v
    DOUBLE PRECISION, DIMENSION(1:3), INTENT(inout) :: x, v
    DOUBLE PRECISION, DIMENSION(1:4, 1:3) :: kp = 0, kv = 0

    INTERFACE
        PURE FUNCTION dvdt(p, c, a, d)
            DOUBLE PRECISION, DIMENSION(1:3), INTENT(in) :: p, c
            DOUBLE PRECISION, DIMENSION(1:3) :: dvdt
            DOUBLE PRECISION, INTENT(in) :: a, d
        END FUNCTION dvdt
    END INTERFACE

    ht = h/3d+0; hs = h/6d+0; hh = h/2d+0    !third, sixth and half of h

    kp(1,:) = v                 !runge-kutta 4th order, change in pos, then vel
    kv(1,:) = dvdt(x, v, q, w)

    kp(2,:) = v + hh*kv(1,:)
    kv(2,:) = dvdt((x + hh*kp(1,:)), kp(2,:), q, w)

    kp(3,:) = v + hh*kv(2,:)
    kv(3,:) = dvdt((x + hh*kp(2,:)), kp(3,:), q, w)

    kp(4,:) = v + h*kv(3,:)
    kv(4,:) = dvdt((x + h*kp(3,:)), kp(4,:), q, w)

!   ____Assign h*k1-4 to velocity and position_________________________________
    v = v + hs*(kv(1,:) + kv(4,:)) + ht*(kv(2,:) + kv(3,:))
    x = x + hs*(kp(1,:) + kp(4,:)) + ht*(kp(2,:) + kp(3,:))

END SUBROUTINE rk4

!   ____Acceleration___________________________________________________________
PURE FUNCTION dvdt(p, c, a, d)
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(1:3), INTENT(in) :: p, c
    DOUBLE PRECISION, DIMENSION(1:3) :: dvdt
    DOUBLE PRECISION, INTENT(in) :: a, d
    DOUBLE PRECISION :: r

    r = SQRT(SUM(p*p))

    dvdt = -(p/r)*((a*a + r*r)**(-1.5))   !plummer

    dvdt(1) = dvdt(1) + d*2d+0*c(2)     !coriolis x
    dvdt(2) = dvdt(2) - d*2d+0*c(1)     !coriolis y

    dvdt(1) = dvdt(1) + 3d+0*p(1)*d*d   !centrifugal

    dvdt(3) = dvdt(3) - p(3)*d*d        !tidal

!    dvdt = dvdt*(1d+0 - 1d+0/(1d+3*r + 1))
END FUNCTION dvdt

!   ____Energy_________________________________________________________________
DOUBLE PRECISION PURE FUNCTION energy(x, v, w, a)
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(1:3), INTENT(in) :: x, v
    DOUBLE PRECISION, INTENT(in) :: w, a
    DOUBLE PRECISION :: r, s

    r = SUM(x*x) !both squares since they're only used once as squares, so...
    s = SUM(v*v)

    energy = s/2d+0 - (a*a + r)**(-0.5) + .5d+0*w*w*(x(3)**2 - 3d+0*x(1)**2)

END FUNCTION energy
