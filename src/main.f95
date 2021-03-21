program molecular_dynamics_2d
use omp_lib

implicit none
integer N,i,j,p
real*8 Lx,Ly,Vxmax,Vymax,tmax,dt,svx,svy,fx,fy,scl,KE,Temp,t,xdiff,ydiff,r,tau
real, dimension (100) :: x,y,vx,vy,ux,uy,ax,ay
integer, parameter :: seed = 99999999
call srand(seed)

Temp = 0.010d0

write ( *, '(a,i8)' ) &
    '  The number of processors available = ', omp_get_num_procs ( )
write ( *, '(a,i8)' ) &
    '  The number of threads available    = ', omp_get_max_threads ( )


N = 100

Lx = 1.0d0
Ly = 1.0d0

Vxmax = 1.0d0
Vymax = 1.0d0

tmax = 1.0d0
dt = 0.010d0

svx = 0.0d0
svy = 0.0d0

!=============== Initial Condition ==================

do i = 1,N
  x(i) = (rand())*2.0d0*Lx - Lx
  y(i) = (rand())*2.0d0*Ly - Ly
  vx(i) = (rand())*Vxmax - Vxmax/2.0d0
  vy(i) = (rand())*Vymax - Vymax/2.0d0
  svx = svx + vx(i)
  svy = svy + vy(i)
enddo

do i = 1,N
  vx(i) = vx(i) - svx/dfloat(N)
  vy(i) = vy(i) - svy/dfloat(N)
enddo

do i = 1,N
  ax(i) = 0.0d0
  ay(i) = 0.0d0
    do j = 1,N
      if (i .ne. j) then
        xdiff = ( x(i)-x(j) ) - nint ((x(i)-x(j))/(2.0d0*Lx)) * 2.0d0*Lx
        ydiff = ( y(i)-y(j) ) - nint ((y(i)-y(j))/(2.0d0*Ly)) * 2.0d0*Ly
        r = dsqrt (xdiff*xdiff + ydiff*ydiff)
        fx = xdiff/(r*r*r)
        fy = ydiff/(r*r*r)
        ax(i) = ax(i) + fx
        ay(i) = ay(i) + fy
      endif
    enddo
enddo

!=============== Time Loop =======================


do t = 0.0d0+dt, tmax, dt

KE = 0.0d0
  !$OMP PARALLEL SHARED(Lx,Ly,x,y,ux,uy,vx,vy,N,dt,ax,ay) PRIVATE (i,j,fx,fy,xdiff,ydiff,r)
  !$OMP DO
  do i = 1,N
    ux(i) = vx(i) + ax(i) * dt / 2.0d0
    uy(i) = vy(i) + ay(i) * dt / 2.0d0
    x(i) = x(i) + ux(i) * dt
    y(i) = y(i) + uy(i) * dt
    x(i) = x(i) - (int(x(i)/(Lx)))  * 2.0d0 * Lx      ! Periodic Boundary Condition
    y(i) = y(i) - (int(y(i)/(Ly))) * 2.0d0 * Ly      ! Periodic Boundary Condition
  enddo
  !$OMP END DO

  !$OMP DO
  do i = 1,N
    ax(i) = 0.0d0
    ay(i) = 0.0d0
      do j = 1,N
        if (i .ne. j) then
          xdiff = ( x(i)-x(j) ) - nint ((x(i)-x(j))/(2.0d0*Lx)) * 2.0d0*Lx
          ydiff = ( y(i)-y(j) ) - nint ((y(i)-y(j))/(2.0d0*Ly)) * 2.0d0*Ly
          r = dsqrt (xdiff*xdiff + ydiff*ydiff)
          fx = xdiff/(r*r*r)
          fy = ydiff/(r*r*r)
          ax(i) = ax(i) + fx
          ay(i) = ay(i) + fy
        endif
      enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  !$OMP PARALLEL SHARED(p,x,y,dt,N,KE) PRIVATE (i,vx,vy)
  !$OMP DO
  do i = 1,N
    vx(i) = ux(i) + ax(i) * dt / 2.0d0
    vy(i) = uy(i) + ay(i) * dt / 2.0d0
    KE = KE + ( vx(i)*vx(i) + vy(i)*vy(i) ) / 2.0d0
  enddo
  !$OMP END DO


 !$OMP DO
  do i = 1,N
    p = int(t/dt)
    write (p+100,*) x(i), y(i)
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  tau = 10.0d0 * dt
  scl = dsqrt (1.0d0 + (dt/tau) * ((Temp/(2.0d0*KE/(3.0d0*dfloat(N)) )) -1.0d0))

  if (t .le. tmax/2.0d0) then
    do i = 1,N
      vx(i) = scl * vx(i)
      vy(i) = scl * vy(i)
    enddo
  else
    vx(i) = vx(i)
    vy(i) = vy(i)
  endif

enddo     !time


end program molecular_dynamics_2d
