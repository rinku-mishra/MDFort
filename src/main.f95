program molecular_dynamics_2d

use omp_lib

use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

integer N,i,j,p
real (kind = 8) Lx,Ly,Vxmax,Vymax,tmin,tmax,dt,svx,svy,fx,fy,scl,KE,Temp,t,xdiff,ydiff,r,tau
real (kind = 8), dimension (:), allocatable :: x,y,vx,vy,ux,uy,ax,ay
integer, parameter :: seed = 99999999
character (len=90) :: filename
call srand(seed)

cfg = parse_cfg("input.ini")

call cfg%get("particle","N",N)

allocate(x(N))
allocate(y(N))
allocate(vx(N))
allocate(vy(N))
allocate(ux(N))
allocate(uy(N))
allocate(ax(N))
allocate(ay(N))

call cfg%get("temperature","Temp",Temp)

write ( *, '(a,i8)' ) &
    '  The number of processors available = ', omp_get_num_procs ( )
write ( *, '(a,i8)' ) &
    '  The number of threads available    = ', omp_get_max_threads ( )

call cfg%get("length","Lx",Lx)
call cfg%get("length","Ly",Ly)

call cfg%get("time","tmin",tmin)
call cfg%get("time","tmax",tmax)
call cfg%get("time","dt",dt)

Vxmax = 1.0d0
Vymax = 1.0d0

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


do t = tmin+dt, tmax, dt

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
    write(filename, '("output/fort.",I8.8)') p+100
    open(unit=p+100,file=filename,status='unknown')
    write(p+100,*) x(i),y(i)
    ! write (p+100,*) x(i), y(i)
  enddo
 !$OMP END DO
 !$OMP END PARALLEL
  close (p+100)
  tau = 10.0d0 * dt
  scl = dsqrt (1.0d0 + (dt/tau) * ((Temp/(2.0d0*KE/(3.0d0*dfloat(N)) )) -1.0d0))

  if (t .le. tmax/2.0d0) then
    do i = 1,N
      vx(i) = scl * vx(i)
      vy(i) = scl * vy(i)
    enddo
  else
    do i = 1,N
      vx(i) = vx(i)
      vy(i) = vy(i)
    enddo
  endif

enddo     !time

end program molecular_dynamics_2d
