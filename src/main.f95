program molecular_dynamics_3d

use omp_lib

use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

integer N,i,j,p,iargc,numb
real (kind = 8) Lx,Ly,Lz,Vxmax,Vymax,Vzmax,tmin,tmax,dt,svx,svy,svz,fx,fy,fz,scl,KE,Temp,t,xdiff,ydiff,zdiff,r,k,tau
real (kind = 8), dimension (:), allocatable :: x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az
integer, parameter :: seed = 99999999
character (len=90) :: filename
character (len=32) :: arg
call srand(seed)

numb = iargc()

call getarg(1, arg)

if (numb == 0) then
  arg = "input.ini"
  write(*,*) "No input file provided. Default input.ini is being used."
endif

cfg = parse_cfg(arg)

call cfg%get("particle","N",N)

allocate(x(N))
allocate(y(N))
allocate(z(N))
allocate(vx(N))
allocate(vy(N))
allocate(vz(N))
allocate(ux(N))
allocate(uy(N))
allocate(uz(N))
allocate(ax(N))
allocate(ay(N))
allocate(az(N))

call cfg%get("temperature","Temp",Temp)

call cfg%get("screening","k",k)

write ( *, '(a,i8)' ) &
    '  The number of processors available = ', omp_get_num_procs ( )
write ( *, '(a,i8)' ) &
    '  The number of threads available    = ', omp_get_max_threads ( )

call cfg%get("length","Lx",Lx)
call cfg%get("length","Ly",Ly)
call cfg%get("length","Lz",Lz)

call cfg%get("time","tmin",tmin)
call cfg%get("time","tmax",tmax)
call cfg%get("time","dt",dt)

Vxmax = 1.0d0
Vymax = 1.0d0
Vzmax = 1.0d0

svx = 0.0d0
svy = 0.0d0
svz = 0.0d0

!=============== Initial Condition ==================

do i = 1,N
  x(i) = (rand())*2.0d0*Lx - Lx
  y(i) = (rand())*2.0d0*Ly - Ly
  z(i) = (rand())*2.0d0*Lz - Lz
  vx(i) = (rand())*Vxmax - Vxmax/2.0d0
  vy(i) = (rand())*Vymax - Vymax/2.0d0
  vz(i) = (rand())*Vzmax - Vzmax/2.0d0
  svx = svx + vx(i)
  svy = svy + vy(i)
  svz = svz + vz(i)
enddo

do i = 1,N
  vx(i) = vx(i) - svx/dfloat(N)
  vy(i) = vy(i) - svy/dfloat(N)
  vz(i) = vz(i) - svz/dfloat(N)
enddo

do i = 1,N
  ax(i) = 0.0d0
  ay(i) = 0.0d0
  az(i) = 0.0d0
    do j = 1,N
      if (i .ne. j) then
        xdiff = ( x(i)-x(j) ) !- nint ((x(i)-x(j))/(2.0d0*Lx)) * 2.0d0*Lx
        ydiff = ( y(i)-y(j) ) !- nint ((y(i)-y(j))/(2.0d0*Ly)) * 2.0d0*Ly
        zdiff = ( z(i)-z(j) ) !- nint ((z(i)-z(j))/(2.0d0*Lz)) * 2.0d0*Lz
        r = dsqrt (xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
        fx = xdiff*(1.0d0+k*r)*dexp(-k*r)/(r*r*r)
        fy = ydiff*(1.0d0+k*r)*dexp(-k*r)/(r*r*r)
        fz = zdiff*(1.0d0+k*r)*dexp(-k*r)/(r*r*r)
        ax(i) = ax(i) + fx
        ay(i) = ay(i) + fy
        az(i) = az(i) + fz
      endif
    enddo
enddo

!=============== Time Loop =======================


do t = tmin+dt, tmax, dt

KE = 0.0d0
  !$OMP PARALLEL SHARED(Lx,Ly,Lz,x,y,z,ux,uy,uz,vx,vy,vz,N,dt,ax,ay,az,k) PRIVATE (i,j,fx,fy,fz,xdiff,ydiff,zdiff,r)
  !$OMP DO
  do i = 1,N
    ux(i) = vx(i) + ax(i) * dt / 2.0d0
    uy(i) = vy(i) + ay(i) * dt / 2.0d0
    uz(i) = vz(i) + az(i) * dt / 2.0d0
    x(i) = x(i) + ux(i) * dt
    y(i) = y(i) + uy(i) * dt
    z(i) = z(i) + uz(i) * dt
    !x(i) = x(i) - (int(x(i)/(Lx)))  * 2.0d0 * Lx      ! Periodic Boundary Condition
    !y(i) = y(i) - (int(y(i)/(Ly))) * 2.0d0 * Ly      ! Periodic Boundary Condition
    !z(i) = z(i) - (int(z(i)/(Lz)))  * 2.0d0 * Lz      ! Periodic Boundary Condition
    !Perfectly Reflecting Boundary Conditions...
    if (x(i) .gt. Lx .or. x(i) .lt. -Lx) then
      x(i) = x(i) - dt*ux(i)
      ux(i) = -ux(i)
    endif

    if (y(i) .gt. Ly .or. y(i) .lt. -Ly) then
      y(i) = y(i) - dt*uy(i)
      uy(i) = -uy(i)
    endif

    if (z(i) .gt. Lz .or. z(i) .lt. -Lz) then
      z(i) = z(i) - dt*uz(i)
      uz(i) = -uz(i)
    endif
  enddo
  !$OMP END DO

  !$OMP DO
  do i = 1,N
    ax(i) = 0.0d0
    ay(i) = 0.0d0
    az(i) = 0.0d0
      do j = 1,N
        if (i .ne. j) then
          xdiff = ( x(i)-x(j) ) !- nint ((x(i)-x(j))/(2.0d0*Lx)) * 2.0d0*Lx
          ydiff = ( y(i)-y(j) ) !- nint ((y(i)-y(j))/(2.0d0*Ly)) * 2.0d0*Ly
          zdiff = ( z(i)-z(j) ) !- nint ((z(i)-z(j))/(2.0d0*Lz)) * 2.0d0*Lz
          r = dsqrt (xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
          fx = xdiff*(1.0d0+k*r)*dexp(-k*r)/(r*r*r)
          fy = ydiff*(1.0d0+k*r)*dexp(-k*r)/(r*r*r)
          fz = zdiff*(1.0d0+k*r)*dexp(-k*r)/(r*r*r)
          ax(i) = ax(i) + fx
          ay(i) = ay(i) + fy
          az(i) = az(i) + fz
        endif
      enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  !$OMP PARALLEL SHARED(p,x,y,z,dt,N,KE) PRIVATE (i,vx,vy,vz)
  !$OMP DO
  do i = 1,N
    vx(i) = ux(i) + ax(i) * dt / 2.0d0
    vy(i) = uy(i) + ay(i) * dt / 2.0d0
    vz(i) = uz(i) + az(i) * dt / 2.0d0
    KE = KE + ( vx(i)*vx(i) + vy(i)*vy(i)  + vz(i)*vz(i)) / 2.0d0
  enddo
  !$OMP END DO


 !$OMP DO
  do i = 1,N
    p = int(t/dt)
    write(filename, '("output/fort.",I8.8)') p+100
    open(unit=p+100,file=filename,status='unknown')
    write(p+100,*) x(i),y(i),z(i)
    ! write (p+100,*) x(i), y(i), z(i)
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
      vz(i) = scl * vz(i)
    enddo
  else
    do i = 1,N
      vx(i) = vx(i)
      vy(i) = vy(i)
      vz(i) = vz(i)
    enddo
  endif

enddo     !time

end program molecular_dynamics_3d
