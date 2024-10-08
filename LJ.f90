program repulsion
  use mod_param
  use mod_force 
  use mtmod

  implicit none
  integer*4::ip,jp,inp
  integer*4::iloop,nloop
  real*8::rip,xip,yip,rin,xin,yin
  real*8::rnum,delta0,maxf,tot_t
  real*8::rhotemp,Unp,dPr
  real*8::a_mu,t_mu,t_begin
!! 
  
  call read_input
  mot0 = Df*gm;  ! Single cell motility
  R0 = dble(Npr0)*a;  ! Radius of sphere.
  R1 = dble(Npr1)*a;  ! Radius of sphere (Cell+Matrix).
  rc = 3.0d0*a;   ! Cut-off for LJ, here it is being used to determine the density in the neighbourhood of cells.
  area_c = pi*rc*rc
  Nrs = int((4.0/3.0)*pi*(R1**3-R0**3)*rhos)
  !Npp = int(dble(Nrs)*phi*0.58)  !!!why 0.58
  Npp = int(dble(Nrs)*poro_phi)  !!! AG deleted 0.58
  mu0 = mu;
  mu = mu0;

!!
  write(cnpar,'(g8.0)') 3+Npy*ndata
  formp = '('//trim(adjustl(cnpar))//'es16.3E2)'
  call open_track
!!
  call allocate_arrays
!!
  call sgrnd(iseed)
  xx = -2.0d0; yy = 0.0d0; zz = 0.0d0;
  dtime = 0.0d0;
  gmm = gm;
  ! y dependent friction factor 
  mot = 0;  ! Inititalize the motility as 0
  Aa = Aa0; !Aamin;    ! Inititalize the attraction potential Aa0
  a_mu = delta*dble(Nt)/32.0d0;
  t_mu = delta*dble(Nt)/16.0d0;
  t_begin = delta*dble(Nt)/2.0d0;

  call initial_condition

!  stop
!
!! Time marching 
  rhotemp = rho0;
  delta0 = delta;
!    do it = 1,Nt;
!    tt = dble(it)*delta;
!    mu =  muf + (mu0-muf)/(1+exp((tt-(t_begin-t_mu))/a_mu));
!    write(12,*)mu
!    end do
   ! stop

  do it = 1,Nt;
    tt = dble(it)*delta;
 !   mu =  muf + (mu0-muf)/(1+exp((tt-(t_begin-t_mu))/a_mu));
!
!! Calculate the force acting on the particle. f(rij) = - [\partial U/\partial r]
    fx = 0.0d0;    fy = 0.0d0;    fz = 0.0d0;
    fxrs = 0.0d0;  fyrs = 0.0d0;  fzrs = 0.0d0;
    fxpp = 0.0d0;  fypp = 0.0d0;  fzpp = 0.0d0;
!!
! First for the substrate
!
  call cal_force_matrix
  call cal_force_passive
!!
! Second for the cells
  call cal_force_cell

!! Adaptive time step
   delta = delta0
   maxf = (max(maxval(fx),maxval(fy))*(delta/gm))/a;
   if(maxf.gt.0.1d0) then
     nloop = (ceiling(maxf/0.1d0))
     delta = delta0/dble(nloop)
     do iloop = 1,nloop
    fx = 0.0d0;    fy = 0.0d0;    fz = 0.0d0;
    fxrs = 0.0d0;  fyrs = 0.0d0;  fzrs = 0.0d0;
       call cal_force_matrix
       call cal_force_passive
       call cal_force_cell
       call evolve_matrix
       call evolve_cell
     enddo
  ! write(99,*)delta,delta0
     goto 75
   endif

!! Evolving the particles 
!
!! Evolve the matrix 
    call evolve_matrix
    call evolve_passive
!
!! Evolve the rest of the cells 
   ! Deterministic process
    call evolve_cell
   ! Wiener process
75    call evolve_cell_Weiner


!!-- Cell proliferation
    if(stress_dependent) call cell_proliferation
!!-- Uniform Cell proliferation
    if((uniform_prol).and.((mod(it,int(0.1d0*tauG/delta)).eq.0))) call cell_proliferation_uniform
!!-- Uniform Cell proliferation at the center
    if((uniform_prol_center).and.((mod(it,int(0.1d0*tauG/delta)).eq.0))) call cell_proliferation_uniform_center
!!

!!
  if(mod(it,250).eq.0)call write_track
  !call write_track
!!
  enddo

  call close_track

end program repulsion



