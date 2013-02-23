!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

! Real space interpolation by B-splines

!-----------------------------------------------------------------------
SUBROUTINE bspline_interpolation (nptx, rg, rhor, rhoint)
  !---------------------------------------------------------------------
  !
  ! Use B-spline interpolation instead of Fourier interpolation
  !
  USE kinds,     ONLY : dp
  USE io_global, ONLY : stdout, ionode
  USE fft_base,  ONLY : dfftp
  USE cell_base, ONLY : bg
  USE bspline  
  !---------------------------------------------------------------------
  implicit none
  integer, intent(in) :: nptx
  real(dp), intent(in) :: rg(3,nptx) ! in alat units
  real(dp), intent(in) :: rhor(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x)
  real(dp), intent(out) :: rhoint(nptx)
  !---------------------------------------------------------------------
  real(dp), allocatable :: xv(:), yv(:), zv(:)
  real(dp), allocatable :: xknot(:), yknot(:), zknot(:)
  real(dp), allocatable :: bcoef(:)
  real(dp), allocatable :: rhoext(:,:,:)
  integer, parameter :: kx = 5, ky = 5, kz = 5   ! order of B-spline
  integer :: nx, ny, nz, i, j, k, ii, jj, kk, ierr

  write(stdout,'(5X,"Interpolation by B-splines")') 

  nx = dfftp%nr1
  ny = dfftp%nr2
  nz = dfftp%nr3

  ! extend grid in all directions
  allocate(rhoext(-kx+1:nx+kx,-ky+1:ny+ky,-kx+1:nz+kz))
  do i = -kx+1, nx+kx
     ii = i
     if (i <= 0) ii = i+nx
     if (i > nx) ii = i-nx
     do j = -ky+1, ny+ky
        jj = j
        if (j <= 0) jj = j+ny
        if (j > ny) jj = j-ny
        do k = -kz+1, nz+kz
           kk = k
           if (k <= 0) kk = k+nz
           if (k > nz) kk = k-nz
           rhoext(i,j,k) = rhor(ii, jj, kk)
        enddo
     enddo
  enddo
  nx = nx + 2*kx
  ny = ny + 2*ky
  nz = nz + 2*kz

  ! prepare B-spline interpolation
  allocate (xv(nx), yv(ny), zv(nz) )
  allocate (xknot(nx+kx), yknot(ny+ky), zknot(nz+kz) )
  allocate (bcoef(nx*ny*nz))
 
 ! setup uniform grid along x
  do i = 1, nx
     xv(i) = dble(i-kx-1)/dble(nx-2*kx)
  enddo
  call dbsnak(nx, xv, kx, xknot, ierr)
  if (ierr /= 0) call errore('bspline_interpolation', 'error in dbsnak/x', ierr)

  ! setup uniform grid along y
  do i = 1, ny
     yv(i) = dble(i-ky-1)/dble(ny-2*ky)  
  enddo
  call dbsnak(ny, yv, ky, yknot, ierr)
  if (ierr /= 0) call errore('bspline_interpolation', 'error in dbsnak/y', ierr)

  ! setup uniform grid along z
  do i = 1, nz
     zv(i) = dble(i-kz-1)/dble(nz-2*kz)
  enddo
  call dbsnak(nz, zv, kz, zknot, ierr)
  if (ierr /= 0) call errore('bspline_interpolation', 'error in dbsnak/z', ierr)

  ! setup B-spline coefficients
  call dbs3in(nx,xv,ny,yv,nz,zv,rhoext,nx,ny,kx,ky,kz,xknot,yknot,zknot,bcoef,ierr)
  if (ierr /= 0) call errore('bspline_interpolation', 'error in dbs3in', ierr)

  ! transform grid points in crystal coordinates
  ! TODO: refold positions
  call cryst_to_cart(nptx, rg, bg, -1)

  ! interpolate
  do i = 1, nptx
     rhoint(i) = dbs3vl(rg(1,i),rg(2,i),rg(3,i),kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef,ierr)
     if (ierr /= 0) call errore('bspline_interpolation', 'error in dbs3vl', ierr)
  enddo

  ! we print the charge on output
  write(stdout, '(5x,"Min, Max charge: ",2f12.6)') maxval(rhoint), minval(rhoint)

END SUBROUTINE bspline_interpolation


!-----------------------------------------------------------------------
SUBROUTINE plot_1d_bspline (nptx, m1, x0, e, rhor, alat, iflag, ounit)
  !---------------------------------------------------------------------
  !
  ! Use B-spline interpolation instead of Fourier
  !
  USE kinds,     ONLY : dp
  USE io_global, ONLY : stdout, ionode
  USE fft_base,  ONLY : dfftp
  !---------------------------------------------------------------------
  implicit none
  integer, intent(in) :: nptx, iflag, ounit
  real(dp), intent(in) :: e(3), x0(3), m1, alat
  real(dp), intent(in) :: rhor(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x)
  !---------------------------------------------------------------------
  real(dp), allocatable :: rg(:,:), carica(:)
  real(dp) :: deltax
  integer :: i

  if (iflag == 0) &
    call errore('plot_1d_bsplint', 'spherical average incompatible with B-splines', 1)

  ! grid in cartesian coordinates, in units of alat
  allocate( rg(3,nptx), carica(nptx) )
  deltax = dble(m1) / dble(nptx - 1) 
  do i = 1, nptx
     rg(1,i) = x0(1) + (i-1) * deltax*e(1)
     rg(2,i) = x0(2) + (i-1) * deltax*e(2)
     rg(3,i) = x0(3) + (i-1) * deltax*e(3)
  enddo

  ! interpolate
  call bspline_interpolation(nptx, rg, rhor, carica) 

  ! we print the charge on output
  if (ionode) then
     do i = 1, nptx
        write (ounit, '(2f20.10)') deltax*dble(i-1), carica(i)
     enddo
  endif

END SUBROUTINE plot_1d_bspline


!-----------------------------------------------------------------------
SUBROUTINE plot_2d_bspline (nx, ny, m1, m2, x0, e1, e2, rhor, alat, &
     at, nat, tau, atm, ityp, output_format, ounit)
  !-----------------------------------------------------------------------
  !
  ! Use B-spline interpolation instead of Fourier
  !
  USE kinds,     ONLY : dp
  USE io_global, ONLY : stdout, ionode
  USE fft_base,  ONLY : dfftp
  !---------------------------------------------------------------------
  implicit none
  integer, intent(in) :: nx, ny, nat, ityp (nat), output_format, ounit
  real(dp), intent(in) :: e1(3), e2(3), x0(3), m1, m2, alat, tau(3,nat), at(3,3)
  character(len=3), intent(in) :: atm(*)
  real(dp), intent(in) :: rhor(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x)
  !---------------------------------------------------------------------
  real(dp), allocatable :: rg(:,:,:), carica(:,:)
  real(dp) :: deltax, deltay
  integer :: i, j, nptx

  ! grid in cartesian coordinates, in units of alat
  allocate( rg(3,nx,ny), carica(nx,ny) )
  deltax = dble(m1) / dble(nx - 1) 
  deltay = dble(m2) / dble(ny - 1) 
  do i = 1, nx
     do j = 1, ny
        rg(:,i,j) = x0(:) + (i-1)*deltax*e1(:) + (j-1)*deltay*e2(:)
     enddo
  enddo

  ! interpolate
  nptx = nx*ny
  call bspline_interpolation(nptx, rg(1,1,1), rhor, carica(1,1)) 

  ! and we print the charge on output
  if (ionode) then
     if (output_format == 0) then
        !
        !     gnuplot format
        !
        !         write(ounit,'(2i6)') nx,ny
        do i = 1, nx
           write (ounit, '(e25.14)') (  dble(carica(i,j)), j = 1, ny )
           write (ounit, * )
        enddo
     elseif (output_format == 1) then
        !
        !     contour.x format
        !
        write (ounit, '(3i5,2e25.14)') nx, ny, 1, deltax, deltay
        write (ounit, '(4e25.14)') ( (  dble(carica(i,j)), j = 1, ny ), i = 1, nx )
     elseif (output_format == 2) then
        !
        !     plotrho format
        !
        write (ounit, '(2i4)') nx - 1, ny - 1
        write (ounit, '(8f8.4)') (deltax * (i - 1) , i = 1, nx)
        write (ounit, '(8f8.4)') (deltay * (j - 1) , j = 1, ny)
        write (ounit, '(6e12.4)') ( (  dble(carica(i,j)), i = 1, nx ), j = 1, ny )
        write (ounit, '(3f8.4)') x0
        write (ounit, '(3f8.4)') (m1 * e1 (i) , i = 1, 3)
        write (ounit, '(3f8.4)') (m2 * e2 (i) , i = 1, 3)

     elseif (output_format == 3) then
        !
        ! xcrysden's xsf format
        !
        call xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
        call xsf_datagrid_2d (carica, nx, ny, m1, m2, x0, e1, e2, alat, ounit)
     elseif (output_format == 7) then
        !
        !     gnuplot format : x, y, f(x,y)
        !
        do i=1, nx
           do j=1, ny 
              write (ounit, '(3e20.8)')  alat*deltax * (i - 1), &
                      alat*deltay * (j - 1), dble(carica(i,j))
           enddo
           write(ounit, *)
        enddo
     else
        call errore('plot_2d', 'wrong output_format', 1)
     endif
  endif

END SUBROUTINE plot_2d_bspline






#if 0
!-----------------------------------------------------------------------
SUBROUTINE extend_grid(nx, ny, nz, rhor, kx, ky, kz, rhoout)
  !---------------------------------------------------------------------
  !
  ! Extend the grid
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE fft_base,  ONLY : dfftp
  !---------------------------------------------------------------------
  implicit none
  integer, intent(in) :: nx, ny, nz   ! grid size
  integer, intent(in) :: kx, ky, kz   ! B-spline order
  real(dp), intent(in) :: rhor(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x)
  real(dp), intent(out) :: rhoout(-kx+1:nx+kx,-ky+1:ny+ky,-kx+1:nz+kz)
  !---------------------------------------------------------------------
  integer :: i, j, k, ii, jj, kk

  do i = -kx+1, nx+kx
     ii = i
     if (i <= 0) ii = i+nx
     if (i > nx) ii = i-nx
     do j = -ky+1, ny+ky
        jj = j
        if (j <= 0) jj = j+ny
        if (j > ny) jj = j-ny
        do k = -kz+1, nz+kz
           kk = k
           if (k <= 0) kk = k+nz
           if (k > nz) kk = k-nz
           rhoout(i,j,k) = rhor(ii, jj, kk)
        enddo
     enddo
  enddo

END SUBROUTINE extend_grid


!-----------------------------------------------------------------------
SUBROUTINE prepare_bspline(nx, ny, nz, rho, kx, ky, kz, xv, yv, zv, xknot, yknot, zknot, bcoef)
  !---------------------------------------------------------------------
  !
  ! Preapare B-spline interpolation: call dbsnak, then dbs3in
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE bspline  
  !---------------------------------------------------------------------
  implicit none
  integer, intent(in) :: nx, ny, nz   ! grid size
  integer, intent(in) :: kx, ky, kz   ! B-spline order
  real(dp), intent(in) :: rho(nx,ny,nz)
  real(dp), intent(out) :: xv(nx), yv(ny), zv(nz)
  real(dp), intent(out) :: xknot(nx+kx), yknot(ny+ky), zknot(nz+kz)
  real(dp), intent(out) :: bcoef(nx*ny*nz)
  !---------------------------------------------------------------------
  integer :: i, ierr
 
  ! setup uniform grid along x
  do i = 1, nx
     xv(i) = dble(i-kx-1)/dble(nx-2*kx)
  enddo
  call dbsnak(nx, xv, kx, xknot, ierr)
  if (ierr /= 0) call errore('prepare_bspline', 'error in dbsnak/x', ierr)

  ! setup uniform grid along y
  do i = 1, ny
     yv(i) = dble(i-ky-1)/dble(ny-2*ky)  
  enddo
  call dbsnak(ny, yv, ky, yknot, ierr)
  if (ierr /= 0) call errore('prepare_bspline', 'error in dbsnak/y', ierr)

  ! setup uniform grid along z
  do i = 1, nz
     zv(i) = dble(i-kz-1)/dble(nz-2*kz)
  enddo
  call dbsnak(nz, zv, kz, zknot, ierr)
  if (ierr /= 0) call errore('prepare_bspline', 'error in dbsnak/z', ierr)

  ! setup B-spline coefficients
  call dbs3in(nx,xv,ny,yv,nz,zv,rho,nx,ny,kx,ky,kz,xknot,yknot,zknot,bcoef,ierr)
  if (ierr /= 0) call errore('prepare_bspline', 'error in dbs3in', ierr)

END SUBROUTINE prepare_bspline  
#endif

