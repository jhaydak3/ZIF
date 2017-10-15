subroutine addCountxy(Count,NGrid,thisx,thisy,XX,YY,r_cutoff)
	implicit none
	integer, intent(in) :: NGrid
!f2p intent(in) ::  NGrid
	real(8), intent(inout), dimension(0:NGrid-1, 0:NGrid-1) :: Count, XX, YY
!f2p intent(in,out) :: Count, XX, YY
	real(8), intent(in) :: thisx, thisy, r_cutoff
!f2p intent(in) :: thisx, thisy, r_cutoff
	integer :: i,j,xndx(1),yndx(1)
	integer :: xlowndx, ylowndx, xhighndx, yhighndx
	real(8) :: dist, dx, dy
	real(8), dimension(0:NGrid-1) :: xvec, yvec
	
	dx = XX(0,1) - XX(0,0)  !assumes linspace used to generate spacing
	dy = YY(1,0) - YY(0,0)
	xvec = XX(0,:)
	yvec = YY(:,0)
	xndx = minloc(abs(xvec-thisx))
	yndx = minloc(abs(yvec-thisy))
	!write(*,*) xndx,yndx,thisx,thisy
	xlowndx = nint(xndx(1) - (r_cutoff / dx * 2.5))
	xhighndx = nint(xndx(1) + (r_cutoff / dx * 2.5))
	ylowndx = nint(yndx(1) - (r_cutoff / dy * 2.5))
	yhighndx = nint(yndx(1) + (r_cutoff / dy * 2.5))
	if (xlowndx < 0) xlowndx = 0
	if (ylowndx < 0) ylowndx = 0
	if (xhighndx > NGrid-1) xhighndx = NGrid - 1
	if (yhighndx > NGrid-1) yhighndx = Ngrid - 1
	do j = xlowndx, xhighndx
		do i = ylowndx, yhighndx
			dist = (thisx-XX(i,j))**2 + (thisy-YY(i,j))**2
			dist = sqrt(dist)
			if (dist .LT. r_cutoff) then
				Count(i,j) = Count(i,j) + 1
			end if
			!write(*,*) thisx,thisy, XX(i,j), YY(i,j), dist
		end do
	end do

end subroutine

subroutine addCountxz(Count,NGrid,thisx,thisz,XX,ZZ,r_cutoff)
	implicit none
	integer, intent(in) :: NGrid
!f2p intent(in) ::  NGrid
	real(8), intent(inout), dimension(0:NGrid-1, 0:NGrid-1) :: Count, XX, ZZ
!f2p intent(in,out) :: Count, XX, ZZ
	real(8), intent(in) :: thisx, thisz, r_cutoff
!f2p intent(in) :: thisx, thisz, r_cutoff
	integer :: i,j,xndx(1),zndx(1)
	integer :: xlowndx, zlowndx, xhighndx, zhighndx
	real(8) :: dist, dx, dz
	real(8), dimension(0:NGrid-1) :: xvec, zvec
	
	dx = XX(0,1) - XX(0,0)  !assumes linspace used to generate spacing
	dz = ZZ(1,0) - ZZ(0,0)
	xvec = XX(0,:)
	zvec = ZZ(:,0)
	xndx = minloc(abs(xvec-thisx))
	zndx = minloc(abs(zvec-thisz))
	!write(*,*) xndx,yndx,thisx,thisy
	xlowndx = nint(xndx(1) - (r_cutoff / dx * 2.5))
	xhighndx = nint(xndx(1) + (r_cutoff / dx * 2.5))
	zlowndx = nint(zndx(1) - (r_cutoff / dz * 2.5))
	zhighndx = nint(zndx(1) + (r_cutoff / dz * 2.5))
	if (xlowndx < 0) xlowndx = 0
	if (zlowndx < 0) zlowndx = 0
	if (xhighndx > NGrid-1) xhighndx = NGrid - 1
	if (zhighndx > NGrid-1) zhighndx = Ngrid - 1
	do j = xlowndx, xhighndx
		do i = zlowndx, zhighndx
			dist = (thisx-XX(i,j))**2 + (thisz-ZZ(i,j))**2
			dist = sqrt(dist)
			if (dist .LT. r_cutoff) then
				Count(i,j) = Count(i,j) + 1
			end if
			!write(*,*) thisx,thisy, XX(i,j), YY(i,j), dist
		end do
	end do

end subroutine

subroutine addCountyz(Count,NGrid,thisy,thisz,YY,ZZ,r_cutoff)
	implicit none
	integer, intent(in) :: NGrid
!f2p intent(in) ::  NGrid
	real(8), intent(inout), dimension(0:NGrid-1, 0:NGrid-1) :: Count, YY, ZZ
!f2p intent(in,out) :: Count, YY, ZZ
	real(8), intent(in) :: thisy, thisz, r_cutoff
!f2p intent(in) :: thisy, thisz, r_cutoff
	integer :: i,j,yndx(1),zndx(1)
	integer :: ylowndx, zlowndx, yhighndx, zhighndx
	real(8) :: dist, dy, dz
	real(8), dimension(0:NGrid-1) :: yvec, zvec
	
	dy = YY(0,1) - YY(0,0)  !assumes linspace used to generate spacing
	dz = ZZ(1,0) - ZZ(0,0)
	yvec = YY(0,:)
	zvec = ZZ(:,0)
	yndx = minloc(abs(yvec-thisy))
	zndx = minloc(abs(zvec-thisz))
	!write(*,*) yndx,zndx,thisy,thisz
	ylowndx = nint(yndx(1) - (r_cutoff / dy * 2.5))
	yhighndx = nint(yndx(1) + (r_cutoff / dy * 2.5))
	zlowndx = nint(zndx(1) - (r_cutoff / dz * 2.5))
	zhighndx = nint(zndx(1) + (r_cutoff / dz * 2.5))
	if (ylowndx < 0) ylowndx = 0
	if (zlowndx < 0) zlowndx = 0
	if (yhighndx > NGrid-1) yhighndx = NGrid - 1
	if (zhighndx > NGrid-1) zhighndx = Ngrid - 1
	do j = ylowndx, yhighndx
		do i = zlowndx, zhighndx
			dist = (thisy-YY(i,j))**2 + (thisz-ZZ(i,j))**2
			dist = sqrt(dist)
			if (dist .LT. r_cutoff) then
				Count(i,j) = Count(i,j) + 1
			end if
			!write(*,*) thisy,thisz, YY(i,j), ZZ(i,j), dist
		end do
	end do

end subroutine