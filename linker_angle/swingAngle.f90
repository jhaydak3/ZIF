subroutine makeAngleList(timePlot,fileName,refFileName,numStepstoRead,numAtoms,AnglesVec,ndxArray,numAngles,numLinkers)
	implicit none
	integer, intent(in) :: numAngles, numAtoms, numLinkers, numStepstoRead
!f2p intent(in) :: numAngles, numAtoms, numLinkers, numStepstoRead
	real(8), intent(inout), dimension (0:numAngles-1) :: AnglesVec
!f2p intent(in,out) :: AnglesVec
	real(8), intent (inout), dimension(0:numStepstoRead-1,0:numLinkers) :: timePlot
!f2p intent(in,out) :: timePlot
	integer, intent(in), dimension(0:numLinkers-1,0:2) :: ndxArray
!f2p intent(in) :: ndxArray
	character (LEN=*), intent(in) :: fileName, refFileName
!f2p intent(in) :: fileName, refFileName
	integer :: i,j, k, C1ndx, C2ndx, C3ndx
	real(8), dimension(0:numAtoms) :: xvec, yvec, zvec
	real(8), dimension(0:2) :: C1, C2, C3, C, v1, Cref, v1ref
	real(8), dimension(0:numLinkers-1,0:2) :: C1ref, C2ref, C3ref
	real(8) :: thisAngle
	character (LEN=10) :: junk, junk2
	external getAngle
	
	!read through reference file
	open(2,file = refFileName)
	read(2,*)
		read(2,*)
		do j = 0, numAtoms-1
			read(2,*) junk, xvec(j), yvec(j), zvec(j)
		end do
		
		do j = 0, numLinkers-1
			C1ndx = ndxArray(j,0)
			C1ref(j,0) = xvec(C1ndx)
			C1ref(j,1) = yvec(C1ndx)
			C1ref(j,2) = zvec(C1ndx)
			C2ndx = ndxArray(j,1)
			C2ref(j,0) = xvec(C2ndx)
			C2ref(j,1) = yvec(C2ndx)
			C2ref(j,2) = zvec(C2ndx)
			C3ndx = ndxArray(j,2)
			C3ref(j,0) = xvec(C3ndx)
			C3ref(j,1) = yvec(C3ndx)
			C3ref(j,2) = zvec(C3ndx)
		end do
	close(2)
	!read all coords
	k = 0 ! thetaVec
	open(1, file = fileName)
	do i = 0, numStepstoRead-1
		if (mod(i,500) == 0) write(*,*) i
		read(1,*) junk
		read(1,*) junk, junk2, timePlot(i,0)
		do j = 0, numAtoms-1
			read(1,*) junk, xvec(j), yvec(j), zvec(j)
		end do
		
		do j = 0, numLinkers-1
			C1ndx = ndxArray(j,0)
			C1(0) = xvec(C1ndx)
			C1(1) = yvec(C1ndx)
			C1(2) = zvec(C1ndx)
			C2ndx = ndxArray(j,1)
			C2(0) = xvec(C2ndx)
			C2(1) = yvec(C2ndx)
			C2(2) = zvec(C2ndx)
			C3ndx = ndxArray(j,2)
			C3(0) = xvec(C3ndx)
			C3(1) = yvec(C3ndx)
			C3(2) = zvec(C3ndx)
			
			!do calculations
			C = .5 * (C2 + C3)
			v1 = C1 - C
			
			Cref = .5 * (C2ref(j,:) + C3ref(j,:))
			v1ref = C1ref(j,:) - Cref
			!write(*,*) i,j,k
			!write(*,*) AnglesVec(k), i, j
			call getAngle(v1,v1ref,thisAngle)
			AnglesVec(k) = thisAngle
			timePlot(i,j+1) = thisAngle
			!write(*,*) k
			k = k + 1
		end do
	end do
	!write(*,*) AnglesVec
	close(1)
end subroutine

subroutine getAngle(a,b,myAngle)
	implicit none
	real(8), intent(in), dimension(0:2) :: a,b
!f2p intent(in) :: a,b
	real(8), intent(inout) :: myAngle
!f2p intent(in,out) :: myAngle
	real(8) :: dotProd, na, nb
	real(8) :: PI = 3.14159265359
	
	dotProd = a(0)*b(0) + a(1)*b(1) + a(2)*b(2) 
	na = SQRT(a(0)**2+a(1)**2+a(2)**2)
	nb = SQRT(b(0)**2+b(1)**2+b(2)**2)
	myAngle = dotProd / (na*nb)
	!write(*,*) myAngle, size(myAngle)
	myAngle = DACOS(myAngle)
	myAngle = myAngle * (180/PI)
	!write(*,*) myAngle
end subroutine