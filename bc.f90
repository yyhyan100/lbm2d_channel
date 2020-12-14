subroutine bc()
!call non_eq_extra()
call bouncebk()
end subroutine

!----------------------------------------------------------------------------
subroutine non_eq_extra 
	use vars
	integer i,j,k
	real g2(0:8)
! lower	
	j=1
	do i=2,ied-1
		call getFeqAt(i,j,rho(i,j+1),u(i,j),v(i,j))
		g2(:)=g(:)
		call getFeqAt(i,j,rho(i,j+1),u(i,j+1),v(i,j+1))
		f(:,i,j)=g2(:)+f(:,i,j+1)-g(:)
	enddo
! upper	
	j=jed
	do i=2,ied-1
		call getFeqAt(i,j,rho(i,j-1),u(i,j),v(i,j))
		g2(:)=g(:)
		call getFeqAt(i,j,rho(i,j-1),u(i,j-1),v(i,j-1))
		f(:,i,j)=g2(:)+f(:,i,j-1)-g(:)
	enddo
! left	
	i=1
	do j=2,jed-1
		call getFeqAt(i,j,rho(i+1,j),u(i,j),v(i,j))
		g2(:)=g(:)
		call getFeqAt(i,j,rho(i+1,j),u(i+1,j),v(i+1,j))
		f(:,i,j)=g2(:)+f(:,i+1,j)-g(:)
	enddo
! right
	i=ied
	do j=2,jed-1
		call getFeqAt(i,j,rho(i-1,j),u(i,j),v(i,j))
		g2(:)=g(:)
		call getFeqAt(i,j,rho(i-1,j),u(i-1,j),v(i-1,j))
		f(:,i,j)=g2(:)+f(:,i-1,j)-g(:)
	enddo
	
! corner points, bounce back -- test
	i=1;j=1
	f(5,i,j)=f(7,i+1,j+1)
	i=ied;j=1
	f(6,i,j)=f(8,i-1,j+1)
	
end subroutine
!----------------------------------------------------------------------------
subroutine bouncebk 
	use vars
	integer i,j,k
! lower	
	j=1
	do i=1,ied-1
		f(2,i,j)=f(4,i,j)
		f(5,i,j)=f(7,i,j)
		f(6,i,j)=f(8,i,j)
	enddo
! upper	
	j=jed
	do i=1,ied-1
		f(4,i,j)=f(2,i,j)
		f(7,i,j)=f(5,i,j)
		f(8,i,j)=f(6,i,j)
	enddo
! left	
!	i=1
!	do j=2,jed-1
!		f(1,i,j)=f(3,i,j)
!		f(5,i,j)=f(7,i,j)
!		f(8,i,j)=f(6,i,j)
!	enddo
! right
	i=ied
	do j=1,jed
		do k=0,Q
			f(k,i,j)=2*f(k,i-1,j)-f(k,i-2,j)
		enddo
		!f(1,i,j)=2*f(1,i-1,j)-f(1,i-2,j)
		!f(5,i,j)=2*f(5,i-1,j)-f(5,i-2,j)
		!f(8,i,j)=2*f(8,i-1,j)-f(8,i-2,j)
		v(i,j)=0.0
		u(i,j)=2*u(i-1,j)-u(i-2,j)
	enddo
	
! corner points, bounce back -- test
!	i=1;j=1
!	f(5,i,j)=f(7,i+1,j+1)
!	i=ied;j=1
!	f(6,i,j)=f(8,i-1,j+1)
	
end subroutine
