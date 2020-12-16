program main
use vars
integer kstep
call init()
tt=0.0
kstep=1
do while (tt<=t_end)
	call streaming()
	call bc()
	call getMacro()
	call getFeq()
	call collision()
	call output(kstep)
	kstep=kstep+1
	tt=tt+dt
enddo
call output(kstep-1)
call deallocateField()
print *, "done!"
end program
!----------------------------------------------------------------------------
subroutine getMacro() 
use vars
integer i,j,k

do i=1,ied
do j=1,jed
	rho(i,j)=0.0
	u(i,j)=0.0
	v(i,j)=0.0
	do k=0,Q
		rho(i,j)=rho(i,j)+f(k,i,j)
		u(i,j)=u(i,j)+ei(k,1)*f(k,i,j)
		v(i,j)=v(i,j)+ei(k,2)*f(k,i,j)
	enddo
	u(i,j)=u(i,j)/rho(i,j)
	v(i,j)=v(i,j)/rho(i,j)
enddo
enddo
i=1
do j=2,jed-1
		rho(i,j)=(f(0,i,j)+f(2,i,j)+f(4,i,j)+2*(f(3,i,j)+f(6,i,j)+f(7,i,j)))/(1-u(i,j))
		u(i,j)=0.1
		v(i,j)=0.0
enddo
end subroutine
