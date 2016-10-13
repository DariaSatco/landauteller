program LandauTeller

use landaumod
use intlib

implicit none

integer,parameter:: neq=21 !system dimension
real(8):: Ri,Pin,actin, t, t0, Yi(16), A0(neq), energy, eti
!Ri,Pin initial values of distance, momentum
! actin - initial value of action Hamiltonian variable
! phii - initial value of the action
! [-t0, t0] - the time interval we explore
! t - the current time variable
! Yi - initial values of derivatives
! A0 - the initial vector of the system of equations
! A - the vector developing in time
! q - the phase, is changing from 0 to 2pi
! energy - the whole energy of the system
! eti - transitional energy (the energy of free particle)

integer:: inst, finst
! inst - initial state of oscillator
! fins - final state of oscillator
! actin=0.5*(inst+finst+1)

integer:: n1,n2
! n1 - number of steps in Runge-Kutta method
! n2 - number of steps in integration of S-matrix element

real(8):: h
! h=2*t0/n1


real(8), dimension(4,4):: Jacobian
! Jacobian - Jacobian for variables P,R,q,I

real(8):: Prob, Prob1, Prob2, Prob21, Prob3, Prob31  !probability variables
!P=1/4*pi^2*(ReS^2+ImS^2)
real(8):: ReS, ImS, ReS1, ImS1, ReS2, ImS2, ReS3, ImS3, ReS4, ImS4, ReS5, ImS5, err1, err2
! ReS, ImS - the real and imaginary part of S-matrix element

real(8), dimension(:), allocatable:: qgen, divqgen !generalized phase vector
real(8), dimension(:), allocatable:: partS, partS1
! partS - so call the expression uder sin and cos
real(8), dimension(:,:), allocatable:: A0matrix
! write in A0matrix array the ode system resolution values


!variables for integration subroutine
real(8), dimension(:), allocatable:: phase !vector of integration interval points
real(8), dimension(:), allocatable:: intfuncRE, intfuncIM, intfuncRE1, intfuncIM1
real(8), dimension(:), allocatable:: intfuncRE2, intfuncIM2, intfuncRE3, intfuncIM3
real(8), dimension(:), allocatable:: intfuncRE4, intfuncIM4, intfuncRE5, intfuncIM5
!vectors of function values


integer:: i,k,l,j,r !counters
real(8):: det !determinant variable
real(8):: ham, aven, deltaimax
!ham - total energy value
!aven - average energy loss
!deltaimax - maximum change of vibrational energy
real(8), allocatable:: Radius(:), deltai(:)

!dlsode variables
integer:: itol, itask, istate, iopt, lrw,liw,mf
integer, allocatable:: iwork(:)
real(8):: rtol, atol(neq)
real(8),allocatable:: rwork(:)

!semiclassical second oder perturbation theory variables
real(8):: q, omega, dqdq, arg, const, qgenpth, constantq, actoscil
real(8):: tau
real(8):: C1,C2,c2s, sumc2s, c2sfactor, Aif, A2c, A2s, psi0, func, funcind, funcdepend, func1st

!****************************************
!PROGRAM

! describe initial values
inst=0
finst=2
energy=10.0D+00
actin=(inst+finst+1)*0.5D+00
!eti=energy-actin
!Pin=-sqrt(2*lambda2*eti)
Ri=600.0D+00
Yi=(/1.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00,0.0D+00,0.0D+00,0.0D+00, &
0.0D+00,1.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00/)
U=1.0D+00 !potential factor

t0=100.0D+00 !trial initial time

n1=400 !number of time steps
n2=800 !number of phase steps

open(unit=23, file='analytic-inverse.txt', status='replace')
open(unit=30, file='action.txt', status='replace')

open(unit=17, file='functioncheck.txt', status='replace')
! functioncheck.txt - file with phase dependence of functions
open(unit=28, file='analytic.txt', status='replace')
! analytic.txt - file with calucalations using analytical formulas
open(unit=10, file='prob02(test).txt', status='replace')
! prob01.txt - probability data

!do r=1,5

!energy=6.0D+00+(r-1)*1.0D+00
eti=energy-actin
Pin=-sqrt(2*lambda2*eti)


!want to find the correct initial time
t=-t0 !trial initial time

open(unit=19, file='test.txt', status='replace')
! test.txt - file with the time evolution of system variables

A0=(/Ri,Pin,0.0D+00+lambda2*Ri/Pin,actin,0.0D+00, &
1.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00,0.0D+00,0.0D+00,0.0D+00, &
0.0D+00,1.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00/)
!A0 - the initial values vector

!variables for subroutine dlsode (see the explanation of arguments in odelib)
itol=2
rtol=0.0D+00
    do i=1,4
    atol(i)=1.0D-10
    end do
    atol(5)=1.0D-14
    do i=6,neq
    atol(i)=1.0D-10
    end do
itask=1
iopt=0
lrw=22+9*neq+neq**2
liw=20+neq
mf=22

!cycle for calculating radius development in oder to choose a correct initial time
        h=2*t0/n1 !h -step size in time
        allocate(Radius(n1))

        do j=1,n1
        allocate(rwork(lrw), iwork(liw))
        istate=1
        call dlsode(derivfunc,neq,A0,t,t+h,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
             Radius(j)=A0(1)

        deallocate(rwork,iwork)
        end do
        
!Look for correct initial time to get the simmetry
k=0
t=0.0D+00

	do 
	k=k+1
		if ((Radius(k+1)-Radius(k)) .ge. 0.0D+00) exit !look for minimum value
	t=t+h !time <---> Radius(k)
	end do

deallocate(Radius)

!******
t0=t !t - the correct initial time
!******

!the cycle for checking the time dependence of different functions at fixed initial values
!use this cycle to find necessary values for perturbation formulas
t=-t0
h=2*t0/n1 !h -step size in time


A0=(/Ri,Pin,0.0D+00+lambda2*Ri/Pin,actin, 0.0D+00, &
1.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00,0.0D+00,0.0D+00,0.0D+00, &
0.0D+00,1.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00/)

!control the value of potential (need to be very small)
write(19, *), U*exp(-alpha*(A0(1)-sqrt(2*A0(4))*cos(A0(3))))
!initial total energy
ham=hamiltonian(A0)

write(19,'(f10.5 f15.7 f15.7 f15.7 f15.7 f15.7 f15.7 f15.7)') t, A0(1), A0(2), A0(3), A0(4), A0(5),&
actionfunc(A0), ham

        do
 		allocate(rwork(lrw), iwork(liw))
        istate=1
        call dlsode(derivfunc,neq,A0,t,t+h,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
        !calculate the total energy to control stability
        ham=hamiltonianred(A0)

    !make from vector a Jacobi matrix
    Jacobian=reshape(A0(6:21), (/4,4/))
    !calculate the determinant
    det=FindDet(Jacobian,4)
                
        		
	!write the information in test.txt
    ! the oder is time, radius, momentum, phase, oscil momentum, action, function under action integral,
    ! total energy, Jacobian
    write(19,'(f10.5 f15.7 f15.7 f15.7 f15.7 f15.7 f15.7 f15.7 f15.7)') t, A0(1), A0(2), A0(3), A0(4), A0(5),&
    actionfunc(A0), ham, det

        deallocate(rwork,iwork)

            if (t .ge. t0) exit

        end do


close(19)

!*****************************************************************
!semiclassical second oder perturbation formulas

omega=alpha*sqrt(2*eti/lambda2)
tau=2.0D+00/omega
arg=pi/omega

!----------
C1=2*pi*lambda2*sqrt(2*actin)/(alpha*sinh(arg))
C2=4*pi**2*lambda2*(eti-actin*arg/tanh(arg))/(omega*(sinh(arg)))**2
const=C2/C1
!----------

c2sfactor=4*pi*lambda2**2*actin/(alpha**2*eti*sinh(2*arg))
sumc2s=0.0D+00
    do k=1,400
        sumc2s=sumc2s+real(k)/(real(k**2)+1.0D+00/(2*omega**2))**2
    end do
print*, 'sumc2s=', sumc2s

c2s=c2sfactor*(-1.0D+00+sumc2s/omega**2)
!----------

A2c=-0.5*c2s
A2s=0.5*C2
psi0=A0(5)
!----------

!********************************************************************

!start to make calculations for different initial phase
!********************************************************************

allocate(phase(n2+1))
allocate(qgen(n2+1),divqgen(n2+1))
allocate(partS(n2+1), partS1(n2+1))
allocate(intfuncRE(n2+1), intfuncIM(n2+1) )
allocate(intfuncRE1(n2+1), intfuncIM1(n2+1) )
allocate(intfuncRE2(n2+1), intfuncIM2(n2+1) )
allocate(intfuncRE3(n2+1), intfuncIM3(n2+1) )
allocate(intfuncRE4(n2+1), intfuncIM4(n2+1) )
allocate(intfuncRE5(n2+1), intfuncIM5(n2+1) )
allocate(A0matrix(neq,n2+1))

!************
	do i=1,n2+1 !phase change cycle
	
	phase(i)=0.0D+00+2*pi*(i-1)/n2

!************

A0=(/Ri,Pin,phase(i)+lambda2*Ri/Pin,actin,0.0D+00, &
1.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00,0.0D+00,0.0D+00,0.0D+00, &
0.0D+00,1.0D+00,0.0D+00,0.0D+00,0.0D+00,0.0D+00,1.0D+00/)
!A0 - the initial values vector

!calculate the initial total energy of the system
ham=hamiltonian(A0)
print*, 'hamiltonian0=', ham

!ode resolution
!---------------------------------------------
t=-t0 !initial time
h=2*t0/n1 !h -step size in time

	    !calculate the final state of the system
        do l=1,n1

            allocate(rwork(lrw), iwork(liw))

        istate=1
        call dlsode(derivfunc,neq,A0,t,t+h,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)

           deallocate(rwork,iwork)

	    end do
!---------------------------------------------
	!make from vector a Jacobi matrix
	Jacobian=reshape(A0(6:21), (/4,4/))
	!calculate the determinant
	det=FindDet(Jacobian,4)

    !calculate the final total energy
    ham=hamiltonian(A0)

!print total energy and jacobian values to control the correctness of calculations
print*, 'hamiltonian final=', ham
print*, 'det=', det

!****************************
!calculate perturbation theory values

        !functions do not depend on ni,nf
        q=phase(i)
        dqdq=1.0D+00-const*sin(q)
        actoscil=actin+C1*sin(q)+C2+c2s*sin(2*q)

        funcind=A2s*sin(2*q)+A2c*cos(2*q)-C1*cos(q)+psi0+omega*tau*C2+0.5*c2s+C1

        qgenpth=q+2*const*cos(q)!+omega*tau

        !function depends on ni, nf
        funcdepend=qgenpth*(inst-finst)

        func=funcind+funcdepend

        func1st=psi0+C1*(1.0D+00-cos(q))+q*(inst-finst)

write(28,'(f15.7 f15.7 f15.7 f15.7 f15.7 f15.7)') phase(i), dqdq, func, func1st, qgenpth-q, actoscil

!calculate vectors to integrate

!perturbation theory 2nd oder
intfuncRE1(i)=cos(func)!*dqdq
intfuncIM1(i)=sin(func)!*dqdq
!perturbation theory 1st oder
intfuncRE4(i)=cos(func1st)
intfuncIM4(i)=sin(func1st)

!***************************************************
!change the places of initial states
! ni <---> nf

        !functions depend on ni, nf
        funcdepend=qgenpth*(finst-inst)

        func=funcind+funcdepend

        func1st=psi0+C1*(1.0D+00-cos(q))+q*(finst-inst)

write(23, '(f15.7 f15.7 f15.7 f15.7 f15.7)') phase(i), func, func1st

!calculate vectors to integrate

!perturbation theory 2nd oder
intfuncRE3(i)=cos(func)!*dqdq
intfuncIM3(i)=sin(func)!*dqdq
!perturbation theory 1st oder
intfuncRE5(i)=cos(func1st)
intfuncIM5(i)=sin(func1st)

!****************************************************

!calculate the generalized phase at the final time t0
qgen(i)=A0(3)-lambda2*A0(1)/A0(2)

!calculate the derivative of generalized phase
divqgen(i)=A0(6)-lambda2*(A0(8)/A0(2)-A0(1)*A0(9)/A0(2)**2)

!save resolution values
    do k=1,neq
        A0matrix(k,i)=A0(k)
    end do
!***************************************************

!write the information in functioncheck.txt file
!ordered: initial phase, vibrational energy, action, dqdqi
!write(17,'(f10.3 f15.7 f15.7 f15.7 f15.7 f15.7 f15.7 f15.7)') phase(i), A0(4), A0(5),&
!A0(6), qgen(i)-phase(i), sqrt(abs(divqgen(i)))


	end do !phase change


!*********************************
!subtract the constant

 call cubint(n2+1, phase, qgen-phase, 1, n2+1, constantq, err1)
 print*, constantq/(2*pi)


    do i=1,n2+1
        qgen(i)=qgen(i)-constantq/(2*pi)

        !calculate under sin and cos expressions
        !ni --> nf
        partS(i)=exponentfunc(A0matrix(1:neq,i),qgen(i),inst,finst,actin)
        partS(i)=partS(i)+(A0matrix(4,i))*constantq/(2*pi)
        !nf -->ni
        partS1(i)=exponentfunc(A0matrix(1:neq,i),qgen(i),finst,inst,actin)
        partS1(i)=partS1(i)+(A0matrix(4,i))*constantq/(2*pi)

        write(17,'(f10.3 f15.7 f15.7 f15.7 f15.7 f15.7 f15.7 f15.7)') phase(i), qgen(i)-phase(i),&
        partS(i), partS1(i), divqgen(i)


        intfuncRE(i)=cos(partS(i))!*sqrt(abs(divqgen(i)))
        intfuncIM(i)=sin(partS(i))!*sqrt(abs(divqgen(i)))

        intfuncRE2(i)=cos(partS1(i))!*sqrt(abs(divqgen(i)))
        intfuncIM2(i)=sin(partS1(i))!*sqrt(abs(divqgen(i)))

    end do

!control if the constant ws subtracted in a correct way, want to get 0
call cubint(n2+1, phase, qgen-phase, 1, n2+1, constantq, err1)
print*, constantq/(2*pi)

!calculate the integrals
!real part
 call cubint(n2+1, phase, intfuncRe, 1, n2+1, ReS, err1)
 call cubint(n2+1, phase, intfuncRe1, 1, n2+1, ReS1, err1)
 call cubint(n2+1, phase, intfuncRe2, 1, n2+1, ReS2, err1)
 call cubint(n2+1, phase, intfuncRe3, 1, n2+1, ReS3, err1)
 call cubint(n2+1, phase, intfuncRe4, 1, n2+1, ReS4, err1)
 call cubint(n2+1, phase, intfuncRe5, 1, n2+1, ReS5, err1)

!imaginary part
 call cubint(n2+1, phase, intfuncIm, 1, n2+1, ImS, err2)
 call cubint(n2+1, phase, intfuncIm1, 1, n2+1, ImS1, err2)
 call cubint(n2+1, phase, intfuncIm2, 1, n2+1, ImS2, err2)
 call cubint(n2+1, phase, intfuncIm3, 1, n2+1, ImS3, err2)
 call cubint(n2+1, phase, intfuncIm4, 1, n2+1, ImS4, err2)
 call cubint(n2+1, phase, intfuncIm5, 1, n2+1, ImS5, err2)

!Probability

Prob=((ReS+ReS2)**2+(ImS+ImS2)**2)/(4*pi)**2
Prob1=(ReS**2+ImS**2+ReS2**2+ImS2**2)/(8*pi**2)
Prob2=((ReS1+ReS3)**2+(ImS1+ImS3)**2)/(4*pi)**2
Prob21=(ReS1**2+ImS1**2+ReS3**2+ImS3**2)/(8*pi**2)
Prob3=((ReS4+ReS5)**2+(ImS4+ImS5)**2)/(4*pi)**2
Prob31=(ReS4**2+ImS4**2+ReS5**2+ImS5**2)/(8*pi**2)


!calculate average energy loss of particle
 !call cubint(n2+1, phase, deltai, 1, n2+1, aven, err3)

! write(23,'(f15.7 f15.7)') actin, aven/(2*pi)

!deltaimax=maxval(deltai)

deallocate(phase)
deallocate(qgen,divqgen)
deallocate(partS, partS1)
deallocate(intfuncRE, intfuncIM)
deallocate(intfuncRE1, intfuncIM1)
deallocate(intfuncRE2, intfuncIM2)
deallocate(intfuncRE3, intfuncIM3)
deallocate(intfuncRE4, intfuncIM4)
deallocate(intfuncRE5, intfuncIM5)
deallocate(A0matrix)


print*, inst, '-->',  finst, ' ', '(numerical) Probability=', Prob
print*, inst, '-->',  finst, ' ', '(numerical) Probability1=', Prob1
print*, inst, '-->',  finst, ' ', '(analytical 1st) Probability=', Prob3
print*, inst, '-->',  finst, ' ', '(analytical 1st) Probability=', Prob31
print*, inst, '-->',  finst, ' ', '(analytical 2nd) Probability=', Prob2
print*, inst, '-->',  finst, ' ', '(analytical 2nd) Probability=', Prob21


write(10,'(f15.7 g15.4 g15.4 g15.4 g15.4 g15.4 g15.4)') energy-actin, Prob, Prob2, Prob3, &
Prob1, Prob21, Prob31

!end do !energy change
!end do !actin change

close(23)
close(17)
close(28)
close(30)
close(10)

end program
