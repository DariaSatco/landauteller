module landaumod

real(8), parameter:: alpha=0.3  !0.114 !0.314  !0.1278   !0.3
real(8), parameter:: lambda2=2.0D+00/3     !0.5D+00  !1.0D+00/3  !0.006268  !0.5D+00 !2.0D+00/3   !0.667 !0.2D+00
real(8), parameter:: pi=3.1415926
real(8), public:: U

!this module is written in oder to use it in resolvation of differential equations system
!we will use further a vector A, wich has the following components
!A=(A1,A2)
!A1=(R,P,q,I,phi) - vector of Hamiltomian variables and action
!A2=(dqdqi, dIdqi, dRdqi, dPdqi, dqdIi, dIdIi, dRdIi, dPdIi, dqdRi, dIdRi, dRdRi, dPdRi, ...
!... dqdPi, dIdPi, dRdPi, dPdPi)
!vector A consists of 20 components, all of them depends on time t

contains

subroutine derivfunc(m,t,vector,DifEq)
! subroutine  derivfunc( t, m, vector, DifEq) 
!    which evaluates the derivative vector(1:m) given the time t and
!    solution vector DifEq(1:m).

implicit none

integer:: m !the system dimension
real(8):: vector(m) !input vector
real(8):: DifEq(m) !output vector
real(8):: M1(4,4), M2(16,16) !matrices for the system
real(8):: t !time variable
real(8):: commonvar !an expression which we use a lot of times
integer:: i,j,k

!initialize first 5 equations for R, P, q, I, phi
commonvar=U*alpha*exp(-alpha*(vector(1)-sqrt(2*vector(4))*cos(vector(3))))

DifEq(1)=vector(2)/lambda2
DifEq(2)=commonvar
DifEq(3)=1.0D+00+commonvar*cos(vector(3))/sqrt(2*vector(4))
DifEq(4)=commonvar*sqrt(2*vector(4))*sin(vector(3))
DifEq(5)=-(vector(3)*DifEq(4)+vector(1)*DifEq(2))

!next 16 equation we will initialize in terms of a matrix

!note, that matrix M 16x16 consists of blocks, so initialize first a 4x4 matrix M1
 M1(1,1)=-cos(vector(3))*(alpha*sqrt(2*vector(4))*sin(vector(3))+tan(vector(3)))/sqrt(2*vector(4))
 M1(1,2)=cos(vector(3))*(alpha*cos(vector(3))/sqrt(2*vector(4))-1._8/(2*vector(4)))/sqrt(2*vector(4))
 M1(1,3)=-alpha*cos(vector(3))/sqrt(2*vector(4))
 M1(1,4)=0._8
 M1(2,1)=(-alpha)*(sqrt(2*vector(4))*sin(vector(3)))**2+sqrt(2*vector(4))*cos(vector(3))
 M1(2,2)=sqrt(2*vector(4))*sin(vector(3))*(alpha*cos(vector(3))/sqrt(2*vector(4))+ 1._8/(2*vector(4)))
 M1(2,3)=-alpha*sqrt(2*vector(4))*sin(vector(3))
 M1(2,4)=0._8
 M1(3,1)=0._8
 M1(3,2)=0._8
 M1(3,3)=0._8
 M1(3,4)=1._8/(lambda2*commonvar)
 M1(4,1)=-alpha*sin(vector(3))*sqrt(2*vector(4))
 M1(4,2)=alpha*cos(vector(3))/sqrt(2*vector(4))
 M1(4,3)=-alpha
 M1(4,4)=0._8

!initialize M2 matrix
	do i=1,16
		do j=1,16
		M2(i,j)=0._8
		end do
	end do

do k=0,3
	do i=1,4
		do j=1,4
		M2(4*k+i,4*k+j)=M1(i,j)
		end do
	end do
end do

M2=commonvar*M2

!initialize DifEq 6-21 components

do k=1,16

DifEq(k+5)=0._8

	do i=1,16
	DifEq(k+5)=DifEq(k+5)+M2(k,i)*vector(i+5)
	end do
end do

return

end subroutine

!------------------------------------------------------------------------
!subroutine calculates the right part of ode system for first five equations
!of my problem
subroutine derivfunc5(m,t,vector,DifEq)

implicit none

integer:: m !the system dimension
real(8):: vector(5) !input vector
real(8):: DifEq(5) !output vector
real(8):: t !time variable
real(8):: commonvar !an expression which we use a lot of times


!initialize first 5 equations for R, P, q, I, phi
commonvar=U*alpha*exp(-alpha*(vector(1)-sqrt(2*vector(4))*cos(vector(3))))

DifEq(1)=vector(2)/lambda2
DifEq(2)=commonvar
DifEq(3)=1.0D+00+commonvar*cos(vector(3))/sqrt(2*vector(4))
DifEq(4)=commonvar*sqrt(2*vector(4))*sin(vector(3))
DifEq(5)=-(vector(3)*DifEq(4)+vector(1)*DifEq(2))

return

end subroutine
!------------------------------------------------------------------------
!chenged hamiltonian for check program
function hamiltonianred(vector)

implicit none
real(8):: vector(5)
real(8):: hamiltonianred

hamiltonianred=vector(2)**2/(2*lambda2)+U*exp(-alpha*(vector(1)-sqrt(2*vector(4))*cos(vector(3))))+vector(4)

end function

!------------------------------------------------------------------------
SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
!a dummy subroutine for ode library use
                   INTEGER  NEQ, ML, MU, NROWPD
                   real(8):: T
                   real(8), allocatable:: Y(:), PD(:,:)

return
end subroutine
!------------------------------------------------------------------------

function exponentfunc(vector,qgen,n1,n2,actin)
!the function staying under exponent in equation for S-matrix element

implicit none
real(8):: vector(21)
real(8):: exponentfunc
real(8):: qgen, actin
integer:: n1,n2

exponentfunc=vector(5)+qgen*(vector(4) - actin + real(n1-n2))

end function exponentfunc

!----------------------------------------------------------------------
function actionfunc(vector)
!the function staying under the action integral

implicit none

real(8):: vector(21), DifEq(2)
real(8):: actionfunc, commonvar

commonvar=U*alpha*exp(-alpha*(vector(1)-sqrt(2*vector(4))*cos(vector(3))))

DifEq(1)=commonvar*sqrt(2*vector(4))*sin(vector(3))
DifEq(2)=commonvar

actionfunc=-(vector(3)*DifEq(1)+vector(1)*DifEq(2))

end function



!-----------------------------------------------------------------------
function hamiltonian(vector)

implicit none
real(8):: vector(21)
real(8):: hamiltonian

hamiltonian=vector(2)**2/(2*lambda2)+U*exp(-alpha*(vector(1)-sqrt(2*vector(4))*cos(vector(3))))+vector(4)

end function
!-----------------------------------------------------------------------

!Function to find the determinant of a square matrix
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements


FUNCTION FindDet(matrix, n)

    IMPLICIT NONE
    
    real(8):: FindDet
    REAL(8), DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL(8) :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0._8) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0._8) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0._8
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
   
END FUNCTION FindDet

end module landaumod
