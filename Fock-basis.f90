MODULE casi
  INTEGER, ALLOCATABLE :: numbe(:,:)!,numbe0(:)
  INTEGER :: n,m,jmax,bb,cmax,maxnzr

  CHARACTER(1) :: onlyje,ordbas
  !--------------ex-casi-------------------------------------
  INTEGER, ALLOCATABLE :: basis(:,:),msf(:),seedvf(:,:),seedinfof(:,:)
  INTEGER ::cont
END MODULE casi

PROGRAM DUMMYs
USE casi
IMPLICIT NONE

n=4 !5
m=5 !6

ordbas='y'

CALL findbasisn

END PROGRAM DUMMYs

SUBROUTINE findbasisn

  USE casi

  IMPLICIT NONE
  INTEGER:: sn2,n2
  INTEGER :: k,contfbasis,l,jmaxm1,bl,i,sn,ssn,ne,ne0
  INTEGER :: acum,acum0f,acum0,acumf
  INTEGER, ALLOCATABLE :: basis0(:,:)

  n2=n

  IF(ABS(n).GE.128.OR.ABS(m).GE.128)THEN
     WRITE(*,*)'abs(n).ge.128.or.abs(m).ge.128',n,m
     WRITE(*,*)'Some integers here are of the type I1B = selected_int_kind(2),'
     WRITE(*,*)'i.e. <128. PROGRAM FINISHED'
     STOP
  END IF

  !calculo de jmax=(n+m-1)!/(n!(m-1)!)
  ALLOCATE(numbe(1:m,0:n))
  numbe(:,0)=1
  DO sn=1,n
     numbe(m,sn)=numbe(m,sn-1)*(m+sn-1)/sn
  END DO
  jmaxm1=numbe(m,n-1);jmax=numbe(m,n)
  WRITE(*,*)'jmax,jmaxm1',jmax,jmaxm1
  ALLOCATE(basis(jmax,m))
  basis=0
  !l=number of cells in the optical lattice
  !sn=number of bosons in the lattice
  !construction of the basis for l=1 cells and sn=0,1,..,n bosons
  !numbe(sn)=contains the number of basis elements in the
  !          for l cells and sn particles.
  !numbe0(sn)=number of basis elemmala vidaents fo l-1 ans sn bosons
  l=1
  DO sn=0,n
     sn2=sn  !integer*2
     !bl=position in the array basis(:,:)
     bl=m-l+1
     numbe(l,sn)=1
     basis(sn+1,bl)=n2-sn2   !integer*2
  END DO
  !Empieza la induccion:
  DO l=2,m-1
     bl=m-l+1
     !construction of numbe(m,sn)
     numbe(l,0)=1
     DO sn=1,n
        numbe(l,sn)=numbe(l,sn-1)*(l+sn-1)/(sn)
     END DO
     contfbasis=0
     !construction of the basis for n bosons and l cells
     sn=0  !it it not necessary to iterate (basis=0)
     acum=numbe(l-1,n)
     contfbasis=acum
     DO sn=1,n
        sn2=sn
        ne0=numbe(l-1,n-sn)
        !recorremos los elementos de la base para l-1 y sn
        DO k=1,ne0
           contfbasis=contfbasis+1
           basis(contfbasis,bl)=sn2
        END DO
        acum=acum+ne0 !=contfbasis
     END DO
     !acum0=numbe0(n)
     acum0=numbe(l-1,n)
     !construction of the basis for n-1,..,0 bosons and l cells
     DO ssn=n-1,1,-1
        !        write(*,*)
        !        write(*,*)ssn,numbe
        !        ne=numbe(ssn)
        ne=numbe(l,ssn)
        acumf=acum+ne
        acum0f=acum0+ne
        basis(acum+1:acumf,bl+1:m)=basis(acum0+1:acum0f,bl+1:m)
        !sn=0  !it it not necessary to iterate (basis=0)
        !acum=numbe0(n)
        !        contfbasis=contfbasis+numbe0(ssn)
        contfbasis=contfbasis+numbe(l-1,ssn)
        DO sn=1,ssn
           sn2=sn
           !           ne0=numbe0(ssn-sn)
           ne0=numbe(l-1,ssn-sn)
           !recorremos los elementos de la base para l-1 y sn
           DO k=1,ne0
              contfbasis=contfbasis+1
              basis(contfbasis,bl)=sn2
           END DO
        END DO
        acum=acumf
        !        acum0=acum0+numbe0(ssn)
        acum0=acum0+numbe(l-1,ssn)
     END DO
     ssn=0
     acum=acumf
     !     ne=numbe(ssn)
     ne=numbe(l,ssn)
     acumf=acum+ne
     contfbasis=contfbasis+1
  END DO

  !  finally l=m
  l=m
  bl=m-l+1
  !construction of numbe(sn)
  !  numbe0=numbe
  !  numbe(0)=1
  numbe(l,0)=1
  DO sn=1,n
     !     numbe(sn)=numbe(sn-1)*(l+sn-1)/(sn)
     numbe(l,sn)=numbe(l,sn-1)*(l+sn-1)/(sn)
  END DO
  acum=0
  contfbasis=0
  !construction of the basis for n bosons and l cells
  !sn=0-->it isn't necessary to iterate
  !  contfbasis=contfbasis+numbe0(n)
  contfbasis=contfbasis+numbe(l-1,n)
  DO sn=1,n
     sn2=sn
     !     ne0=numbe0(n-sn)
     ne0=numbe(l-1,n-sn)
     !recorremos los elementos de la base para l-1 y sn
     DO k=1,ne0
        contfbasis=contfbasis+1
        basis(contfbasis,bl)=sn2
     END DO
     acum=acum+ne0
  END DO

  !write(*,*)'acumf,contfbasis',acumf,contfbasis
  !do i=1,acumf
  !   write(*,*)basis(i,:),'  sum=',sum(basis(i,:))
  !end do

  IF(contfbasis.NE.jmax) THEN
     WRITE(*,*)'contfbasis.ne.jmax',contfbasis,jmax
     WRITE(*,*)'PROGRAM FINISHED'
     STOP
  END IF
  WRITE(*,*)'basis size',jmax
  !BASIS AS IN MATLAB
  !  open(4,file='basisn.dat')
  !  do i=jmax,1,-1
  !     write(4,1)basis(i,:)
  !1    format(24i3)
  !  end do
  !  close(4)

  IF(ordbas.EQ.'y') THEN
     ALLOCATE(basis0(jmax,m))
     DO i=1,jmax
        basis0(i,:)=basis(jmax-i+1,:)
     END DO
     basis=basis0
     DEALLOCATE(basis0)
     WRITE(*,*)'BASIS AS IN MATLAB'

     OPEN(4,file='basisprg.dat')
     DO i=1,jmax
        WRITE(4,1)basis(i,:)
1       FORMAT(24i3)
     END DO
     CLOSE(4)

  END IF

  WRITE(*,*)'numbe(l,:) for l=1,..,m'
  DO l=1,m
     WRITE(*,*)numbe(l,:)
  END DO

  RETURN

END SUBROUTINE findbasisn
