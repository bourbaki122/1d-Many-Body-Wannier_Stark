MODULE nrtype
  INTEGER, PARAMETER :: ITT=SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RT=SELECTED_REAL_KIND(P=15)
  INTEGER, PARAMETER :: RT2=SELECTED_REAL_KIND(P=31)
  INTEGER, PARAMETER :: RTC=SELECTED_REAL_KIND(P=15)
  INTEGER, PARAMETER :: CT=RT
  REAL(rt),PARAMETER :: PI=3.14159265359
  REAL(rt),PARAMETER :: Euler=2.71828182846
  REAL(rt) :: resultado !look a beter place for put this one
END MODULE nrtype

MODULE single_particleP
USE nrtype
INTEGER,PARAMETER :: SITES=5 !6
REAL(rt),PARAMETER :: hopping=0.033_rt !0.033definite like positive see SP_eigensystem
REAL(rt),DIMENSION(sites,sites):: SP_eigenvectors
REAL(rt),DIMENSION(SITES) :: Energylevel


END MODULE single_particleP

MODULE many_particleP
USE nrtype
INTEGER,PARAMETER :: PARTICLES=4 !5
INTEGER(ITT) :: HDIMENSION
REAL(rt) :: Real_HDIMENSION
INTEGER, ALLOCATABLE :: NBase(:,:)
REAL(rt), ALLOCATABLE :: spectrum(:)
REAL(rt),PARAMETER :: Wpot=0.033_rt!0.0038_rt !0.3074!  (0.085_rt ! good one) 3p-5s
END MODULE many_particleP

MODULE HamilMatrix
USE nrtype
REAL(rt),ALLOCATABLE :: Ma_trix(:,:),Md_trix(:,:),My_trix(:,:)
REAL(rt),ALLOCATABLE :: Mb_trix(:,:),Mc_trix(:,:),Mx_trix(:,:),Mz_trix(:,:)
REAL(rt),ALLOCATABLE :: BH_trix(:,:),ergodic(:,:)
REAL(rt),ALLOCATABLE :: Diag_local(:), Diag_nonlocal(:)
REAL(rt),ALLOCATABLE :: Vij(:,:)
END MODULE HamilMatrix

!===============================================================================
! steps:
!===============================================================================
!*remember always calling SP_eigensystem before whatever change in the hopping
!or field is made, but try not put it inside others subroutines unless that is
!necesary
!
!*All the s¿matrix search and matrixx storage process is in Smart_decom.f90
!
!*remember that the sum of the ergodic and integrable part is in the subroutine
! hamiltonian
!
!*finally call hamiltonian_diagonalization
!
!*TAKE ACCOUNT THAT LAPACK CHANGE THE INPUT MATRIX INTO A UPPER DIAGONAL
!*AND THE EIGENVECTORS ARE THE COLUMS OF THE OUTPUN MATRIX I.E EIGENVALUE(1)
!*CORREPONDS TO EIGENVECTOR(:,1)
!
!*IMPORTANT: Cuando termine el trabajo implementar la fuencion de busqueda f(n1,n2,..)
!la cual me lleva de las coordenanas en el espacio de fock a un valor concreto
!en la base, esto ahorraria demasiado tiempo de computo y alocacion de ram.
!===============================================================================

PROGRAM kernelDynamics

USE nrtype
USE integrator
USE combinatories
USE single_particleP
USE many_particleP
USE HamilMatrix
USE class_info
USE Interactions
IMPLICIT NONE

INTEGER :: i,j
REAL(rt), EXTERNAL :: M_local,M_global

!HDIMENSION=factorial(PARTICLES+SITES-1)/(factorial(PARTICLES)*factorial(SITES-1))
Real_HDIMENSION=GAMMA(REAL(PARTICLES+SITES-1+1))/(GAMMA(REAL(particles+1))*GAMMA(REAL(SITES-1+1)))
HDIMENSION=INT(Real_HDIMENSION)
WRITE(*,*)"====IMPORTANT INFO===="
WRITE(*,*)"DIMENSION:", HDIMENSION

ALLOCATE(NBase(HDIMENSION,SITES))

!!------------reads the basis result from dummy ----------!!
  OPEN(UNIT=11, FILE="basisprg.dat")
  DO i = 1,HDIMENSION
      READ(11,*) (NBase(i,j),j=1,SITES)
  END DO
!!---------------------------------------------!!

  ALLOCATE(BH_trix(HDIMENSION,HDIMENSION))
  ALLOCATE(ergodic(HDIMENSION,HDIMENSION))
  ALLOCATE(Diag_local(HDIMENSION))
  ALLOCATE(Diag_nonlocal(HDIMENSION))
  ALLOCATE(Vij(HDIMENSION,HDIMENSION))


  ALLOCATE(Ma_trix(HDIMENSION,HDIMENSION))
  ALLOCATE(Md_trix(HDIMENSION,HDIMENSION))
  ALLOCATE(My_trix(HDIMENSION,HDIMENSION))
  ALLOCATE(Mb_trix(HDIMENSION,HDIMENSION))
  ALLOCATE(Mx_trix(HDIMENSION,HDIMENSION))
  ALLOCATE(Mc_trix(HDIMENSION,HDIMENSION))
  ALLOCATE(Mz_trix(HDIMENSION,HDIMENSION))


  WRITE(*,*)"Check Moves"

  CALL two_moves
  CALL Three_moves    !call all the smart ones
  CALL four_moves
  !Esta subrutina esta mal optimizada y para N/L=7/7 ocupa 16 gb de ram
  !ademas es reemplazable por la descomposicion en movimientos de sitios
  !CALL Complete_twobody !no se debe de desmarcar al menos de que Vij=0
  WRITE(*,*)"Alright"
  WRITE(*,*)"====End IMPORTANT INFO===="

!WRITE(*,*)"CHECK THE DECOMPOSITION"
!WRITE(*,*)"===========My_infor(it,SITES+1,1)================================"
!DO i=1,full_index!maxindex_y!full_index
!WRITE(*,*) Mcomplete_infor(i,:),i!My_infor(i,:,1), i
!END DO
!WRITE(*,*)"WITH THEIR RESPECTIVE COUPLES"
!DO i=1,maxindex_A
!WRITE(*,*) MA_infor(i,:,2), i
!END DO
!WRITE(*,*)"============================================"

CALL hamiltonian_diagonalization
!call EXPORT

END PROGRAM kernelDynamics

!===========================================================================
!subroutine SP_eigensystem:
!solves the single-body wannier stark sytem for a field f and a hopping j.
!lapack stores the eingnevector in the colums of the output Ma_trix
!also take account that the eigenvalues are sorting in acending order
!===========================================================================

SUBROUTINE SP_eigensystem(j,f)

  USE nrtype
  USE single_particleP
  USE many_particleP
  IMPLICIT NONE

  INTEGER :: i,k
  REAL(rt) :: j,f,j1,f1

  !Lapack Parametes
  REAL(rt),ALLOCATABLE :: work (:,:) !is a double type originally
  INTEGER :: lwork
  INTEGER :: info, LDA

  SP_eigenvectors=0.0_rt

  j1=j/Wpot !use this if you want to reescale H/W
  f1=f/Wpot !use this if you want to reescale H/W

!====here the single particle hamiltonian is defined...

  DO i=1,SITES
    DO k=1,SITES

      IF(i.EQ.k)THEN
        SP_eigenvectors(i,k)=i*f1
      END IF

      IF(ABS(i-k).EQ.1)THEN
        SP_eigenvectors(i,k)=(-1.0_rt)*j1/2.0_rt
      END IF

    END DO
  END DO
!====here the single particle hamiltonian was defined...

  LDA=SITES
  lwork = MAX(1,3*SITES-1)
  ALLOCATE(work(1:lwork,1:lwork))

  !WRITE(*,*)"CHECK THE MATRIX BEFORE"
  !DO i=1,SITES
  !  WRITE(*,*) SP_eigenvectors(i,:)
  !END DO

  CALL dsyev("V","U",SITES,SP_eigenvectors,LDA,Energylevel,work,lwork,info)
  WRITE(*,"(3X,A,I3)") 'Diagonalization performed, info equals ',info

  !WRITE(*,*)"CHECK THE MATRIX After"
  !DO i=1,SITES
!    WRITE(*,*) SP_eigenvectors(i,:)
  !END DO

 END SUBROUTINE SP_eigensystem

 !===========================================================================
 !functions M-type coefficient:
 !computes de m-type coeficients for a "n" level in the first band, with a
 !hopping term "j" and a external field "f" and a nearest neighbor "delta"
 !===========================================================================
 REAL(rt) FUNCTION M_global(a,b,c,d)
   USE nrtype
   USE single_particleP
   IMPLICIT NONE

   INTEGER :: i,a,b,c,d
   REAL(rt) :: sum,aux

   sum=0.0_rt

   DO i=1,SITES
      aux=SP_eigenvectors(i,a)*SP_eigenvectors(i,b)
      aux=aux*SP_eigenvectors(i,c)*SP_eigenvectors(i,d)
      sum=sum+aux
   END DO

   M_global=sum
 END FUNCTION M_global


 REAL(rt) FUNCTION M_local(tipo,n,delta)
  USE nrtype
  USE single_particleP

  CHARACTER (len=1) :: tipo
  INTEGER :: delta,n
  REAL(rt) :: aux
  REAL(rt), EXTERNAL :: M_global

  aux=0.0_rt

  IF(tipo.EQ."O") THEN
     aux=M_global(n,n,n,n) !M0
  ELSE IF(tipo.EQ."A")THEN
     aux=M_global(n,n,n,n+delta) !MA
  ELSE IF(tipo.EQ."D")THEN
     aux=M_global(n,n+delta,n+delta,n+delta) !MD
  ELSE IF(tipo.EQ."B")THEN
     aux=M_global(n,n,n+delta,n+(2*delta)) !MB
  ELSE IF(tipo.EQ."C")THEN
     aux=M_global(n,n+delta,n+(2*delta),n+(2*delta)) !MC
  ELSE IF(tipo.EQ."X")THEN
     aux=M_global(n,n+delta,n+delta,n+(2*delta)) !MX
  ELSE IF(tipo.EQ."Y")THEN
     aux=M_global(n,n,n+delta,n+delta) !MY
  ELSE IF(tipo.EQ."Z")THEN
     aux=M_global(n,n+delta,n+(2*delta),n+(3*delta)) !Mz
  ELSE
     WRITE(*,*)"========================="
     WRITE(*,*)"ERROR"
     WRITE(*,*)"========================="
  END IF

  M_local=aux

END FUNCTION M_local

!============================================================
!Subroutine hamiltonian
!computes the hamiltonian of the system
!===========================================================

SUBROUTINE hamiltonian(job)

USE nrtype
USE single_particleP
USE many_particleP
USE class_info
USE HamilMatrix
USE combinatories
USE Interactions
IMPLICIT NONE

INTEGER :: it,jt,delta,gamma,rho
REAL(rt),EXTERNAL :: M_local
REAL(rt) :: job

ergodic=0.0_rt
BH_trix=0.0_rt
Diag_local=0.0_rt
Diag_nonlocal=0.0_rt
Vij=0.0_rt
My_trix=0.0_rt
Ma_trix=0.0_rt
Md_trix=0.0_rt
Mb_trix=0.0_rt
Mc_trix=0.0_rt
Mx_trix=0.0_rt
Mz_trix=0.0_rt

  !CALL Two_body(job)

 DO it=1,HDIMENSION
   DO jt=1,SITES
   Diag_local(it)=Diag_local(it)+((job/2.0_rt)*(M_local("O",jt,1)*NBase(it,jt)*&
   (NBase(it,jt)-1)))+(Energylevel(jt)*NBase(it,jt))
   END DO
 END DO

 DO it=1,HDIMENSION
   DO jt=1,SITES
     DO delta=1,SITES-jt
       Diag_nonlocal(it)=Diag_nonlocal(it)+((job/2.0_rt)*4*M_local("Y",jt,delta)*&
       NBase(it,jt)*NBase(it,jt+delta))
     END DO
   END DO
 END DO

  DO delta=1,SITES
     CALL MA_matrix_smart(job,delta)
     CALL MD_matrix_smart(job,delta)
     CALL MY_matrix_smart(job,delta)
     !ergodic=ergodic+Ma_trix+Md_trix+My_trix
  END DO

DO delta=1,SITES
  DO gamma=1,sites
     CALL MC_matrix_smart(job,delta,gamma)
     CALL MB_matrix_smart(job,delta,gamma)
     CALL MX_matrix_smart(job,delta,gamma)
     !ergodic=ergodic+Mb_trix+Mc_trix+Mx_trix
  END DO
END DO

DO delta=1,sites
  DO gamma=1,sites
    DO rho=1,sites
      CALL Mz_matrix_smart(job,delta,gamma,rho)
       !ergodic=ergodic+Mz_trix
    END DO
  END DO
END DO

  DO it=1,HDIMENSION
   BH_trix(it,it)=Diag_local(it)+Diag_nonlocal(it)
  END DO

  BH_trix=BH_trix+ergodic

WRITE(*,*) " "
WRITE(*,*) "MAX/MIN VALUES OF THE INTEGRABLE PART"
WRITE(*,*) MAXVAL(Diag_local+Diag_nonlocal),MINVAL(Diag_local+Diag_nonlocal)

WRITE(*,*) "MAX/MIN VALUES OF THE Ergodic PART"
WRITE(*,*) MAXVAL(ergodic),MINVAL(ergodic)
WRITE(*,*) " "

END SUBROUTINE hamiltonian

!============================================================
!Subroutine hamiltonian-diagonalization
!makes the diagonalization of the hamiltonian for different values of the
!field
!===========================================================

SUBROUTINE hamiltonian_diagonalization

  USE nrtype
  USE integrator
  USE combinatories
  USE Thermalization_routines
  USE single_particleP
  USE many_particleP
  USE HamilMatrix

  IMPLICIT NONE

  INTEGER :: it,u=11
  INTEGER, PARAMETER :: field_iter=1000!700
  REAL(rt), PARAMETER :: delta_it=0.0001_rt
  REAL(rt) :: f
  REAL(rt),EXTERNAL :: CHAIN,M_local
  REAL(rt), ALLOCATABLE :: spacings(:), Initial_survival(:)
  REAL(rt), DIMENSION(HDIMENSION) :: dummy_vec

  !Lapack Parametes
  REAL(rt),ALLOCATABLE :: worktwo(:,:)
  INTEGER :: lworktwo
  INTEGER :: infotwo, LDAtwo

  ALLOCATE(spectrum(HDIMENSION))
  ALLOCATE(Initial_survival(HDIMENSION))
  ALLOCATE(spacings(HDIMENSION-1))

  dummy_vec=0.0_rt
  Initial_survival=0.0_rt
  Initial_survival(1597)=1.0_rt  !definition of the initial state

  LDAtwo=HDIMENSION
  lworktwo = MAX(1,3*HDIMENSION-1)

  ALLOCATE(worktwo(1:lworktwo,1:lworktwo))

  f=0.0_rt

  OPEN(newunit=u, file="Spectrum_vs_field.dat", status="replace")

  DO it=1,field_iter
    f=0.0_rt+(it*delta_it)   !0.22 rare

  CALL SP_eigensystem(hopping,f)
  CALL hamiltonian(1.0_rt) !put 1 if you are using H/W

  CALL dsyev("V","U",HDIMENSION,BH_trix,LDAtwo,spectrum,worktwo,lworktwo,infotwo)
  WRITE(*,"(3X,A,I3)") 'Diagonalization performed, info equals ',infotwo

  WRITE(u,*) f,spectrum*wpot
  write(*,*) it
  END DO

  CLOSE(u)

  !From here and further i compute the density of eignestates and then
  !starcaise function from a single vale of the parameter space

  spectrum=0.0_rt

  write(*,*)"  "
  write(*,*)"¡=====SPECTRAL PART=====!"
  write(*,*)"  "
  write(*,*)"Single Particle Diagonalization"
  CALL SP_eigensystem(hopping,0.0013_rt) !0.02 rare !this controls the histogram
  CALL hamiltonian(1.0_rt) !put 1 if you are using H/W
  write(*,*)"  "
  write(*,*)"Many Particle Diagonalization"
  CALL dsyev("V","U",HDIMENSION,BH_trix,LDAtwo,spectrum,worktwo,lworktwo,infotwo)
  WRITE(*,"(3X,A,I3)") 'Diagonalization performed, info equals ',infotwo

  CALL quicksort(spectrum)
  dummy_vec=spectrum
  WRITE(*,*) spectrum(1),SPECTRUM(HDIMENSION)
  CALL density_of_states(dummy_vec,HDIMENSION,35)
  CALL spacing_distribution(dummy_vec,HDIMENSION,0.10_rt,5,30,0.1_rt)
  write(*,*)"  "
  write(*,*)"=====¡DYNAMICAL PART!====="
  !call Spectral_measures(spectrum,HDIMENSION)
  !call lorentz_density(spectrum,HDIMENSION,spectrum(1),spectrum(HDIMENSION),10,"T")
  !CALL survival_p(Initial_survival,0.0_rt,10.0_rt,spectrum,BH_trix,200)
END SUBROUTINE

!===============================================================
!Export
!*this subroutine contains the matrix that im going to export to mathematica
! for make a matrix plot
!==============================================================

SUBROUTINE Export
  USE nrtype
  USE single_particleP
  USE many_particleP
  USE class_info
  USE HamilMatrix

  IMPLICIT NONE

  INTEGER :: i,j
  INTEGER, ALLOCATABLE :: matplot(:,:)
  ALLOCATE(matplot(HDIMENSION,HDIMENSION))
  matplot=0

  CALL SP_eigensystem(hopping,0.0013_rt)
  CALL hamiltonian(wpot) !here it dosent matter to put 1
  !CALL Two_body(wpot) !here it dosent matter to put 1

  !write(*,*)"EXPORT"

  DO i=1,HDIMENSION
    DO j=1,HDIMENSION
      !write(*,*)Vij(i,j)
      IF(ABS(Vij(i,j)).GT.0.0_rt)THEN
        matplot(i,j)=1
      END IF
  END DO
END DO

  OPEN(UNIT=12, FILE="Total_matplot.txt", ACTION="write", STATUS="replace")
 DO i=1,HDIMENSION
   WRITE(12,*) (matplot(i,j), j=1,HDIMENSION)
 END DO

 matplot=0

 DO i=1,HDIMENSION
   DO j=1,HDIMENSION
     !write(*,*) BH_trix(i,j)
     IF(ABS(BH_trix(i,j)).GT.0.0_rt)THEN
       matplot(i,j)=1
     END IF
   END DO
 END DO

OPEN(UNIT=13, FILE="Partial_matplot.txt", ACTION="write", STATUS="replace")
DO i=1,HDIMENSION
 WRITE(13,*) (matplot(i,j), j=1,HDIMENSION)
END DO

END SUBROUTINE

!==============================================================
!FUNCTION lorentz_density:
!Gives the density of levels by taking the lorenz aproximation
!p(E)=sum(y/(E-E_{n})^2 + y^2)) , whit y= mean value of the spacings
!check Parra Phd Thesis for more INFO.
!VARIABLES:
!input: spectrum of the hamiltonian
!in_dim : dimension of input
!sigma : mean value of the original spacings (y)
!energy : E (independent variable), p(E)
!==================================================================

SUBROUTINE lorentz_density(input,in_dim,a0,b0,N0,add_plot)

  USE nrtype
  USE integrator
  USE many_particleP
  IMPLICIT NONE

  CHARACTER (len=1) :: add_plot
  INTEGER ::it,u=3
  INTEGER, INTENT(in) :: N0,in_dim
  REAL(rt), INTENT(in) :: a0,b0
  REAL(rt) :: sigma0,x
  REAL(rt), DIMENSION(in_dim) :: input


  sigma0=0.0_rt

  DO it=1,in_dim-1
  sigma0=input(it+1)-input(it)+sigma0
  END DO

  sigma0=sigma0/(0.8*(in_dim-1)) !modified sigma0
  resultado=Simpson_peak(rho,a0,b0,N0)


IF(add_plot.EQ."T")THEN

OPEN(newunit=u, file="density_levels.dat", status="replace")

  DO it=1,10000
    x=input(1)+((input(in_dim)-input(1))/10000.0_rt)*it
    WRITE(u,*) x,rho(x)
  END DO

  CLOSE(u)
END IF

  CONTAINS

  REAL(rt) FUNCTION rho(x) RESULT(p)
    REAL(rt), INTENT(in) :: x
    INTEGER :: it

    p=0.0_rt
    DO it=1,in_dim
    p=p+sigma0/(((x-input(it))**2)+(sigma0**2))
    END DO
    p=(p/(PI))   !this definition doesnt
                 !have the normalization dim(hilbert)
  END FUNCTION

END SUBROUTINE lorentz_density


SUBROUTINE Spectral_measures(vector,vec_dim)
USE nrtype
USE integrator
USE combinatories
USE Thermalization_routines
IMPLICIT NONE

REAL(rt) :: vector(*)
REAL(rt) :: error,f_cost,aux1,aux2,x,delta_rms
INTEGER ::  vec_dim, adap_step, max_iter,it,u,control_count,up
REAL(rt), ALLOCATABLE :: neighbor_levels(:)

ALLOCATE(neighbor_levels(vec_dim-1))
!CALL quicksort(vector)

control_count=1
delta_rms=0.0_rt
up=0
neighbor_levels=0.0_rt
adap_step=vec_dim+3  !initial number of steps, must be a multiple of 3
f_cost=0.0_rt  !cost function
error=1.0E-9_rt !< 1E-7 t>>100000
aux1=0.0_rt
aux2=0.0_rt
u=2
max_iter=1000000 ! max bound of iterations if the algorithm doesnt converge

DO it=1,vec_dim-1

  DO WHILE(adap_step < max_iter)

   CALL lorentz_density(vector,vec_dim,vector(it),vector(it+1),adap_step,"F")
   aux1=resultado
   CALL lorentz_density(vector,vec_dim,vector(it),vector(it+1),2*adap_step,"F")
   aux2=resultado

   adap_step=2*adap_step
   f_cost=ABS(aux2-aux1)
   !write(*,*)adap_step,f_cost,error
   IF(f_cost<error) EXIT
   !write(*,*)adap_step,f_cost,error
  END DO

neighbor_levels(it)=resultado
adap_step=3+vec_dim

END DO

aux1=0
aux2=0
adap_step=3+vec_dim
f_cost=0.0_rt

CALL easy_comulative(neighbor_levels)

OPEN(newunit=u, file="staircase_func.dat", status="replace")
DO it=1,20000

  x=vector(1)+((vector(vec_dim)-vector(1))/20000.0_rt)*it

  DO WHILE(adap_step < max_iter)
   CALL lorentz_density(vector,vec_dim,vector(1),x,adap_step,"F")
   aux1=resultado
   CALL lorentz_density(vector,vec_dim,vector(1),x,2*adap_step,"F")
   aux2=resultado
   adap_step=2*adap_step
   f_cost=ABS(aux2-aux1)
   !write(*,*)adap_step,f_cost,error
   IF(f_cost<1.0E-5_rt) EXIT
   WRITE(*,*)adap_step,f_cost!,error
  END DO

  adap_step=3+vec_dim

   IF(x.LE.vector(control_count))THEN
     up=control_count
   ELSE IF(x.GT.vector(control_count))THEN
     control_count=control_count+1
     up=control_count
   END IF

WRITE(u,*) x,resultado,up-1
delta_rms=delta_rms+((resultado-(up-1))**2)
END DO

CLOSE(u)

delta_rms=delta_rms/20000

WRITE(*,*)"========la desviacion es==================="
WRITE(*,*)delta_rms
WRITE(*,*)"==========================="

END SUBROUTINE Spectral_measures



SUBROUTINE density_of_states(espectro,espectro_dim,nbins)

  USE nrtype
  USE many_particleP

  INTEGER :: nbins,espectro_dim
  REAL(rt),DIMENSION(espectro_dim) :: espectro
  REAL(rt), ALLOCATABLE ::  midenergy(:),countbin(:)
  REAL(rt) :: bin


  ALLOCATE(midenergy(nbins+1))
  ALLOCATE(countbin(nbins))
  midenergy=0.0_rt
  countbin=0

  bin=(espectro(SIZE(espectro))-espectro(1))/REAL(nbins)

  DO it=1,nbins+1
    midenergy(it)=espectro(1)-bin-(bin/2.0_rt)+(bin*REAL(it))
  END DO

  DO it=1,SIZE(espectro)
    DO jt=1,nbins+1
      IF(espectro(it)>=midenergy(jt).AND.espectro(it)<midenergy(jt+1))THEN
      countbin(jt) = countbin(jt)+1
      END IF
    END DO
  END DO

WRITE(30,130) midenergy(1),0.0d0
  DO kt=1,nbins
   WRITE(30,130) midenergy(kt),REAL(countbin(kt))/(bin*REAL(SIZE(espectro)))
   WRITE(30,130) midenergy(kt+1),REAL(countbin(kt))/(bin*REAL(SIZE(espectro)))
  END DO
130 FORMAT (2E18.9)

END SUBROUTINE density_of_states


SUBROUTINE spacing_distribution(espectro,espectro_dim,porcentaje,nblocks,nbins,binsize)

USE nrtype
USE many_particleP
USE Thermalization_routines

INTEGER :: nbins,espectro_dim,nblocks,half,it,jt
REAL(rt),DIMENSION(espectro_dim) :: espectro
REAL(rt) :: efective_dim,boundary,binsize,mean_energy,dist_norm,porcentaje
REAL(rt), ALLOCATABLE :: space_dis(:),spgrid(:)
INTEGER, ALLOCATABLE :: spcounts(:)

boundary=porcentaje*SIZE(espectro)
half=INT(boundary/2.0_rt)
efective_dim=REAL(SIZE(espectro))-boundary

ALLOCATE(space_dis(INT(efective_dim)))


efective_dim=efective_dim/nblocks
DO it=1,INT(efective_dim)
mean_energy=(espectro(half+(nblocks*it))-espectro(half+(nblocks*(it-1))))/nblocks
DO jt=1+(nblocks*(it-1)),nblocks*it
  space_dis(jt)=(espectro(half+jt)-espectro(half+jt-1))/mean_energy
END DO
END DO

ALLOCATE(spgrid(nbins+1))
ALLOCATE(spcounts(nbins))

spgrid=0.0_rt
DO it=1,nbins
  spgrid(it+1)=spgrid(it)+binsize
END DO
spcounts=0

DO it=1,nblocks*INT(efective_dim)
  DO jt=1,nbins
    IF((space_dis(it).GE.spgrid(jt)).AND.(space_dis(it).LT.spgrid(jt+1)))THEN
      spcounts(jt)=spcounts(jt)+1
    END IF
  END DO
END DO

dist_norm=0.0_rt

DO it=1,nbins
dist_norm=dist_norm+(binsize*spcounts(it))
END DO

WRITE(40,130) spgrid(1),0.0_rt
DO it=1,Nbins
  WRITE(40,130) spgrid(it),spcounts(it)/dist_norm
  WRITE(40,130) spgrid(it+1),spcounts(it)/dist_norm
END DO
130 FORMAT (2E18.9)

CALL easy_comulative(space_dis)

END SUBROUTINE spacing_distribution
