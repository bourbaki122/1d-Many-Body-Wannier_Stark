!========================================================================
!In this files i gonna put all the usual mathematical tools that i ususlly use
!*Integrators
!*diferentiators
!*Evolutors
!*polynomials
!========================================================================

MODULE Integrator
  USE nrtype
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: simpson,simpson_38,simpson_peak

  INTERFACE
     REAL(rt) FUNCTION func(x)
       USE nrtype
       IMPLICIT NONE
       REAL(rt), INTENT(in) :: x
     END FUNCTION func
  END INTERFACE

CONTAINS

  REAL(rt) FUNCTION Simpson_peak(f,a,b,steps) RESULT(s)! must be a multple of 3

    INTEGER, INTENT(in) :: steps
    REAL(rt), INTENT(in) :: a,b
    INTEGER :: it
    REAL(rt) :: h,x,s_aux
    PROCEDURE(func) :: f

    h=(b-a)*1.0_rt/steps
    s=0.0_rt
    s_aux=0.0_rt

    DO it=0,steps

       x=a+(it*h)
       IF((it.EQ.0).OR.(it.EQ.steps))THEN
          s_aux=9.0_rt*f(x)
       ELSE IF((it.EQ.1).OR.(it.EQ.(steps-1)))THEN
          s_aux=28.0_rt*f(x)
       ELSE IF((it.EQ.2).OR.(it.EQ.(steps-2)))THEN
          s_aux=23.0_rt*f(x)
       ELSE
          s_aux=24.0_rt*f(x)
       END IF

       s=s+s_aux
    END DO

    s=(h/24.0_rt)*s
  END FUNCTION Simpson_peak

  REAL(rt) FUNCTION simpson(f,a,b,steps) RESULT(s) !steps must be even

    INTEGER, INTENT(in) :: steps
    REAL(rt), INTENT(in) :: a,b
    PROCEDURE(func) :: f
    INTEGER :: it
    REAL(rt) :: h,x,s_aux


    h=(b-a)*1.0_rt/steps
    !write(*,*)h
    s=0.0_rt
    s_aux=0.0_rt

    DO it=0,steps

       x=a+(it*h)
       IF(it.EQ.0) s_aux=f(a)
       IF(MOD(it,2).EQ.1) s_aux=4.0_rt*f(x)
       IF(MOD(it,2).EQ.0) s_aux=2.0_rt*f(x)
       IF(it.EQ.steps) s_aux=f(b)
       s=s+s_aux

    END DO

    s=(h/3.0_rt)*s

  END FUNCTION simpson

  REAL(rt) FUNCTION simpson_38(f,a,b,steps) RESULT(s) !steps must be multiple of 3

    INTEGER, INTENT(in) :: steps
    REAL(rt), INTENT(in) :: a,b
    PROCEDURE(func) :: f
    INTEGER :: it
    REAL(rt) :: h,x,s_aux


    h=(b-a)*1.0_rt/steps
    s=0.0_rt
    s_aux=0.0_rt

    DO it=0,steps
       x=a+(it*h)
       IF(it.EQ.0) s_aux=f(a)
       IF(MOD(it,3).GT.0) s_aux=3.0_rt*f(x)
       IF(MOD(it,3).EQ.0) s_aux=2.0_rt*f(x)
       IF(it.EQ.steps) s_aux=f(b)
       s=s+s_aux
    END DO

    s=(3.0_rt*h/8.0_rt)*s

  END FUNCTION simpson_38

  SUBROUTINE test_func(a,b,n)
    USE nrtype

    REAL(rt), INTENT(in) :: a,b
    INTEGER, INTENT(in) :: n
    WRITE(*,*) Simpson_peak(g,a,b,n)

  CONTAINS

    REAL(rt) FUNCTION g(x) RESULT(y)
      REAL(rt), INTENT(in) ::x
      y=DSin(x)**2
    END FUNCTION g
  END SUBROUTINE test_func

END MODULE integrator

!=========================================================================
!Module combinatories:
!here im going to put all the things related with factorial, sorting .and.
!stuff like those
!========================================================================0

MODULE combinatories

  USE nrtype
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: factorial, quicksort

CONTAINS

  !===========================================================
  !SUBROUTINE quicksort:
  !MAKES A SIMPLE SORTING ALGORITHM ON A VECTOR OF DIM=LAST
  !REMEMBER THAT Lapack SORTS THE SPECTRUM
  !===========================================================

  RECURSIVE SUBROUTINE quicksort(a)
    USE nrtype
    IMPLICIT NONE
    REAL(rt) :: a(:)
    REAL(rt) :: x, t
    INTEGER :: first = 1, last
    INTEGER :: i, j

    last = SIZE(a, 1)
    x = a( (first+last) / 2 )
    i = first
    j = last

    DO
       DO WHILE (a(i) < x)
          i=i+1
       END DO
       DO WHILE (x < a(j))
          j=j-1
       END DO
       IF (i >= j) EXIT
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    END DO

    IF (first < i - 1) CALL quicksort(a(first : i - 1))
    IF (j + 1 < last)  CALL quicksort(a(j + 1 : last))
  END SUBROUTINE quicksort

  INTEGER(ITT) FUNCTION factorial(n)
    USE nrtype
    IMPLICIT NONE
    INTEGER :: n,m
    INTEGER(ITT) :: l

    l=1
    DO m=1,n
       l=m*l
    END DO
    factorial=l
  END FUNCTION factorial

END MODULE combinatories

!==============================================================
!MODULE: MATRIX MANIPULATION
!here im going to have all the matrix things,
!if its better if we latter implement lapack on this NONE
!==============================================================

MODULE Matrix_manipulation
  USE nrtype
  IMPLICIT NONE
  !PRIVATE
  !PUBLIC :: Exp_mat
  COMPLEX(rtc), ALLOCATABLE :: matExp(:,:) ! FOR THE EXP OF THE MATRIX

CONTAINS

  FUNCTION Exp_mat(A,U,t0,tf) RESULT(matexp)

    USE nrtype
    USE HamilMatrix
    IMPLICIT NONE

    REAL(rt), INTENT(in) :: t0,tf
    REAL(rt),DIMENSION(:,:), INTENT(in) :: U
    REAL(rt),DIMENSION(:), INTENT(in) :: A
    COMPLEX(rtc), DIMENSION(SIZE(A),SIZE(A)) :: matexp
    INTEGER :: dimA,it

    dimA=SIZE(A)
    matExp=DCmplx(0.0_RT,0.0_RT)

    DO it=1,dimA
       matExp(it,:)=DCmplx(DCos(A(it)*(tf-t0))*U(:,it),DSin(A(it)*(tf-t0))*U(:,it))
    END DO

    matexp=MATMUL(U,matexp)

  END FUNCTION Exp_mat

END MODULE Matrix_manipulation


MODULE Thermalization_routines
  USE nrtype
  USE Matrix_manipulation
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: survival_p,easy_comulative

CONTAINS

  !==========================================================================
  !Subroutine survival:
  !takes an vector evola it in time an see what is the probability of finding
  !on lattice sites, also computes the number of principal componets of the
  !vector, and their corresponding number distribution
  !==========================================================================

  SUBROUTINE survival_p(psi0,t0,tf,A,V,partitions)
    USE nrtype
    USE Matrix_manipulation
    USE single_particleP
    USE many_particleP
    IMPLICIT NONE

    REAL(rt), DIMENSION(:) :: psi0
    REAL(rt) :: t0,tf,t,NPC,NPC_stationary,graduar
    INTEGER :: partitions,it,jt,kt,zt,psi0_dim,u,u2,u3,u4
    REAL(rt),DIMENSION(:,:) :: V
    REAL(rt),DIMENSION(:) :: A
    REAL(rt), ALLOCATABLE :: psit(:),psi_diag(:),psi_temp(:)
    REAL(rt), ALLOCATABLE :: occupation_sites(:),occup_sites_fluc(:)
    REAL(rt), ALLOCATABLE :: occupation_sites_temp(:)
    COMPLEX(rtc) :: aux_sum
    COMPLEX(rtc),ALLOCATABLE :: Exponetial_mat(:,:)

    t=0.0_rt
    graduar=0.0_rt
    u=22
    u2=23
    u3=24
    u4=25
    psi0_dim=SIZE(psi0)

    ALLOCATE(psit(psi0_dim))
    ALLOCATE(psi_diag(psi0_dim))
    ALLOCATE(psi_temp(psi0_dim))
    ALLOCATE(Exponetial_mat(psi0_dim,psi0_dim))
    ALLOCATE(occupation_sites(sites))
    ALLOCATE(occup_sites_fluc(sites))
    ALLOCATE(occupation_sites_temp(sites))

    psit=DCmplx(0.0_rt,0.0_rt)
    psi_diag=0.0_rt
    psi_temp=0.0_rt
    aux_sum=DCmplx(0.0_rt,0.0_rt)

    OPEN(newunit=u, file="Survival_probablity.dat", status="replace")
    OPEN(newunit=u2, file="Number_P_Componets.dat", status="replace")
    OPEN(newunit=u3, file="Occupation_distribution.dat", status="replace")
    OPEN(newunit=u4, file="Occupation_distribution_fluctuations.dat", status="replace")


    DO jt=1,Psi0_dim
       DO kt=1,Psi0_dim
          psi_diag(jt)=psi_diag(jt)+((V(jt,kt)**(2.0_rt))&
               *(Dot_PRODUCT(V(:,kt),psi0)**(2.0_rt)))
       END DO
    END DO

    WRITE(*,*)"Infinity time average: ..."
    DO it=0,200
       t=(it*50.0_rt/200.0_rt)
       Exponetial_mat=Exp_mat(A,V,t0,t)
       DO kt=1,psi0_dim
          DO zt=1,Psi0_dim
             aux_sum=aux_sum+(Exponetial_mat(kt,zt)*psi0(zt))
          END DO
          psit(kt)=(ABS(aux_sum)**2.0_rt)

          IF(it.EQ.0)THEN
             psi_temp(kt)=psi_temp(kt)+(psit(kt)/(2.0_rt))
          ELSE IF(it.EQ.200)THEN
             psi_temp(kt)=psi_temp(kt)+(psit(kt)/(2.0_rt))
          ELSE
             psi_temp(kt)=psi_temp(kt)+psit(kt)
          END IF
          aux_sum=DCmplx(0.0_rt,0.0_rt)

       END DO
    END DO
    psi_temp=(psi_temp/200.0_rt)

    NPC_stationary=0.0_rt
    DO jt=1,psi0_dim
       NPC_stationary=NPC_stationary+(psi_temp(jt)**(2.0_rt))
    END DO


    WRITE(*,*)"Infinity time average: ... done"
    WRITE(*,*)"max\min:", MAXVAL(psi_temp),MINVAL(psi_temp),NPC_stationary

    NPC_stationary=NPC_stationary**(-1.0_rt)


    psit=DCmplx(0.0_rt,0.0_rt)
    aux_sum=DCmplx(0.0_rt,0.0_rt)


    DO it=0,partitions

       t=t0+(it*(tf-t0)/partitions)
       Exponetial_mat=Exp_mat(A,V,t0,t)

       DO kt=1,psi0_dim
          DO zt=1,Psi0_dim
             aux_sum=aux_sum+(Exponetial_mat(kt,zt)*psi0(zt))
          END DO
          psit(kt)=(ABS(aux_sum)**2.0_rt)
          aux_sum=DCmplx(0.0_rt,0.0_rt)
          !write(*,*) psit(kt),psi_diag(kt),psi_temp(kt)
       END DO
       WRITE(*,*)"=========================="

       NPC=0.0_rt

       DO jt=1,psi0_dim
          NPC=NPC+(psit(jt)**(2.0_rt))
       END DO

       NPC=NPC**(-1.0_rt)
       occupation_sites=0.0_rt

       DO kt=1,Sites
          DO jt=1,psi0_dim
             occupation_sites(kt)=occupation_sites(kt)+(NBase(jt,kt)*psit(jt))
          END DO
       END DO


       DO jt=1,psi0_dim
          WRITE(u,*)  t, jt, psit(jt)/MAXVAL(psit)
       END DO
       WRITE(u,*)"  "
       WRITE(u,*)"  "

       WRITE(u2,*) t, NPC!,NPC_stationary
       WRITE(u3,*) t, occupation_sites

       WRITE(*,*) it
    END DO

    occup_sites_fluc=0.0_rt
    occupation_sites_temp=0.0_rt

    DO kt=1,Sites
       DO jt=1,psi0_dim
          occupation_sites_temp(kt)=occupation_sites_temp(kt)+(NBase(jt,kt)*psi_temp(jt))
          occup_sites_fluc(kt)=occup_sites_fluc(kt)+((NBase(jt,kt)**(2.0_rt))*psi_temp(jt))
       END DO
    END DO

    DO jt=1,SITES
      graduar=(occup_sites_fluc(jt)-(occupation_sites_temp(jt)**(2.0_rt)))/(occupation_sites_temp(jt)**2)
       WRITE(u4,*) graduar-1.0_rt, (1.0_rt/occupation_sites_temp(jt)), occupation_sites_temp(jt)
    END DO

    CLOSE(u)
    CLOSE(u2)
    CLOSE(u3)
    CLOSE(u4)

    DEALLOCATE(Exponetial_mat)
    DEALLOCATE(occupation_sites)
    DEALLOCATE(occup_sites_fluc)
    DEALLOCATE(psit)
    DEALLOCATE(psi_diag)

  END SUBROUTINE survival_p


  SUBROUTINE easy_comulative(vec)

    USE nrtype
    USE combinatories
    USE many_particleP
    IMPLICIT NONE

    INTEGER :: it,u
    REAL(rt) :: sum
    REAL(rt), DIMENSION(:) :: vec

    CALL quicksort(vec)

    sum=0.0_rt
    u=1

    OPEN(newunit=u, file="comulative_dist_lea.dat", status="replace")

    DO it=1,SIZE(vec)
       sum=1.0_rt+sum
       WRITE(u,*) vec(it), sum/SIZE(vec)
    END DO

    CLOSE(u)

    sum=0.0_rt


  END SUBROUTINE easy_comulative


  REAL(rt) FUNCTION delta_rms(vec_exp,vec_theo) RESULT(del)

    USE nrtype
    IMPLICIT NONE

    REAL(rt), DIMENSION(:) :: vec_exp,vec_theo
    INTEGER :: it

    del=0.0_rt

    DO it=1,SIZE(vec_exp)
       del=del+((vec_exp(it)-vec_theo(it))**2)
    END DO

    del=del/SIZE(vec_exp)

  END FUNCTION delta_rms

END MODULE Thermalization_routines
