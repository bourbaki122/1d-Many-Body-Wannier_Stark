MODULE class_info

USE nrtype
INTEGER :: maxindex_a,maxindex_d,maxindex_y
INTEGER :: maxindex_bf,maxindex_xf,maxindex_cf
INTEGER :: maxindex_b,maxindex_x,maxindex_c
INTEGER :: maxindex_z,maxindex_zm, maxindex_zf
INTEGER, ALLOCATABLE :: Ma_infor(:,:,:),Md_infor(:,:,:),MY_infor(:,:,:)
INTEGER, ALLOCATABLE :: Mbf_infor(:,:,:),Mcf_infor(:,:,:),Mxf_infor(:,:,:)
INTEGER, ALLOCATABLE :: Mb_infor(:,:,:),Mc_infor(:,:,:),Mx_infor(:,:,:)
INTEGER, ALLOCATABLE :: Mz_infor(:,:,:),Mzm_infor(:,:,:),Mzf_infor(:,:,:)
!=================for complete two-body========================!
INTEGER(ITT) :: full_index
INTEGER, ALLOCATABLE :: Mcomplete_infor(:,:)


CONTAINS

!==================================================================
!Subroutine Complete_twobody
!
!in this subroutine y put the complete two body interaction....
!im going to create a matrix that has the form
! M(n,m,o,p,q,s) where n--|n>, m--|m>, and they are related by
!|m>---a^t_{o}a^t_{p}a_{q}a_{s}|n> beeing a^{t}_{s} the creation and
!anihilation operators
!===============================================================

SUBROUTINE Complete_twobody

  USE nrtype
  USE single_particleP
  USE many_particleP

  IMPLICIT NONE

  INTEGER :: it,jt,uy
  INTEGER :: o,p,q,s
  INTEGER,DIMENSION(SITES) :: Aux1
  INTEGER,ALLOCATABLE :: Mcomplete_inforaux(:,:)


  ALLOCATE(Mcomplete_inforaux(HDIMENSION*HDIMENSION*(SITES**4),6))

  Mcomplete_inforaux=0
  full_index=0
  uy=0

  DO it=1,HDIMENSION
  !  write(*,*)"el it"
  !  write(*,*) it

     DO o=1,SITES
        DO p=1,SITES
           DO q=1,SITES
              DO s=1,SITES

                 !write(*,*) "BLOQUE"
                 Aux1=NBase(it,:)
                 IF(Aux1(s).EQ.0)CYCLE
                 !write(*,*) Aux1(:)
                 Aux1(s)=Aux1(s)-1
                 !write(*,*) Aux1(:)
                 IF(Aux1(q).EQ.0)CYCLE
                 Aux1(q)=Aux1(q)-1
                 !write(*,*) Aux1(:)
                 Aux1(p)=Aux1(p)+1
                 !write(*,*) Aux1(:)
                 Aux1(o)=Aux1(o)+1
                 !write(*,*) Aux1(:)
                 !write(*,*) "BLOQUE"

                 IF(MAXVAL(Aux1(1:SITES)).GT.PARTICLES)CYCLE
                 IF(MINVAL(Aux1(1:SITES)).LT.0)CYCLE
                 uy=uy+1

                 !WRITE(*,*) Aux1(:), uy
                 !write(*,*)  full_index

                 DO jt=1,HDIMENSION

                    IF(ALL(Aux1(1:sites).EQ.NBase(jt,1:Sites)))THEN

                       full_index=full_index+1
                       !write(*,*) s,q,p,o, full_index

                       Mcomplete_inforaux(full_index,1)=it
                       Mcomplete_inforaux(full_index,2)=jt

                       Mcomplete_inforaux(full_index,3)=s
                       Mcomplete_inforaux(full_index,4)=q
                       Mcomplete_inforaux(full_index,5)=p
                       Mcomplete_inforaux(full_index,6)=o

                    END IF

                 END DO

              END DO
           END DO
        END DO
     END DO

  END DO



  ALLOCATE(Mcomplete_infor(full_index,6))
  DO it=1,full_index
     Mcomplete_infor(it,:)=Mcomplete_inforaux(it,:)
  END DO
  DEALLOCATE(Mcomplete_inforaux)


END SUBROUTINE Complete_twobody


SUBROUTINE Two_body(w)
  USE nrtype
  USE single_particleP
  USE many_particleP
  USE HamilMatrix

  IMPLICIT NONE

  INTEGER :: it,numbe,a,b,c,d
  REAL(rt) :: w,aux,before
  REAL(rt),EXTERNAL :: M_global
  INTEGER, DIMENSION(SITES) :: playground
  before=0.0_rt
  aux=0.0_rt
  Vij=0.0_rt

  playground=0

  DO it=1,full_index

    a=Mcomplete_infor(it,3)
    b=Mcomplete_infor(it,4)
    c=Mcomplete_infor(it,5)
    d=Mcomplete_infor(it,6)

    playground=NBase(Mcomplete_infor(it,1),:)
    numbe=playground(a)
    playground(a)=playground(a)-1
    numbe=numbe*playground(b)
    playground(b)=playground(b)-1
    numbe=numbe*(playground(c)+1)
    playground(c)=playground(c)+1
    numbe=numbe*(playground(d)+1)
    playground(d)=playground(d)+1
    aux=SQRT(REAL(numbe))
    !WRITE(*,*) M_global(a,b,c,d)
    aux=aux*M_global(a,b,c,d)
    before=Vij(Mcomplete_infor(it,2),Mcomplete_infor(it,1))
    Vij(Mcomplete_infor(it,2),Mcomplete_infor(it,1))=before+aux
  END DO

  Vij=Vij*(w/2.0_rt)

END SUBROUTINE Two_body

END MODULE class_info
