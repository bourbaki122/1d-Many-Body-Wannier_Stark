MODULE Interactions
  use nrtype
  use many_particleP
  use single_particleP
  use HamilMatrix

  integer :: maxnumber_y,Maxnumber_a,Maxnumber_d
  INTEGER, ALLOCATABLE :: Ma_matches(:,:),Md_matches(:,:),My_matches(:,:)


  integer :: Maxnumber_b,Maxnumber_bf,Maxnumber_x,Maxnumber_xf
  integer :: Maxnumber_c,Maxnumber_cf
  INTEGER, ALLOCATABLE :: Mc_matches(:,:),Mcf_matches(:,:)
  INTEGER, ALLOCATABLE :: Mx_matches(:,:),Mxf_matches(:,:)
  INTEGER, ALLOCATABLE :: Mb_matches(:,:),Mbf_matches(:,:)

  INTEGER :: Maxnumber_z,Maxnumber_zm,Maxnumber_zf
  INTEGER, ALLOCATABLE :: Mz_matches(:,:),Mzm_matches(:,:),Mzf_matches(:,:)


  contains

   subroutine two_moves

     IMPLICIT none

     integer :: it,kt,zt
     integer :: delta
     INTEGER,DIMENSION(SITES) :: inicial_aux,mavec_aux,mdvec_aux,myvec_aux
     integer,ALLOCATABLE :: MA_temp_storage(:,:),MD_temp_storage(:,:)
     integer,ALLOCATABLE :: MY_temp_storage(:,:)

     ALLOCATE(MA_temp_storage(HDIMENSION*HDIMENSION*Sites,4))
     ALLOCATE(MD_temp_storage(HDIMENSION*HDIMENSION*Sites,4))
     ALLOCATE(MY_temp_storage(HDIMENSION*HDIMENSION*Sites,4))

     MA_temp_storage=0
     MD_temp_storage=0
     MY_temp_storage=0
     maxnumber_a=0
     maxnumber_d=0
     maxnumber_y=0

     do it=1,HDIMENSION

       inicial_aux=NBase(it,:)

       do kt=1,sites
         do delta=1,sites-kt

           mavec_aux=inicial_aux
           mdvec_aux=inicial_aux
           myvec_aux=inicial_aux

           IF(inicial_aux(kt).GT.0)THEN
              IF(inicial_aux(kt+delta).GT.0)THEN
                 mavec_aux(kt+delta)=inicial_aux(kt+delta)-1
                 mavec_aux(kt)=inicial_aux(kt)+1
              END IF
           END IF

           IF(inicial_aux(kt+delta).GT.0)THEN
              IF(inicial_aux(kt).GT.0)THEN
                 mdvec_aux(kt)=inicial_aux(kt)-1
                 mdvec_aux(kt+delta)=inicial_aux(kt+delta)+1
              END IF
           END IF

           IF(inicial_aux(kt+delta).GT.1)THEN
              myvec_aux(kt+delta)=inicial_aux(kt+delta)-2
              myvec_aux(kt)=inicial_aux(kt)+2
           END IF

           DO zt=1,HDIMENSION

              IF(ALL(mavec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

              ELSE
                 IF(ALL(mavec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                    maxnumber_a=maxnumber_a+1
                    MA_temp_storage(maxnumber_a,1)=kt
                    MA_temp_storage(maxnumber_a,2)=delta
                    MA_temp_storage(maxnumber_a,3)=zt
                    MA_temp_storage(maxnumber_a,4)=it
                 END IF
              END IF

              IF(ALL(mdvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

              ELSE
                 IF(ALL(mdvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                    maxnumber_d=maxnumber_d+1
                    MD_temp_storage(maxnumber_d,1)=kt
                    MD_temp_storage(maxnumber_d,2)=delta
                    MD_temp_storage(maxnumber_d,3)=zt
                    MD_temp_storage(maxnumber_d,4)=it
                 END IF
              END IF

              IF(ALL(myvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

              ELSE
                 IF(ALL(myvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                    maxnumber_y=maxnumber_y+1
                    MY_temp_storage(maxnumber_y,1)=kt
                    MY_temp_storage(maxnumber_y,2)=delta
                    MY_temp_storage(maxnumber_y,3)=zt
                    MY_temp_storage(maxnumber_y,4)=it
                 END IF
             END IF
           END DO

        end do
      end do
    end do

    ALLOCATE(MA_matches(1:maxnumber_a,1:4))
    ALLOCATE(MD_matches(1:maxnumber_d,1:4))
    ALLOCATE(MY_matches(1:maxnumber_y,1:4))
    DO it=1,maxnumber_a
       MA_matches(it,:)=MA_temp_storage(it,:)
    END DO
    DO it=1,maxnumber_d
       MD_matches(it,:)=MD_temp_storage(it,:)
    END DO
    DO it=1,maxnumber_y
       MY_matches(it,:)=MY_temp_storage(it,:)
    END DO

    DEALLOCATE(MA_temp_storage,MD_temp_storage,MY_temp_storage)

    write(*,*) "Two Moves: CHECK"

   end subroutine two_moves


  subroutine Three_moves

   IMPLICIT none

   integer :: it,kt,zt
   integer :: gamma,delta
   INTEGER,DIMENSION(SITES) :: inicial_aux,mcvec_aux,mbvec_aux,mxvec_aux
   INTEGER,DIMENSION(SITES) :: mcfvec_aux,mbfvec_aux,mxfvec_aux
   integer,ALLOCATABLE :: MC_temp_storage(:,:),MCF_temp_storage(:,:)
   integer,ALLOCATABLE :: MX_temp_storage(:,:),MXF_temp_storage(:,:)
   integer,ALLOCATABLE :: MB_temp_storage(:,:),MBF_temp_storage(:,:)

   ALLOCATE(MC_temp_storage(HDIMENSION*HDIMENSION*Sites,5))
   ALLOCATE(MX_temp_storage(HDIMENSION*HDIMENSION*Sites,5))
   ALLOCATE(MB_temp_storage(HDIMENSION*HDIMENSION*Sites,5))
   ALLOCATE(MCF_temp_storage(HDIMENSION*HDIMENSION*Sites,5))
   ALLOCATE(MBF_temp_storage(HDIMENSION*HDIMENSION*Sites,5))
   ALLOCATE(MXF_temp_storage(HDIMENSION*HDIMENSION*Sites,5))

   Maxnumber_b=0
   Maxnumber_bf=0
   Maxnumber_x=0
   Maxnumber_xf=0
   Maxnumber_c=0
   Maxnumber_cf=0

   MC_temp_storage=0
   MCF_temp_storage=0
   MB_temp_storage=0
   MBF_temp_storage=0
   MX_temp_storage=0
   MXF_temp_storage=0


   do it=1,HDIMENSION

     inicial_aux=NBase(it,:)

     do kt=1,SITES

      do delta=1,sites-kt

       do gamma=1,sites-delta-kt

         mcvec_aux=inicial_aux
         mcfvec_aux=inicial_aux
         mbvec_aux=inicial_aux
         mbfvec_aux=inicial_aux
         mxvec_aux=inicial_aux
         mxfvec_aux=inicial_aux

         IF(inicial_aux(kt).GT.0)THEN
            IF(inicial_aux(kt+delta).GT.0)THEN
               mbfvec_aux(kt+delta)=inicial_aux(kt+delta)-1
               mbfvec_aux(kt+delta+gamma)=inicial_aux(kt+delta+gamma)+1
            END IF
         END IF

         IF(inicial_aux(kt+delta).GT.0)THEN
            IF(inicial_aux(kt).GT.0)THEN
               mxfvec_aux(kt)=inicial_aux(kt)-1
               mxfvec_aux(kt+delta+gamma)=inicial_aux(kt+delta+gamma)+1
            END IF
         END IF

         IF(inicial_aux(kt+delta+gamma).GT.0)THEN
            IF(inicial_aux(kt).GT.0)THEN
               mcfvec_aux(kt)=inicial_aux(kt)-1
               mcfvec_aux(kt+delta)=inicial_aux(kt+delta)+1
            END IF
         END IF

         IF(inicial_aux(kt).GE.2)THEN
            mbvec_aux(kt)=inicial_aux(kt)-2
            mbvec_aux(kt+ delta)=inicial_aux(kt+ delta)+1
            mbvec_aux(kt+ delta+gamma)=inicial_aux(kt+delta+gamma)+1
         END IF

         IF(inicial_aux(kt+delta).GE.2)THEN
            mxvec_aux(kt+delta)=inicial_aux(kt+delta)-2
            mxvec_aux(kt)=inicial_aux(kt)+1
            mxvec_aux(kt+delta+gamma)=inicial_aux(kt+delta+gamma)+1
         END IF

         IF(inicial_aux(kt+delta+gamma).GE.2)THEN
            mcvec_aux(kt+delta+gamma)=inicial_aux(kt+delta+gamma)-2
            mcvec_aux(kt+delta)=inicial_aux(kt+delta)+1
            mcvec_aux(kt)=inicial_aux(kt)+1
         END IF


         DO zt=1,HDIMENSION

            IF(ALL(mbfvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

            ELSE
               IF(ALL(mbfvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                  Maxnumber_bf=Maxnumber_bf+1
                  MBF_temp_storage(Maxnumber_bf,1)=kt
                  MBF_temp_storage(Maxnumber_bf,2)=delta
                  MBF_temp_storage(Maxnumber_bf,3)=gamma
                  MBF_temp_storage(Maxnumber_bf,4)=zt
                  MBF_temp_storage(Maxnumber_bf,5)=it
               END IF
            END IF

            IF(ALL(mxfvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

            ELSE
               IF(ALL(mxfvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                  Maxnumber_xf=Maxnumber_xf+1
                  MXF_temp_storage(Maxnumber_xf,1)=kt
                  MXF_temp_storage(Maxnumber_xf,2)=delta
                  MXF_temp_storage(Maxnumber_xf,3)=gamma
                  MXF_temp_storage(Maxnumber_xf,4)=zt
                  MXF_temp_storage(Maxnumber_xf,5)=it
                END IF
            END IF

            IF(ALL(mcfvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

            ELSE
               IF(ALL(mcfvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                  Maxnumber_cf=Maxnumber_cf+1
                  MCF_temp_storage(Maxnumber_cf,1)=kt
                  MCF_temp_storage(Maxnumber_cf,2)=delta
                  MCF_temp_storage(Maxnumber_cf,3)=gamma
                  MCF_temp_storage(Maxnumber_cf,4)=zt
                  MCF_temp_storage(Maxnumber_cf,5)=it
               END IF
            END IF

            IF(ALL(mbvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

            ELSE
               IF(ALL(mbvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                  Maxnumber_b=Maxnumber_b+1
                  MB_temp_storage(Maxnumber_b,1)=kt
                  MB_temp_storage(Maxnumber_b,2)=delta
                  MB_temp_storage(Maxnumber_b,3)=gamma
                  MB_temp_storage(Maxnumber_b,4)=zt
                  MB_temp_storage(Maxnumber_b,5)=it
               END IF
            END IF

            IF(ALL(mxvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

            ELSE
               IF(ALL(mxvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                  Maxnumber_x=Maxnumber_x+1
                  MX_temp_storage(Maxnumber_x,1)=kt
                  MX_temp_storage(Maxnumber_x,2)=delta
                  MX_temp_storage(Maxnumber_x,3)=gamma
                  MX_temp_storage(Maxnumber_x,4)=zt
                  MX_temp_storage(Maxnumber_x,5)=it
               END IF
            END IF

            IF(ALL(mcvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

            ELSE
               IF(ALL(mcvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                  Maxnumber_c=Maxnumber_c+1
                  MC_temp_storage(Maxnumber_c,1)=kt
                  MC_temp_storage(Maxnumber_c,2)=delta
                  MC_temp_storage(Maxnumber_c,3)=gamma
                  MC_temp_storage(Maxnumber_c,4)=zt
                  MC_temp_storage(Maxnumber_c,5)=it
               END IF
            END IF

        end do
      end do
    end do
  end do
end do

ALLOCATE(Mbf_matches(1:Maxnumber_bf,1:5))
ALLOCATE(Mxf_matches(1:Maxnumber_xf,1:5))
ALLOCATE(Mcf_matches(1:Maxnumber_cf,1:5))
ALLOCATE(Mb_matches(1:Maxnumber_b,1:5))
ALLOCATE(Mx_matches(1:Maxnumber_x,1:5))
ALLOCATE(Mc_matches(1:Maxnumber_c,1:5))

DO it=1,Maxnumber_bf
   Mbf_matches(it,:)=MBF_temp_storage(it,:)
END DO
DO it=1,Maxnumber_xf
   Mxf_matches(it,:)=MXF_temp_storage(it,:)
END DO
DO it=1,Maxnumber_cf
   Mcf_matches(it,:)=MCF_temp_storage(it,:)
END DO
DO it=1,Maxnumber_b
   Mb_matches(it,:)=MB_temp_storage(it,:)
END DO
DO it=1,Maxnumber_x
   Mx_matches(it,:)=MX_temp_storage(it,:)
END DO
DO it=1,Maxnumber_c
   Mc_matches(it,:)=MC_temp_storage(it,:)
END DO

DEALLOCATE(MBF_temp_storage,MXF_temp_storage,MCF_temp_storage)
DEALLOCATE(MB_temp_storage,MC_temp_storage,MX_temp_storage)


write(*,*) "Three Moves: CHECK"

end subroutine Three_moves


subroutine four_moves
   implicit none

   integer :: it,kt,zt
   integer :: gamma,delta,rho
   INTEGER,DIMENSION(SITES) :: inicial_aux,mzvec_aux,mzmvec_aux,mzfvec_aux
   integer,ALLOCATABLE :: MZ_temp_storage(:,:),MZM_temp_storage(:,:)
   integer,ALLOCATABLE :: MZF_temp_storage(:,:)

   ALLOCATE(MZ_temp_storage(HDIMENSION*HDIMENSION*Sites,6))
   ALLOCATE(MZM_temp_storage(HDIMENSION*HDIMENSION*Sites,6))
   ALLOCATE(MZF_temp_storage(HDIMENSION*HDIMENSION*Sites,6))

   Maxnumber_z=0
   Maxnumber_zm=0
   Maxnumber_zf=0

   MZ_temp_storage=0
   MZm_temp_storage=0
   MZf_temp_storage=0

   do it=1,HDIMENSION

     inicial_aux=NBase(it,:)

     do kt=1,SITES
       do delta=1,sites-kt
         do gamma=1,sites-kt-delta
           do rho=1,sites-kt-delta-gamma

            mzvec_aux=inicial_aux
            mzmvec_aux=inicial_aux
            mzfvec_aux=inicial_aux

            IF(Mzvec_aux(kt).GT.0)THEN
               IF(Mzvec_aux(kt+delta).GT.0)THEN
                  Mzvec_aux(kt)=inicial_aux(kt)-1
                  Mzvec_aux(kt+(delta))=inicial_aux(kt+(delta))-1
                  Mzvec_aux(kt+(delta+gamma+rho))= &
                  inicial_aux(kt+(delta+gamma+rho))+1
                  Mzvec_aux(kt+(delta+gamma))=inicial_aux(kt+(delta+gamma))+1
               END IF
            END IF

            IF(Mzmvec_aux(kt).GT.0)THEN
               IF(Mzmvec_aux(kt+(delta+gamma)).GT.0)THEN
                  Mzmvec_aux(kt)=inicial_aux(kt)-1
                  Mzmvec_aux(kt+(delta+gamma))=inicial_aux(kt+(delta+gamma))-1
                  Mzmvec_aux(kt+(delta+gamma+rho))= &
                  inicial_aux(kt+(delta+gamma+rho))+1
                  Mzmvec_aux(kt+(delta))=inicial_aux(kt+(delta))+1
               END IF
            END IF

            IF(Mzfvec_aux(kt+delta).GT.0)THEN
               IF(Mzfvec_aux(kt+(delta+gamma)).GT.0)THEN
                  Mzfvec_aux(kt+delta)=inicial_aux(kt+delta)-1
                  Mzfvec_aux(kt+(delta+gamma))=inicial_aux(kt+(delta+gamma))-1
                  Mzfvec_aux(kt+(delta+gamma+rho))= &
                  inicial_aux(kt+(delta+gamma+rho))+1
                  Mzfvec_aux(kt)=inicial_aux(kt)+1
               END IF
            END IF

            DO zt=1,HDIMENSION

               IF(ALL(mzvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

               ELSE
                  IF(ALL(mzvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                     Maxnumber_z=Maxnumber_z+1
                     MZ_temp_storage(Maxnumber_z,1)=kt
                     MZ_temp_storage(Maxnumber_z,2)=delta
                     MZ_temp_storage(Maxnumber_z,3)=GAMMA
                     MZ_temp_storage(Maxnumber_z,4)=rho
                     MZ_temp_storage(Maxnumber_z,5)=zt
                     MZ_temp_storage(Maxnumber_z,6)=it
                  END IF
               END IF


               IF(ALL(mzmvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

               ELSE
                  IF(ALL(mzmvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                     Maxnumber_zm=Maxnumber_zm+1
                     MZM_temp_storage(Maxnumber_zm,1)=kt
                     MZM_temp_storage(Maxnumber_zm,2)=delta
                     MZM_temp_storage(Maxnumber_zm,3)=Gamma
                     MZM_temp_storage(Maxnumber_zm,4)=rho
                     MZM_temp_storage(Maxnumber_zm,5)=zt
                     MZM_temp_storage(Maxnumber_zm,6)=it
                  END IF
               END IF

               IF(ALL(mzfvec_aux(1:SITES).EQ.NBase(it,1:SITES)))THEN

               ELSE
                  IF(ALL(mzfvec_aux(1:SITES).EQ.NBase(zt,1:SITES)))THEN
                     Maxnumber_zf=Maxnumber_zf+1
                     MzF_temp_storage(Maxnumber_zf,1)=kt
                     MzF_temp_storage(Maxnumber_zf,2)=delta
                     MzF_temp_storage(Maxnumber_zf,3)=GAMMA
                     MzF_temp_storage(Maxnumber_zf,4)=rho
                     MzF_temp_storage(Maxnumber_zf,5)=zt
                     MzF_temp_storage(Maxnumber_zf,6)=it
                  END IF
               END IF

            END DO

       end do
     end do
   end do
 end do
end do

ALLOCATE(Mz_matches(1:Maxnumber_z,1:6))
ALLOCATE(Mzm_matches(1:Maxnumber_zm,1:6))
ALLOCATE(Mzf_matches(1:Maxnumber_zf,1:6))

DO it=1,Maxnumber_z
   Mz_matches(it,:)=Mz_temp_storage(it,:)
END DO
DO it=1,Maxnumber_zm
   Mzm_matches(it,:)=MZm_temp_storage(it,:)
END DO
DO it=1,Maxnumber_zf
   Mzf_matches(it,:)=MzF_temp_storage(it,:)
END DO

DEALLOCATE(MZ_temp_storage,MZm_temp_storage,MzF_temp_storage)

write(*,*) "Four Moves: CHECK"

end subroutine four_moves


!======================================================
!Here and further are the subroutines matrix !
!======================================================

SUBROUTINE MA_matrix_smart(w,delta)

  use single_particleP
  IMPLICIT NONE

  INTEGER :: it,delta
  REAL(rt) :: w,numbe,aux
  REAL(rt),EXTERNAL :: M_local,M_global

  Ma_trix=0.0_rt

  DO it=1,maxnumber_a
     IF(Ma_matches(it,2).eq.delta)THEN
        numbe=NBase(Ma_matches(it,4),ma_matches(it,1))
        numbe=numbe*SQRT(REAL(NBase(Ma_matches(it,4),ma_matches(it,1)+delta)))
        numbe=numbe*SQRT(REAL(NBase(Ma_matches(it,4),Ma_matches(it,1))+1))
        numbe=numbe*M_global(ma_matches(it,1),ma_matches(it,1), &
        ma_matches(it,1),ma_matches(it,1)+delta)*2.0_rt
        !write(*,*) numbe
        aux=MA_trix(Ma_matches(it,3),Ma_matches(it,4))
        MA_trix(Ma_matches(it,3),Ma_matches(it,4))=numbe+aux
     END IF
  END DO

  Ma_trix=(Ma_trix+TRANSPOSE(Ma_trix))*(w/2.0_RT)
END SUBROUTINE MA_matrix_smart


SUBROUTINE MD_matrix_smart(w,delta)

use single_particleP
IMPLICIT NONE

INTEGER :: it,delta
REAL(rt) :: w,numbe,aux
REAL(rt),EXTERNAL :: M_local,M_global

Md_trix=0.0_rt

DO it=1,maxnumber_d
 IF(MD_matches(it,2).eq.delta)THEN
  numbe=NBase(MD_matches(it,4),MD_matches(it,1)+delta)
  numbe=numbe*SQRT(REAL(NBase(MD_matches(it,4),MD_matches(it,1))))
  numbe=numbe*SQRT(REAL(NBase(MD_matches(it,4),MD_matches(it,1)+delta)+1))
  numbe=numbe*M_global(MD_matches(it,1)+delta,MD_matches(it,1)+delta, &
  MD_matches(it,1)+delta,MD_matches(it,1))*2.0_rt
  !write(*,*) numbe
  aux=Md_trix(MD_matches(it,3),MD_matches(it,4))
  Md_trix(MD_matches(it,3),MD_matches(it,4))=numbe+aux
 END IF
END DO

Md_trix=(Md_trix+TRANSPOSE(Md_trix))*(w/2.0_RT)
END SUBROUTINE MD_matrix_smart


SUBROUTINE MY_matrix_smart(w,delta)

USE single_particleP

IMPLICIT NONE

INTEGER :: it,delta
REAL(rt) :: w,numbe,aux
REAL(rt),EXTERNAL :: M_local,M_global

MY_trix=0.0_rt

DO it=1,maxnumber_y
 IF(MY_matches(it,2).eq.delta)THEN
  numbe=NBase(MY_matches(it,4),MY_matches(it,1)+delta)
  numbe=SQRT(REAL(numbe*(numbe-1)))
  numbe=numbe*SQRT(REAL(NBase(MY_matches(it,4),MY_matches(it,1))+1))
  numbe=numbe*SQRT(REAL(NBase(MY_matches(it,4),MY_matches(it,1))+2))
  numbe=numbe*M_global(MY_matches(it,1),MY_matches(it,1),&
  MY_matches(it,1)+delta,MY_matches(it,1)+delta)
  aux=My_trix(MY_matches(it,3),MY_matches(it,4))
  My_trix(MY_matches(it,3),MY_matches(it,4))=numbe+aux
 END IF
END DO

My_trix=(My_trix+TRANSPOSE(My_trix))*(w/2.0_rt)
END SUBROUTINE MY_matrix_smart



SUBROUTINE MB_matrix_smart(w,delta,gamma)

USE nrtype
USE single_particleP

IMPLICIT NONE

INTEGER :: it,delta,gamma
REAL(rt) :: w,numbe,aux
REAL(rt),EXTERNAL :: M_local,M_global

Mb_trix=0.0_rt

DO it=1,Maxnumber_b
 IF((Mb_matches(it,2).eq.delta).and.(Mb_matches(it,3).eq.gamma))THEN
  numbe=NBase(Mb_matches(it,5),Mb_matches(it,1))
  numbe=SQRT(numbe*(numbe-1))
  numbe=numbe*SQRT(REAL(NBase(Mb_matches(it,5),Mb_matches(it,1)+delta)+1))
  numbe=numbe*SQRT(REAL(NBase(Mb_matches(it,5),Mb_matches(it,1) &
  +delta+gamma)+1))
  numbe=numbe*M_global(Mb_matches(it,1),Mb_matches(it,1),&
  Mb_matches(it,1)+delta,Mb_matches(it,1)+delta+gamma)*2.0_rt
  !write(*,*) numbe
  aux=Mb_trix(Mb_matches(it,4),Mb_matches(it,5))
  Mb_trix(Mb_matches(it,4),Mb_matches(it,5))=numbe+aux
 END IF
END DO

DO it=1,Maxnumber_bf
 IF((Mbf_matches(it,2).eq.delta).and.(Mbf_matches(it,3).eq.gamma))THEN
  numbe=NBase(Mbf_matches(it,5),Mbf_matches(it,1))
  numbe=numbe*SQRT(REAL(NBase(Mbf_matches(it,5),Mbf_matches(it,1)+delta)))
  numbe=numbe*SQRT(REAL(NBase(Mbf_matches(it,5),Mbf_matches(it,1) &
  +delta+gamma)+1))
  numbe=numbe*M_global(Mbf_matches(it,1),Mbf_matches(it,1),&
  Mbf_matches(it,1)+delta,Mbf_matches(it,1)+delta+gamma)*4.0_RT
  !write(*,*) numbe
  aux=Mb_trix(Mbf_matches(it,4),Mbf_matches(it,5))
  Mb_trix(Mbf_matches(it,4),Mbf_matches(it,5))=aux+numbe
 END IF
END DO

Mb_trix=(Mb_trix+TRANSPOSE(Mb_trix))*(w/2.0_RT)
END SUBROUTINE MB_matrix_smart


SUBROUTINE Mx_matrix_smart(w,delta,gamma)

USE nrtype
USE single_particleP

IMPLICIT NONE

INTEGER :: it,delta,gamma
REAL(rt) :: w,numbe,aux
REAL(rt),EXTERNAL :: M_local,M_global

Mx_trix=0.0_rt

DO it=1,Maxnumber_x
 IF((Mx_matches(it,2).eq.delta).and.(Mx_matches(it,3).eq.gamma))THEN
  numbe=NBase(Mx_matches(it,5),Mx_matches(it,1)+delta)
  numbe=SQRT(numbe*(numbe-1))
  numbe=numbe*SQRT(REAL(NBase(Mx_matches(it,5),Mx_matches(it,1))+1))
  numbe=numbe*SQRT(REAL(NBase(Mx_matches(it,5),Mx_matches(it,1) &
  + delta+gamma)+1))
  numbe=numbe*M_global(Mx_matches(it,1)+delta,Mx_matches(it,1)+delta, &
  Mx_matches(it,1),Mx_matches(it,1)+delta+gamma)*2.0_rt
  !write(*,*) numbe
  aux=Mx_trix(Mx_matches(it,4),Mx_matches(it,5))
  Mx_trix(Mx_matches(it,4),Mx_matches(it,5))=aux+numbe
 END IF
END DO

DO it=1,maxnumber_xf
 IF((Mxf_matches(it,2).eq.delta).and.(Mxf_matches(it,3).eq.gamma))THEN
  numbe=NBase(Mxf_matches(it,5),Mxf_matches(it,1)+delta)
  numbe=numbe*SQRT(REAL(NBase(Mxf_matches(it,5),Mxf_matches(it,1))))
  numbe=numbe*SQRT(REAL(NBase(Mxf_matches(it,5),Mxf_matches(it,1) &
  + delta+gamma)+1))
  numbe=numbe*M_global(Mxf_matches(it,1)+delta,Mxf_matches(it,1)+delta, &
  Mxf_matches(it,1),Mxf_matches(it,1)+delta+gamma)*4.0_rt
  !write(*,*) numbe
  aux=Mx_trix(Mxf_matches(it,4),Mxf_matches(it,5))
  Mx_trix(Mxf_matches(it,4),Mxf_matches(it,5))=aux+numbe
 END IF
END DO

Mx_trix=(Mx_trix+TRANSPOSE(Mx_trix))*(W/2.0_RT)
END SUBROUTINE Mx_matrix_smart

SUBROUTINE Mc_matrix_smart(w,delta,gamma)

USE nrtype
USE single_particleP

IMPLICIT NONE

INTEGER :: it,delta,gamma
REAL(rt) :: w,numbe,aux
REAL(rt),EXTERNAL :: M_local,M_global

Mc_trix=0.0_rt

DO it=1,maxnumber_c
 IF((Mc_matches(it,2).eq.delta).and.(Mc_matches(it,3).eq.gamma))THEN
  numbe=NBase(Mc_matches(it,5),Mc_matches(it,1)+delta+gamma)
  numbe=SQRT(numbe*(numbe-1))
  numbe=numbe*SQRT(REAL(NBase(Mc_matches(it,5),Mc_matches(it,1))+1))
  numbe=numbe*SQRT(REAL(NBase(Mc_matches(it,5),Mc_matches(it,1)+delta)+1))
  numbe=numbe*M_global(Mc_matches(it,1),Mc_matches(it,1)+delta, &
  Mc_matches(it,1)+delta+gamma,Mc_matches(it,1)+delta+gamma)*2.0_rt
  !write(*,*) numbe
  aux=Mc_trix(Mc_matches(it,4),Mc_matches(it,5))
  Mc_trix(Mc_matches(it,4),Mc_matches(it,5))=numbe+aux
 END IF
END DO


DO it=1,Maxnumber_cf
 IF((Mcf_matches(it,2).eq.delta).and.(Mcf_matches(it,3).eq.gamma))THEN
  numbe=NBase(Mcf_matches(it,5),Mcf_matches(it,1)+delta+gamma)
  numbe=numbe*SQRT(REAL(NBase(Mcf_matches(it,5),Mcf_matches(it,1))))
  numbe=numbe*SQRT(REAL(NBase(Mcf_matches(it,5),Mcf_matches(it,1) &
  +delta)+1))
  numbe=numbe*M_global(Mcf_matches(it,1),Mcf_matches(it,1)+delta, &
  Mcf_matches(it,1)+delta+gamma,Mcf_matches(it,1)+delta+gamma)*4.0_rt
  !write(*,*) numbe
  aux=Mc_trix(Mcf_matches(it,4),Mcf_matches(it,5))
  Mc_trix(Mcf_matches(it,4),Mcf_matches(it,5))=numbe+aux
 END IF
END DO

Mc_trix=(Mc_trix+TRANSPOSE(Mc_trix))*(w/2.0_RT)
END SUBROUTINE Mc_matrix_smart


SUBROUTINE Mz_matrix_smart(w,delta,gamma,rho)

USE nrtype
USE single_particleP
USE many_particleP
USE HamilMatrix

IMPLICIT NONE

INTEGER :: it,delta,gamma,rho
REAL(rt) :: w,numbe,aux
REAL(rt),EXTERNAL :: M_local, M_global

Mz_trix=0.0_rt

DO it=1,Maxnumber_z
 IF((Mz_matches(it,2).eq.delta).and.(Mz_matches(it,3).eq.gamma) &
 .and.(Mz_matches(it,4).eq.rho))THEN
  numbe=NBase(Mz_matches(it,6),Mz_matches(it,1))
  numbe=numbe*NBase(Mz_matches(it,6),Mz_matches(it,1)+delta)
  numbe=SQRT(REAL(numbe))
  numbe=numbe*SQRT(REAL(NBase(Mz_matches(it,6),Mz_matches(it,1) &
  +delta+gamma)+1))
  numbe=numbe*SQRT(REAL(NBase(Mz_matches(it,6),Mz_matches(it,1) &
  +delta+gamma+rho)+1))
  numbe=numbe*M_global(Mz_matches(it,1),Mz_matches(it,1)+delta, &
  Mz_matches(it,1)+delta+gamma,Mz_matches(it,1)+delta+gamma+rho)*4.0_rt
  !write(*,*) numbe
  aux=Mz_trix(Mz_matches(it,5),Mz_matches(it,6))
  Mz_trix(Mz_matches(it,5),Mz_matches(it,6))=numbe+aux
 END IF
END DO

DO it=1,Maxnumber_zm
 IF((Mzm_matches(it,2).eq.delta).and.(Mzm_matches(it,3).eq.gamma) &
 .and.(Mzm_matches(it,4).eq.rho))THEN
  numbe=NBase(Mzm_matches(it,6),Mzm_matches(it,1))
  numbe=numbe*NBase(Mzm_matches(it,6),Mzm_matches(it,1)+delta+gamma)
  numbe=SQRT(REAL(numbe))
  numbe=numbe*SQRT(REAL(NBase(Mzm_matches(it,6),Mzm_matches(it,1) &
  +delta)+1))
  numbe=numbe*SQRT(REAL(NBase(Mzm_matches(it,6),Mzm_matches(it,1) &
  +delta+gamma+rho)+1))
  numbe=numbe*M_global(Mzm_matches(it,1),Mzm_matches(it,1)+delta, &
  Mzm_matches(it,1)+delta+gamma,Mzm_matches(it,1)+delta+gamma+rho)*4.0_rt
  !write(*,*) numbe
  aux=Mz_trix(Mzm_matches(it,5),Mzm_matches(it,6))
  Mz_trix(Mzm_matches(it,5),Mzm_matches(it,6))=numbe+aux
 END IF
END DO

DO it=1,maxnumber_zf
 IF((Mzf_matches(it,2).eq.delta).and.(Mzf_matches(it,3).eq.gamma) &
 .and.(Mzf_matches(it,4).eq.rho))THEN
  numbe=NBase(Mzf_matches(it,6),Mzf_matches(it,1)+delta)
  numbe=numbe*NBase(Mzf_matches(it,6),Mzf_matches(it,1)+delta+gamma)
  numbe=SQRT(REAL(numbe))
  numbe=numbe*SQRT(REAL(NBase(Mzf_matches(it,6),Mzf_matches(it,1))+1))
  numbe=numbe*SQRT(REAL(NBase(Mzf_matches(it,6),Mzf_matches(it,1) &
  +delta+gamma+rho)+1))
  numbe=numbe*M_global(Mzf_matches(it,1),Mzf_matches(it,1)+delta, &
  Mzf_matches(it,1)+delta+gamma,Mzf_matches(it,1)+delta+gamma+rho)*4.0_rt
  aux=Mz_trix(Mzf_matches(it,5),Mzf_matches(it,6))
  Mz_trix(Mzf_matches(it,5),Mzf_matches(it,6))=numbe+aux
 END IF
END DO

Mz_trix=(Mz_trix+TRANSPOSE(Mz_trix))*(w/2.0_rt)
END SUBROUTINE Mz_matrix_smart


End MODULE Interactions
