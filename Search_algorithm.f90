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



Integer function Smart_coordinates(Fock_b, particles, sites)

Integer(Itt) :: Fock_b(:)
INTEGER :: it,jt,kt
Integer, allocatable :: numbe_fact(:,:)

ALLOCATE(numbe_fact(1:Sites,0:particles))
numbe(:,0)=1
DO sn=1,n
   numbe(m,sn)=numbe(m,sn-1)*(m+sn-1)/sn
END DO
