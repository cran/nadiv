subroutine makea_ml(pedn, ped, T, F, Dii, Aout)
 implicit none
!define variables passed to subroutine
 integer, intent(in) :: pedn, ped(3, pedn)
 real, intent(inout) :: T(pedn, pedn), F(pedn), Dii(pedn), Aout(pedn, pedn)
!define variables derived in subroutine
 integer :: i, j, sj, dj
 real :: ai, D(pedn, pedn), L(pedn, pedn), tL(pedn, pedn)
 integer, allocatable :: ANi(:) 
!the program, following Meuwissen & Luo's 1992 algorithm 

 Dii = 1.0
 do i=1, pedn
   ai = 0.0
   allocate(ANi(1:i))
   ANi = 0
   ANi(i) = i
   T(i,i) = 1.0
   
   if (ped(3,i) /= -998) then
        Dii(i) = Dii(i) - 0.25*(1.0 + F(ped(3,i)))
   end if

   if (ped(2,i) /= -998) then 
        Dii(i) = Dii(i) - 0.25*(1.0 + F(ped(2,i)))
   end if


   do
        if (maxval(ANi) == 0) exit
        j = maxval(ANi)
        sj = ped(3,j)
        dj = ped(2,j)

        if (sj /= -998) then
	   ANi(sj) = sj
	   T(i, sj) =  T(i, sj) + 0.5*T(i, j)
        end if

        if (dj /= -998) then
	   ANi(dj) = dj
	   T(i, dj) =  T(i, dj) + 0.5*T(i,j)
        end if
        ai = ai + T(i, j)*T(i, j)*Dii(j)
        ANi(j) = 0
   end do
 
   if(allocated(ANi)) then
     deallocate(ANi)
   end if
   F(i) = ai - 1.0
 end do

end subroutine makea_ml

