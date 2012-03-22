subroutine a_ml(N, ped, Tn, T, Dii)
 implicit none
!define variables passed to subroutine
 integer, intent(in) :: N, ped(3, N), Tn
 real, intent(inout) :: T(Tn), Dii(N)
!define variables derived in subroutine
 integer :: i, j, k, sj, dj, ind
 real :: ai, F(N)
 integer, allocatable :: ANi(:) 
 type row
    real, pointer :: R(:)
 end type
 type(row) :: Ttmp(N)

 do j=1,N
    allocate (Ttmp(j)%R(1:j))
    Ttmp(j)%R(1:j)=0.0
 end do



!the program, following Meuwissen & Luo's 1992 algorithm 

 Dii = 1.0
 do i=1, N
   ai = 0.0
   allocate(ANi(1:i))
   ANi = 0
   ANi(i) = i
   Ttmp(i)%R(i) = 1.0
   
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
	   Ttmp(i)%R(sj) =  Ttmp(i)%R(sj) + 0.5*Ttmp(i)%R(j)
        end if

        if (dj /= -998) then
	   ANi(dj) = dj
	   Ttmp(i)%R(dj) =  Ttmp(i)%R(dj) + 0.5*Ttmp(i)%R(j)
        end if
        ai = ai + Ttmp(i)%R(j)*Ttmp(i)%R(j)*Dii(j)
        ANi(j) = 0
   end do
 
   if(allocated(ANi)) then
     deallocate(ANi)
   end if
   F(i) = ai - 1.0
 end do

 T(1)=Ttmp(1)%R(1)
 ind=1
 do k=2, N
   T=(/T(1:ind), Ttmp(k)%R(1:k)/)
   ind=ind+k
 end do

end subroutine a_ml

