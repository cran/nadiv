subroutine dij(iparents, iparents_r, Adense, dim_a, answer)
 implicit none
!define inputs
 integer, intent(in) :: iparents(4, iparents_r), iparents_r, dim_a
 real, intent(in) :: Adense(dim_a, dim_a)
 real, intent(out) :: answer(iparents_r)  
!define variables derived in subroutine
 integer :: i
 real, dimension(1:20) :: rmmp, rffp, rmfp, rfmp, d

!now the program that computes Dij for each j individual related to individual i
 do i=1,iparents_r, 20
     rmmp=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
     rffp=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
     rmfp=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
     rfmp=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
!11111
   rmmp=Adense(iparents(1,i), iparents(3,i))
   rffp=Adense(iparents(2,i), iparents(4,i))
   rmfp=Adense(iparents(1,i), iparents(4,i))
   rfmp=Adense(iparents(2,i), iparents(3,i))
!222222
   rmmp=(/rmmp(1), Adense(iparents(1,i+1), iparents(3,i+1)), rmmp(3:20)/)
   rffp=(/rffp(1), Adense(iparents(2,i+1), iparents(4,i+1)), rffp(3:20)/)
   rmfp=(/rmfp(1), Adense(iparents(1,i+1), iparents(4,i+1)), rmfp(3:20)/)
   rfmp=(/rfmp(1), Adense(iparents(2,i+1), iparents(3,i+1)), rfmp(3:20)/)
!333333
   rmmp=(/rmmp(1:2), Adense(iparents(1,i+2), iparents(3,i+2)), rmmp(4:20)/)
   rffp=(/rffp(1:2), Adense(iparents(2,i+2), iparents(4,i+2)), rffp(4:20)/)
   rmfp=(/rmfp(1:2), Adense(iparents(1,i+2), iparents(4,i+2)), rmfp(4:20)/)
   rfmp=(/rfmp(1:2), Adense(iparents(2,i+2), iparents(3,i+2)), rfmp(4:20)/)
!444444
   rmmp=(/rmmp(1:3), Adense(iparents(1,i+3), iparents(3,i+3)), rmmp(5:20)/)
   rffp=(/rffp(1:3), Adense(iparents(2,i+3), iparents(4,i+3)), rffp(5:20)/)
   rmfp=(/rmfp(1:3), Adense(iparents(1,i+3), iparents(4,i+3)), rmfp(5:20)/)
   rfmp=(/rfmp(1:3), Adense(iparents(2,i+3), iparents(3,i+3)), rfmp(5:20)/)
!555555
   rmmp=(/rmmp(1:4), Adense(iparents(1,i+4), iparents(3,i+4)), rmmp(6:20)/)
   rffp=(/rffp(1:4), Adense(iparents(2,i+4), iparents(4,i+4)), rffp(6:20)/)
   rmfp=(/rmfp(1:4), Adense(iparents(1,i+4), iparents(4,i+4)), rmfp(6:20)/)
   rfmp=(/rfmp(1:4), Adense(iparents(2,i+4), iparents(3,i+4)), rfmp(6:20)/)
!666666
   rmmp=(/rmmp(1:5), Adense(iparents(1,i+5), iparents(3,i+5)), rmmp(7:20)/)
   rffp=(/rffp(1:5), Adense(iparents(2,i+5), iparents(4,i+5)), rffp(7:20)/)
   rmfp=(/rmfp(1:5), Adense(iparents(1,i+5), iparents(4,i+5)), rmfp(7:20)/)
   rfmp=(/rfmp(1:5), Adense(iparents(2,i+5), iparents(3,i+5)), rfmp(7:20)/)
!777777
   rmmp=(/rmmp(1:6), Adense(iparents(1,i+6), iparents(3,i+6)), rmmp(8:20)/)
   rffp=(/rffp(1:6), Adense(iparents(2,i+6), iparents(4,i+6)), rffp(8:20)/)
   rmfp=(/rmfp(1:6), Adense(iparents(1,i+6), iparents(4,i+6)), rmfp(8:20)/)
   rfmp=(/rfmp(1:6), Adense(iparents(2,i+6), iparents(3,i+6)), rfmp(8:20)/)
!888888
   rmmp=(/rmmp(1:7), Adense(iparents(1,i+7), iparents(3,i+7)), rmmp(9:20)/)
   rffp=(/rffp(1:7), Adense(iparents(2,i+7), iparents(4,i+7)), rffp(9:20)/)
   rmfp=(/rmfp(1:7), Adense(iparents(1,i+7), iparents(4,i+7)), rmfp(9:20)/)
   rfmp=(/rfmp(1:7), Adense(iparents(2,i+7), iparents(3,i+7)), rfmp(9:20)/)
!999999
   rmmp=(/rmmp(1:8), Adense(iparents(1,i+8), iparents(3,i+8)), rmmp(10:20)/)
   rffp=(/rffp(1:8), Adense(iparents(2,i+8), iparents(4,i+8)), rffp(10:20)/)
   rmfp=(/rmfp(1:8), Adense(iparents(1,i+8), iparents(4,i+8)), rmfp(10:20)/)
   rfmp=(/rfmp(1:8), Adense(iparents(2,i+8), iparents(3,i+8)), rfmp(10:20)/)
!10
   rmmp=(/rmmp(1:9), Adense(iparents(1,i+9), iparents(3,i+9)), rmmp(11:20)/)
   rffp=(/rffp(1:9), Adense(iparents(2,i+9), iparents(4,i+9)), rffp(11:20)/)
   rmfp=(/rmfp(1:9), Adense(iparents(1,i+9), iparents(4,i+9)), rmfp(11:20)/)
   rfmp=(/rfmp(1:9), Adense(iparents(2,i+9), iparents(3,i+9)), rfmp(11:20)/)
!11
   rmmp=(/rmmp(1:10), Adense(iparents(1,i+10), iparents(3,i+10)), rmmp(12:20)/)
   rffp=(/rffp(1:10), Adense(iparents(2,i+10), iparents(4,i+10)), rffp(12:20)/)
   rmfp=(/rmfp(1:10), Adense(iparents(1,i+10), iparents(4,i+10)), rmfp(12:20)/)
   rfmp=(/rfmp(1:10), Adense(iparents(2,i+10), iparents(3,i+10)), rfmp(12:20)/)
!12
   rmmp=(/rmmp(1:11), Adense(iparents(1,i+11), iparents(3,i+11)), rmmp(13:20)/)
   rffp=(/rffp(1:11), Adense(iparents(2,i+11), iparents(4,i+11)), rffp(13:20)/)
   rmfp=(/rmfp(1:11), Adense(iparents(1,i+11), iparents(4,i+11)), rmfp(13:20)/)
   rfmp=(/rfmp(1:11), Adense(iparents(2,i+11), iparents(3,i+11)), rfmp(13:20)/)
!13
   rmmp=(/rmmp(1:12), Adense(iparents(1,i+12), iparents(3,i+12)), rmmp(14:20)/)
   rffp=(/rffp(1:12), Adense(iparents(2,i+12), iparents(4,i+12)), rffp(14:20)/)
   rmfp=(/rmfp(1:12), Adense(iparents(1,i+12), iparents(4,i+12)), rmfp(14:20)/)
   rfmp=(/rfmp(1:12), Adense(iparents(2,i+12), iparents(3,i+12)), rfmp(14:20)/)
!14
   rmmp=(/rmmp(1:13), Adense(iparents(1,i+13), iparents(3,i+13)), rmmp(15:20)/)
   rffp=(/rffp(1:13), Adense(iparents(2,i+13), iparents(4,i+13)), rffp(15:20)/)
   rmfp=(/rmfp(1:13), Adense(iparents(1,i+13), iparents(4,i+13)), rmfp(15:20)/)
   rfmp=(/rfmp(1:13), Adense(iparents(2,i+13), iparents(3,i+13)), rfmp(15:20)/)
!15
   rmmp=(/rmmp(1:14), Adense(iparents(1,i+14), iparents(3,i+14)), rmmp(16:20)/)
   rffp=(/rffp(1:14), Adense(iparents(2,i+14), iparents(4,i+14)), rffp(16:20)/)
   rmfp=(/rmfp(1:14), Adense(iparents(1,i+14), iparents(4,i+14)), rmfp(16:20)/)
   rfmp=(/rfmp(1:14), Adense(iparents(2,i+14), iparents(3,i+14)), rfmp(16:20)/)
!16
   rmmp=(/rmmp(1:15), Adense(iparents(1,i+15), iparents(3,i+15)), rmmp(17:20)/)
   rffp=(/rffp(1:15), Adense(iparents(2,i+15), iparents(4,i+15)), rffp(17:20)/)
   rmfp=(/rmfp(1:15), Adense(iparents(1,i+15), iparents(4,i+15)), rmfp(17:20)/)
   rfmp=(/rfmp(1:15), Adense(iparents(2,i+15), iparents(3,i+15)), rfmp(17:20)/)
!17
   rmmp=(/rmmp(1:16), Adense(iparents(1,i+16), iparents(3,i+16)), rmmp(18:20)/)
   rffp=(/rffp(1:16), Adense(iparents(2,i+16), iparents(4,i+16)), rffp(18:20)/)
   rmfp=(/rmfp(1:16), Adense(iparents(1,i+16), iparents(4,i+16)), rmfp(18:20)/)
   rfmp=(/rfmp(1:16), Adense(iparents(2,i+16), iparents(3,i+16)), rfmp(18:20)/)
!18
   rmmp=(/rmmp(1:17), Adense(iparents(1,i+17), iparents(3,i+17)), rmmp(19:20)/)
   rffp=(/rffp(1:17), Adense(iparents(2,i+17), iparents(4,i+17)), rffp(19:20)/)
   rmfp=(/rmfp(1:17), Adense(iparents(1,i+17), iparents(4,i+17)), rmfp(19:20)/)
   rfmp=(/rfmp(1:17), Adense(iparents(2,i+17), iparents(3,i+17)), rfmp(19:20)/)
!19
   rmmp=(/rmmp(1:18), Adense(iparents(1,i+18), iparents(3,i+18)), rmmp(20)/)
   rffp=(/rffp(1:18), Adense(iparents(2,i+18), iparents(4,i+18)), rffp(20)/)
   rmfp=(/rmfp(1:18), Adense(iparents(1,i+18), iparents(4,i+18)), rmfp(20)/)
   rfmp=(/rfmp(1:18), Adense(iparents(2,i+18), iparents(3,i+18)), rfmp(20)/)
!20
   rmmp=(/rmmp(1:19), Adense(iparents(1,i+19), iparents(3,i+19))/)
   rffp=(/rffp(1:19), Adense(iparents(2,i+19), iparents(4,i+19))/)
   rmfp=(/rmfp(1:19), Adense(iparents(1,i+19), iparents(4,i+19))/)
   rfmp=(/rfmp(1:19), Adense(iparents(2,i+19), iparents(3,i+19))/)

   d=((rmmp*rffp) + (rmfp*rfmp)) * 0.25
   answer=(/answer(1:i-1), d/)
 end do
end subroutine dij

