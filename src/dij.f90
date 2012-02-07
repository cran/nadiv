subroutine dij(iparents, iparents_r, dim_aij, rcA, aij, N, answer)

 implicit none
!define inputs
 integer, intent(in) :: iparents(4, iparents_r), iparents_r, dim_aij, rcA(2,dim_aij), N
 real, intent(in) :: aij(dim_aij)
 real, intent(out) :: answer(iparents_r) 
!define variables derived in subroutine
 integer :: i, itmp, j, k, mm(2), ff(2), mf(2), fm(2)
 real, dimension(1:20) :: rmmp, rffp, rmfp, rfmp, d

 type row
    real, pointer :: R(:)
 end type
 type(row) :: A(N)

 do j=1,N
    allocate (A(j)%R(1:j))
    A(j)%R(1:j)=0
 end do

 do k=1,dim_aij
    A(rcA(1,k))%R(rcA(2,k))=aij(k)
 end do

  
!now the program that computes Dij for each j individual related to individual i
 do i=1,iparents_r,20
     rmmp=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
     rffp=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
     rmfp=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
     rfmp=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
!1
    mm=(/max(iparents(1,i),iparents(3,i)), min(iparents(1,i),iparents(3,i))/)
    ff=(/max(iparents(2,i),iparents(4,i)), min(iparents(2,i),iparents(4,i))/)
    mf=(/max(iparents(1,i),iparents(4,i)), min(iparents(1,i),iparents(4,i))/)
    fm=(/max(iparents(2,i),iparents(3,i)), min(iparents(2,i),iparents(3,i))/)
   rmmp=(/A(mm(1))%R(mm(2)), rmmp(2:20)/)
   rffp=(/A(ff(1))%R(ff(2)), rffp(2:20)/)
   rmfp=(/A(mf(1))%R(mf(2)), rmfp(2:20)/)
   rfmp=(/A(fm(1))%R(fm(2)), rfmp(2:20)/)
!2
    itmp=i+1.0
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1), A(mm(1))%R(mm(2)), rmmp(3:20)/)
   rffp=(/rffp(1), A(ff(1))%R(ff(2)), rffp(3:20)/)
   rmfp=(/rmfp(1), A(mf(1))%R(mf(2)), rmfp(3:20)/)
   rfmp=(/rfmp(1), A(fm(1))%R(fm(2)), rfmp(3:20)/)
!3
    itmp=i+2
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:2), A(mm(1))%R(mm(2)), rmmp(4:20)/)
   rffp=(/rffp(1:2), A(ff(1))%R(ff(2)), rffp(4:20)/)
   rmfp=(/rmfp(1:2), A(mf(1))%R(mf(2)), rmfp(4:20)/)
   rfmp=(/rfmp(1:2), A(fm(1))%R(fm(2)), rfmp(4:20)/)
!4
    itmp=i+3
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:3), A(mm(1))%R(mm(2)), rmmp(5:20)/)
   rffp=(/rffp(1:3), A(ff(1))%R(ff(2)), rffp(5:20)/)
   rmfp=(/rmfp(1:3), A(mf(1))%R(mf(2)), rmfp(5:20)/)
   rfmp=(/rfmp(1:3), A(fm(1))%R(fm(2)), rfmp(5:20)/)
!5
    itmp=i+4
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:4), A(mm(1))%R(mm(2)), rmmp(6:20)/)
   rffp=(/rffp(1:4), A(ff(1))%R(ff(2)), rffp(6:20)/)
   rmfp=(/rmfp(1:4), A(mf(1))%R(mf(2)), rmfp(6:20)/)
   rfmp=(/rfmp(1:4), A(fm(1))%R(fm(2)), rfmp(6:20)/)
!6
    itmp=i+5
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:5), A(mm(1))%R(mm(2)), rmmp(7:20)/)
   rffp=(/rffp(1:5), A(ff(1))%R(ff(2)), rffp(7:20)/)
   rmfp=(/rmfp(1:5), A(mf(1))%R(mf(2)), rmfp(7:20)/)
   rfmp=(/rfmp(1:5), A(fm(1))%R(fm(2)), rfmp(7:20)/)
!7
    itmp=i+6
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:6), A(mm(1))%R(mm(2)), rmmp(8:20)/)
   rffp=(/rffp(1:6), A(ff(1))%R(ff(2)), rffp(8:20)/)
   rmfp=(/rmfp(1:6), A(mf(1))%R(mf(2)), rmfp(8:20)/)
   rfmp=(/rfmp(1:6), A(fm(1))%R(fm(2)), rfmp(8:20)/)
!8
    itmp=i+7
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:7), A(mm(1))%R(mm(2)), rmmp(9:20)/)
   rffp=(/rffp(1:7), A(ff(1))%R(ff(2)), rffp(9:20)/)
   rmfp=(/rmfp(1:7), A(mf(1))%R(mf(2)), rmfp(9:20)/)
   rfmp=(/rfmp(1:7), A(fm(1))%R(fm(2)), rfmp(9:20)/)
!9
    itmp=i+8
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:8), A(mm(1))%R(mm(2)), rmmp(10:20)/)
   rffp=(/rffp(1:8), A(ff(1))%R(ff(2)), rffp(10:20)/)
   rmfp=(/rmfp(1:8), A(mf(1))%R(mf(2)), rmfp(10:20)/)
   rfmp=(/rfmp(1:8), A(fm(1))%R(fm(2)), rfmp(10:20)/)
!10
    itmp=i+9
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:9), A(mm(1))%R(mm(2)), rmmp(11:20)/)
   rffp=(/rffp(1:9), A(ff(1))%R(ff(2)), rffp(11:20)/)
   rmfp=(/rmfp(1:9), A(mf(1))%R(mf(2)), rmfp(11:20)/)
   rfmp=(/rfmp(1:9), A(fm(1))%R(fm(2)), rfmp(11:20)/)
!11
    itmp=i+10
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:10), A(mm(1))%R(mm(2)), rmmp(12:20)/)
   rffp=(/rffp(1:10), A(ff(1))%R(ff(2)), rffp(12:20)/)
   rmfp=(/rmfp(1:10), A(mf(1))%R(mf(2)), rmfp(12:20)/)
   rfmp=(/rfmp(1:10), A(fm(1))%R(fm(2)), rfmp(12:20)/)
!12
    itmp=i+11
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:11), A(mm(1))%R(mm(2)), rmmp(13:20)/)
   rffp=(/rffp(1:11), A(ff(1))%R(ff(2)), rffp(13:20)/)
   rmfp=(/rmfp(1:11), A(mf(1))%R(mf(2)), rmfp(13:20)/)
   rfmp=(/rfmp(1:11), A(fm(1))%R(fm(2)), rfmp(13:20)/)
!13
    itmp=i+12
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:12), A(mm(1))%R(mm(2)), rmmp(14:20)/)
   rffp=(/rffp(1:12), A(ff(1))%R(ff(2)), rffp(14:20)/)
   rmfp=(/rmfp(1:12), A(mf(1))%R(mf(2)), rmfp(14:20)/)
   rfmp=(/rfmp(1:12), A(fm(1))%R(fm(2)), rfmp(14:20)/)
!14
    itmp=i+13
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:13), A(mm(1))%R(mm(2)), rmmp(15:20)/)
   rffp=(/rffp(1:13), A(ff(1))%R(ff(2)), rffp(15:20)/)
   rmfp=(/rmfp(1:13), A(mf(1))%R(mf(2)), rmfp(15:20)/)
   rfmp=(/rfmp(1:13), A(fm(1))%R(fm(2)), rfmp(15:20)/)
!15
    itmp=i+14
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:14), A(mm(1))%R(mm(2)), rmmp(16:20)/)
   rffp=(/rffp(1:14), A(ff(1))%R(ff(2)), rffp(16:20)/)
   rmfp=(/rmfp(1:14), A(mf(1))%R(mf(2)), rmfp(16:20)/)
   rfmp=(/rfmp(1:14), A(fm(1))%R(fm(2)), rfmp(16:20)/)
!16
    itmp=i+15
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:15), A(mm(1))%R(mm(2)), rmmp(17:20)/)
   rffp=(/rffp(1:15), A(ff(1))%R(ff(2)), rffp(17:20)/)
   rmfp=(/rmfp(1:15), A(mf(1))%R(mf(2)), rmfp(17:20)/)
   rfmp=(/rfmp(1:15), A(fm(1))%R(fm(2)), rfmp(17:20)/)
!17
    itmp=i+16
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:16), A(mm(1))%R(mm(2)), rmmp(18:20)/)
   rffp=(/rffp(1:16), A(ff(1))%R(ff(2)), rffp(18:20)/)
   rmfp=(/rmfp(1:16), A(mf(1))%R(mf(2)), rmfp(18:20)/)
   rfmp=(/rfmp(1:16), A(fm(1))%R(fm(2)), rfmp(18:20)/)
!18
    itmp=i+17
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:17), A(mm(1))%R(mm(2)), rmmp(19:20)/)
   rffp=(/rffp(1:17), A(ff(1))%R(ff(2)), rffp(19:20)/)
   rmfp=(/rmfp(1:17), A(mf(1))%R(mf(2)), rmfp(19:20)/)
   rfmp=(/rfmp(1:17), A(fm(1))%R(fm(2)), rfmp(19:20)/)
!19
    itmp=i+18
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:18), A(mm(1))%R(mm(2)), rmmp(20)/)
   rffp=(/rffp(1:18), A(ff(1))%R(ff(2)), rffp(20)/)
   rmfp=(/rmfp(1:18), A(mf(1))%R(mf(2)), rmfp(20)/)
   rfmp=(/rfmp(1:18), A(fm(1))%R(fm(2)), rfmp(20)/)
!20
    itmp=i+19
    mm=(/max(iparents(1,itmp),iparents(3,itmp)), min(iparents(1,itmp),iparents(3,itmp))/)
    ff=(/max(iparents(2,itmp),iparents(4,itmp)), min(iparents(2,itmp),iparents(4,itmp))/)
    mf=(/max(iparents(1,itmp),iparents(4,itmp)), min(iparents(1,itmp),iparents(4,itmp))/)
    fm=(/max(iparents(2,itmp),iparents(3,itmp)), min(iparents(2,itmp),iparents(3,itmp))/)
   rmmp=(/rmmp(1:19), A(mm(1))%R(mm(2))/)
   rffp=(/rffp(1:19), A(ff(1))%R(ff(2))/)
   rmfp=(/rmfp(1:19), A(mf(1))%R(mf(2))/)
   rfmp=(/rfmp(1:19), A(fm(1))%R(fm(2))/)


   d=(rmmp*rffp) + (rmfp*rfmp)
   answer=(/answer(1:i-1), d/)
 end do
end subroutine dij

