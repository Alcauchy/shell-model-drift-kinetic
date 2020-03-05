program main
  use HDF5
  use Shell
  USE ISO_C_BINDING
  use IO
  
  

  implicit none
  integer :: jk,ik,kk,n
  integer :: tempAr(2)
  real :: coords(3)
  integer, allocatable :: sizeAr(:)
  integer :: hi
  character :: ar(0:2) = (/'n','l','f'/)
  real, allocatable :: neighAr(:,:,:),neighAr1(:,:,:)
  type(Nodes) :: node
  n = 20
  node%numberOfShells = n
  call node%initializeNodes()
!  print *,shape(node%pairIndexList)
!  print *, 'hi there, here is you index list:', 8*n-1
!  do ik = 0,8*n-1
!    print *,node%pairIndexList(ik,:)
! end do
! print *,'neighbours for icosahedron:'
! do ik = 1,5
! print *,icosahedronList(:,ik,1,2)
!enddo
!
!print *,'neighbours for icosahedron:'
!do ik = 1,5
!print *,node%neighboursFirstShell(:,ik,1,2)
!enddo
!!
!!print *,'neighbours for dodecahedron:'
!!do ik = 1,3
!!print *,dodecahedronList(:,ik,2,1)
!!enddo
!print *, shape(node%conjFlagFirstShell)
!print *,'conjugated flag:'
!
!do n = 1,10
!        print *,n
!        print *,''
!        do ik = 1,3
!                print *,node%conjFlagSecondShell(:,ik,2,n)
!        enddo
!        print *,''
!enddo
!print *, icosahedronList(2,:,2,2)
!
!
!print *, node%pairIndexList(8,:)
!do ik = 0,2
!        print *,'look at ',ar(ik)
!        do jk = 0,8
!                tempAr = node%neighbours(8)%vector(ik,:,jk)
!                print *,tempAr
!        enddo
!        print *,''
!enddo
do ik = 0,n
!        !print *,''
        tempAr = node%pairIndexList(ik,:)
!        !print *, tempAr

        allocate(neighAr, SOURCE = node%getNeighboursCoordinates(ik))
        allocate(neighAr1, SOURCE = node%getInterNodesCoords(ik))
        allocate(sizeAr, SOURCE = shape(node%getNeighboursCoordinates(ik)))
        hi = sizeAr(3)
        !print *,shape(neighAr)
        coords = node%getCoordinates(ik)
        
        !print *,''
        do jk = 1,hi

                !do kk = 1,2
                        !print *, neighAr(:,kk,jk)
                !enddo
                !print *, (neighAr(:,1,jk)+neighAr(:,2,jk))+coords
                if (minval(abs(neighAr(:,1,jk)+neighAr(:,2,jk)+coords))>0.01) then
                        print *,ik,jk,'!'
                endif
                do kk = 0,2
                        if (neighAr(kk+1,1,jk)/=neighAr1(kk+1,1,jk)) then
                                print *,ik,jk,'!'
                        endif
                enddo
                !print *,''
        enddo
        deallocate(neighAr)
        deallocate(neighAr1)
        deallocate(sizeAr)
end do

call writeSpaceConfig(node)

end program main
