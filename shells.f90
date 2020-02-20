!
! flattened node - a number which obtained from shell number and vertice number: {n,l} -> flattened node
!
!
!
module Shell
        use Neighbours_lists
        private:: readLinksList,fillNeighbours,fillIndexNumber,fillCoordinates,fillFlags

        type :: irregArray
                integer, allocatable :: vector(:,:)
        end type irregArray



        type Nodes
                integer :: numberOfShells
                real, allocatable :: basisIco(:,:),basisDod(:,:)                                                                    !                                                                                                                           ! list of flattened node numbers of interacting neighbours
                integer, allocatable :: pairIndexList(:,:)                                                                          ! list of {shell, vertice} numbers for each flattened node number
                integer, allocatable :: conjFlagFirstShell(:,:,:,:),conjFlagSecondShell(:,:,:,:)
                integer, allocatable :: neighboursFirstShell(:,:,:,:),neighboursSecondShell(:,:,:,:)
                real, allocatable :: coordinates(:,:)                                                                               ! list of nodes' coordinates
                type(irregArray), allocatable :: neighbours(:)


        contains
                procedure, pass(self) :: initializeNodes
                procedure, pass(self) :: readLinksList
                procedure, pass(self) :: fillNeighbours
                procedure, pass(self) :: fillIndexNumber
                procedure, pass(self) :: fillCoordinates
                procedure, pass(self) :: fillFlags
        end type Nodes

        contains
                subroutine initializeNodes(self)
                        class(Nodes), intent(inout) :: self
                        allocate(self%basisIco(0:2,0:5))
                        allocate(self%basisDod(0:2,0:9))
                        call fillIndexNumber(self)
                        call initializePolyLists()
                        call getArbitaryPolyhedra(self%basisIco,self%basisDod)
                        call fillFlags(self)
                end subroutine initializeNodes

                subroutine readLinksList(self)
                        class(Nodes), intent(inout) :: self
                end subroutine readLinksList

                subroutine fillNeighbours(self)
                        class(Nodes), intent(inout) :: self
                        allocate(self%neighbours(0:self%numberOfShells*16-1))

                        do i=0,self%numberOfShells*16-1
                                if (mod(i,2) == 1) then
                                        allocate(self%neighbours(i)%vector(0:2,0:5))
                                else if (mod(i,2) == 0) then
                                        allocate(self%neighbours(i)%vector(0:2,0:10))
                                end if

                        end do

                end subroutine fillNeighbours

                subroutine fillIndexNumber(self)
                        class(Nodes), intent(inout) :: self
                        integer,parameter :: icoVertices = 6
                        integer,parameter :: dodVertices = 10                                                                       ! number of vertices for each polyhedra devided by 2
                        integer :: verticeIndexArray(0:icoVertices+dodVertices-1)                                                   ! list of polyhedra vertices indecies
                        integer :: dodTile(0:dodVertices-1,0:1)
                        integer :: icoTile(0:icoVertices-1,0:1)

                        verticeIndexArray(:) = (/(i, i = 0,(icoVertices+dodVertices-1)) /)
                        icoTile(:,1) = verticeIndexArray(0:icoVertices-1)
                        dodTile(:,1) = verticeIndexArray(0:dodVertices-1)

                        allocate(self%pairIndexList(0:self%numberOfShells*(icoVertices+dodVertices)/2-1,0:1))

                        do i = 0,self%numberOfShells/2-1
                                icoTile(:,0) = 2*i
                                dodTile(:,0) = 2*i+1

                                self%pairIndexList(icoVertices * i + dodVertices * i : icoVertices * (i + 1) + dodVertices * i - 1, : ) &
                                                                        = icoTile(:,:)
                                self%pairIndexList(icoVertices * (i + 1) + dodVertices * i : icoVertices * (i+1) + dodVertices * (i + 1) - 1, : ) &
                                                                        = dodTile(:,:)
                        end do
                end subroutine fillIndexNumber

                subroutine fillCoordinates(self)
                        class(Nodes), intent(inout) :: self
                end subroutine fillCoordinates

                subroutine fillFlags(self)
                        class(Nodes), intent(inout) :: self
                        integer :: shapeFirstShell(4) = shape(icosahedronList)
                        integer :: shapeSecondShell(4) = shape(dodecahedronList)
                        integer :: cutoff(2)
                        allocate(self%conjFlagFirstShell, MOLD = icosahedronList)
                        allocate(self%conjFlagSecondShell, MOLD = dodecahedronList)
                        allocate(self%neighboursFirstShell, MOLD = icosahedronList)
                        allocate(self%neighboursSecondShell, MOLD = dodecahedronList)
                        do i = 1,3
                                if (i == 1) then
                                        cutoff = (/shapeFirstShell(4),shapeSecondShell(4)/)
                                else if (i == 2) then
                                        cutoff = (/shapeSecondShell(4),shapeSecondShell(4)/)
                                else if (i == 3) then
                                        cutoff = (/shapeSecondShell(4),shapeFirstShell(4)/)
                                end if

                                do j = 1,shapeFirstShell(1)
                                        self%conjFlagFirstShell(j,:,i,:) = merge(1,0,icosahedronList(j,:,i,:) >= cutoff(j))
                                        self%neighboursFirstShell(j,:,i,:) = merge(icosahedronList(j,:,i,:) - cutoff(j),icosahedronList(j,:,i,:),icosahedronList(j,:,i,:) >= cutoff(j))

                                end do
                        end do

                        do i = 1,3
                                if (i == 1) then
                                        cutoff = (/shapeSecondShell(4),shapeFirstShell(4)/)
                                else if (i == 2) then
                                        cutoff = (/shapeFirstShell(4),shapeFirstShell(4)/)
                                else if (i == 3) then
                                        cutoff = (/shapeFirstShell(4),shapeSecondShell(4)/)
                                end if
                                do j = 1,shapeFirstShell(1)
                                        self%conjFlagSecondShell(j,:,i,:) = merge(1,0,dodecahedronList(j,:,i,:) >= cutoff(j))
                                        self%neighboursSecondShell(j,:,i,:) = merge(dodecahedronList(j,:,i,:) - cutoff(j),dodecahedronList(j,:,i,:),dodecahedronList(j,:,i,:) >= cutoff(j))
                                end do
                        end do


                end subroutine fillFlags





end module Shell

program main
  use Shell
  implicit none
  integer :: n,ik
  type(Nodes) :: node
  n = 10
  node%numberOfShells = n
  call node%initializeNodes()
!  print *,shape(node%pairIndexList)
!  print *, 'hi there, here is you index list:', 8*n-1
!  do ik = 0,8*n-1
!    print *,node%pairIndexList(ik,:)
! end do
 print *,'neighbours for icosahedron:'
 do ik = 1,5
 print *,icosahedronList(:,ik,1,2)
enddo

print *,'neighbours for icosahedron:'
do ik = 1,5
print *,node%neighboursFirstShell(:,ik,1,2)
enddo
!
!print *,'neighbours for dodecahedron:'
!do ik = 1,3
!print *,dodecahedronList(:,ik,2,1)
!enddo
print *, shape(node%conjFlagFirstShell)
print *,'conjugated flag:'

do n = 1,10
        print *,n
        print *,''
        do ik = 1,3
                print *,node%conjFlagSecondShell(:,ik,2,n)
        enddo
        print *,''
enddo
print *, icosahedronList(2,:,2,2)
end program main
