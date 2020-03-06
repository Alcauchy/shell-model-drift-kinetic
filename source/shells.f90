!
! flattened node - a number which obtained from shell number and vertice number: {n,l} -> flattened node
!
!
!
module Shell
        use Neighbours_lists
        private :: readLinksList,fillNeighbours,fillIndexNumber,fillCoordinates,fillFlags
        public :: getCoordinates,getNeighboursCoordinates

        type :: irregArray
                integer, allocatable :: vector(:,:,:)
        end type irregArray



        type Nodes
                integer :: numberOfShells
                real, allocatable :: basisIco(:,:),basisDod(:,:)                                                                    !                                                                                                                           ! list of flattened node numbers of interacting neighbours
                integer, allocatable :: pairIndexList(:,:)                                                                          ! list of {shell, vertice} numbers for each flattened node number
                integer, allocatable :: conjFlagFirstShell(:,:,:,:),conjFlagSecondShell(:,:,:,:)
                integer, allocatable :: neighboursFirstShell(:,:,:,:),neighboursSecondShell(:,:,:,:)
                real, allocatable :: coordinates(:,:)                                                                               ! list of nodes' coordinates
                type(irregArray), allocatable :: neighbours(:)
                type(irregArray), allocatable :: interactingNodesFlattened(:)


        contains
                procedure, pass(self) :: initializeNodes
                procedure, pass(self) :: readLinksList
                procedure, pass(self) :: fillNeighbours
                procedure, pass(self) :: fillIndexNumber
                procedure, pass(self) :: fillCoordinates
                procedure, pass(self) :: fillFlags
                procedure, pass(self) :: getCoordinates
                procedure, pass(self) :: getNeighboursCoordinates
                procedure, pass(self) :: getInterNodesFlat
                procedure, pass(self) :: getInterNodesCoords
        end type Nodes

        contains
                subroutine initializeNodes(self)
                        class(Nodes), intent(inout) :: self
                        allocate(self%basisIco(0:5,0:2))
                        allocate(self%basisDod(0:9,0:2))
                        call fillIndexNumber(self)
                        call initializePolyLists()
                        call getArbitaryPolyhedra(self%basisIco,self%basisDod)
                        call fillFlags(self)
                        call fillNeighbours(self)
                        call fillCoordinates(self)
                        call getInterNodesFlat(self)
                end subroutine initializeNodes

                subroutine readLinksList(self)
                        class(Nodes), intent(inout) :: self
                end subroutine readLinksList

                subroutine fillNeighbours(self)
                        class(Nodes), intent(inout) :: self
                        integer,parameter :: shapeFirstShell(4) = shape(icosahedronList)
                        integer,parameter :: shapeSecondShell(4) = shape(dodecahedronList)
                        integer :: numOfPairs(2)
                        integer :: k(2,shapeFirstShell(3))
                        numofPairs = (/shapeFirstShell(3)*shapeFirstShell(2),shapeSecondShell(3)*shapeSecondShell(2)/)
                        allocate(self%neighbours(0:self%numberOfShells*16-1))
                        k = reshape((/-2,-1,-1,1,1,2/),(/2,shapeFirstShell(3)/))
                        do i=0,self%numberOfShells*8-1
                                if (mod(self%pairIndexList(i,0),2) == 0) then

                                        allocate(self%neighbours(i)%vector(0:2,0:1,0:numOfPairs(1)-1))
                                        do j=0,shapeFirstShell(3)-1
                                                self%neighbours(i)%vector(0,0,(shapeFirstShell(2))*j : (shapeFirstShell(2))*(j + 1) - 1) = self%pairIndexList(i,0)+k(1,j+1)
                                                self%neighbours(i)%vector(0,1,(shapeFirstShell(2))*j : (shapeFirstShell(2))*(j + 1) - 1) = self%pairIndexList(i,0)+k(2,j+1)
                                        end do
                                        do j = 0,shapeFirstShell(3)-1
                                                do m = 0, shapeFirstShell(1)-1
                                                        self%neighbours(i)%vector(1,m,(shapeFirstShell(2))*j : (shapeFirstShell(2))*(j + 1) - 1) = self%neighboursFirstShell(m+1,:,j+1,self%pairIndexList(i,1)+1)
                                                        self%neighbours(i)%vector(2,m,(shapeFirstShell(2))*j : (shapeFirstShell(2))*(j + 1) - 1) = self%conjFlagFirstShell(m+1,:,j+1,self%pairIndexList(i,1)+1)
                                                end do
                                        end do
                                else if (mod(self%pairIndexList(i,0),2) == 1) then
                                        allocate(self%neighbours(i)%vector(0:2,0:1,0:numOfPairs(2)-1))
                                        do j=0,shapeSecondShell(3)-1
                                                self%neighbours(i)%vector(0,0,(shapeSecondShell(2))*j : (shapeSecondShell(2))*(j + 1) - 1) = self%pairIndexList(i,0)+k(1,j+1)
                                                self%neighbours(i)%vector(0,1,(shapeSecondShell(2))*j : (shapeSecondShell(2))*(j + 1) - 1) = self%pairIndexList(i,0)+k(2,j+1)
                                        end do
                                        do j = 0,shapeFirstShell(3)-1
                                                do m = 0, shapeFirstShell(1)-1
                                                        self%neighbours(i)%vector(1,m,(shapeSecondShell(2))*j : (shapeSecondShell(2))*(j + 1) - 1) = self%neighboursSecondShell(m+1,:,j+1,self%pairIndexList(i,1)+1)
                                                        self%neighbours(i)%vector(2,m,(shapeSecondShell(2))*j : (shapeSecondShell(2))*(j + 1) - 1) = self%conjFlagSecondShell(m+1,:,j+1,self%pairIndexList(i,1)+1)
                                                end do
                                        end do
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
                        allocate(self%coordinates(3,0:self%numberOfShells*16-1))
                        do i = 0,self%numberOfShells*8-1
                                self%coordinates(:,i) = getCoordinates(self,i)
                        end do
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
                
                subroutine getInterNodesFlat(self)
                        class(Nodes), intent(inout) :: self
                        integer :: extraElem
                        integer, allocatable :: sizeOfNeigh(:)
                        integer,allocatable :: n(:,:),l(:,:),f(:,:),newN(:,:),newL(:,:),newF(:,:)
                        logical,allocatable :: mask(:)
                        integer :: flatNum(2)
                        
                        allocate(self%interactingNodesFlattened(0:self%numberOfShells*16-1))
                        
                        do kNum = 0,self%numberOfShells*16-1
                                allocate(sizeOfNeigh,SOURCE = shape(self%neighbours(kNum)%vector))
                                allocate(n, SOURCE = self%neighbours(kNum)%vector(0,:,:))
                                allocate(l, SOURCE = self%neighbours(kNum)%vector(1,:,:))
                                allocate(f, SOURCE = self%neighbours(kNum)%vector(2,:,:))
                                
                                allocate(mask(sizeOfNeigh(3)))
                                
                                mask = merge(.true.,.false.,(n(1,:) < 0).or.(n(2,:) < 0))
                                extraElem = sizeOfNeigh(3)-count(mask)
                                if (extraElem /= 0) then
                                        allocate(newN(2,extraElem))
                                        allocate(newL(2,extraElem))
                                        allocate(newF(2,extraElem))
                                        newN(1,:) = pack(n(1,:),.not.mask)
                                        newN(2,:) = pack(n(2,:),.not.mask)
                                        newL(1,:) = pack(l(1,:),.not.mask)
                                        newL(2,:) = pack(l(2,:),.not.mask)
                                        newF(1,:) = pack(f(1,:),.not.mask)
                                        newF(2,:) = pack(f(2,:),.not.mask)
                                        deallocate(n)
                                        deallocate(l)
                                        deallocate(f)
                                        allocate(n, SOURCE = newN)
                                        allocate(l, SOURCE = newL)
                                        allocate(f, SOURCE = newF)
                                        sizeOfNeigh(3) = extraElem
                                        deallocate(newN)
                                        deallocate(newL)
                                        deallocate(newF)
                                        
                                
                                end if
                                deallocate(mask)
                                allocate(self%interactingNodesFlattened(kNum)%vector(0:1,0:1,sizeofNeigh(3)))
                                do j = 0,sizeOfNeigh(3)-1
                                        flatNum = floor(real(n(:,j+1)/2.)) * 16 + mod(n(:,j+1),2) * 6 + l(:,j+1)
                                        self%interactingNodesFlattened(kNum)%vector(0,0:1,j+1) = flatNum(:)
                                        self%interactingNodesFlattened(kNum)%vector(1,0:1,j+1) = f(:,j+1)
                                         
                                end do
                                deallocate(sizeOfNeigh)
                                deallocate(n)
                                deallocate(l)
                                deallocate(f)
                        end do
                end subroutine getInterNodesFlat





                function getCoordinates(self,kNum) result (k)
                        class(Nodes), intent(in) :: self
                        integer, intent(in) :: kNum
                        !real, intent(in) :: lambda(2)
                        real :: k(0:2)
                        integer :: n,l
                        n = self%pairIndexList(kNum,0)
                        l = self%pairIndexList(kNum,1)


                        if (mod(n,2) == 0) then
                                k = g**real(n)*lambda(0)*self%basisIco(l,:)
                        else if (mod(n,2) == 1) then
                                k = g**real(n)*lambda(1)*self%basisDod(l,:)
                        end if
                        !k = merge(g**n*lambda(0)*self%basisIco(:,l),g**n*lambda(1)*self%basisDod(:,l),mod(n,2) == 0)
                end function getCoordinates

                function getNeighboursCoordinates(self,kNum) result(k)
                        class(Nodes), intent(in) :: self
                        integer, intent(in) :: kNum
                        integer :: i,j,extraElem
                        integer, allocatable :: sizeOfNeigh(:)
                        integer,allocatable :: n(:,:),l(:,:),f(:,:),newN(:,:),newL(:,:),newF(:,:)
                        logical,allocatable :: mask(:)
                        logical,allocatable :: lmask(:)
                        real, allocatable :: k(:,:,:)

                        allocate(sizeOfNeigh,SOURCE = shape(self%neighbours(kNum)%vector))
                        allocate(n, SOURCE = self%neighbours(kNum)%vector(0,:,:))
                        allocate(l, SOURCE = self%neighbours(kNum)%vector(1,:,:))
                        allocate(f, SOURCE = self%neighbours(kNum)%vector(2,:,:))

                        allocate(mask(sizeOfNeigh(3)))

                        mask = merge(.true.,.false.,(n(1,:) < 0).or.(n(2,:) < 0))
                        extraElem = sizeOfNeigh(3)-count(mask)

                        if (extraElem /= 0) then
                                allocate(newN(2,extraElem))
                                allocate(newL(2,extraElem))
                                allocate(newF(2,extraElem))
                                newN(1,:) = pack(n(1,:),.not.mask)
                                newN(2,:) = pack(n(2,:),.not.mask)
                                newL(1,:) = pack(l(1,:),.not.mask)
                                newL(2,:) = pack(l(2,:),.not.mask)
                                newF(1,:) = pack(f(1,:),.not.mask)
                                newF(2,:) = pack(f(2,:),.not.mask)
                                deallocate(n)
                                deallocate(l)
                                deallocate(f)
                                allocate(n, SOURCE = newN)
                                allocate(l, SOURCE = newL)
                                allocate(f, SOURCE = newF)
                                sizeOfNeigh(3) = extraElem
                                !do i = 1,extraElem
                                !        print *,n(:,i)
                                !enddo
                        end if
                        allocate(k(0:2,0:sizeOfNeigh(2)-1,0:sizeOfNeigh(3)-1))

                        !newN(2,:) = pack(n(2,:),mask)
                        !if (max) ==
                        !print*,'n ',(n(:,4:6))
                        do j = 0,sizeOfNeigh(3)-1
                                !print*,(n(:,j+1))
                                do i = 0,sizeOfNeigh(2)-1

                                        if (mod(n(i+1,j+1),2) == 0) then
                                                k(:,i,j) = g**(n(i+1,j+1))*lambda(0)*self%basisIco(l(i+1,j+1),:)
                                        else if (mod(n(i+1,j+1),2) == 1) then
                                                k(:,i,j) = g**(n(i+1,j+1))*lambda(1)*self%basisDod(l(i+1,j+1),:)
                                        end if
                                        k(:,i,j) = merge(-k(:,i,j),k(:,i,j),f(i+1,j+1) == 1)
                                        !k(:,:,j) = merge(-k(:,:,j),k(:,:,j),f(:,j+1) == 1)
                                end do
                        end do
                end function getNeighboursCoordinates
                        
                
                function getInterNodesCoords(self, kNum) result(k)
                        class(Nodes), intent(in) :: self
                        integer, intent(in) :: kNum
                        real, allocatable :: k(:,:,:)
                        integer, allocatable :: sizeOfNeigh(:)
                        integer, allocatable :: interactingNodes(:,:,:)
                        
                        allocate(sizeOfNeigh,SOURCE = shape(self%interactingNodesFlattened(kNum)%vector))
                        allocate(interactingNodes, SOURCE = self%interactingNodesFlattened(kNum)%vector)
                        allocate(k(0:2,0:sizeOfNeigh(2)-1,0:sizeOfNeigh(3)-1))
                        do i = 0,sizeOfNeigh(3)-1
                                do j = 0,sizeOfNeigh(2)-1
                                        k(:,j,i) = self%coordinates(:,interactingNodes(0,j,i+1))
                                        k(:,j,i) = merge(-k(:,j,i),k(:,j,i),interactingNodes(1,j,i+1) == 1)
                                        !k(:,j,i) = 1.0
                                end do
                                        
                        end do 
                        
                        !do i = 0,sizeOfNeigh(3)-1
                        !        do j = 0,sizeOfNeigh(2)-1
                        !        k(:,j,i) = getCoordinates(self,interactingNodes(0,j,i+1))
                        !        k(:,j,i) = merge(-k(:,j,i),k(:,j,i),interactingNodes(1,j,i+1) == 1)
                        !        end do
                        
                        !end do 
                end function getInterNodesCoords
                


end module Shell


