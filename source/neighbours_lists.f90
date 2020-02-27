module Neighbours_lists
        use Constants
        real,parameter :: g = sqrt((1.+sqrt(5.))/2.)
        real,parameter :: phi = (1.+sqrt(5.))/2.
        real,parameter :: alpha = asin(phi/sqrt(3.)) - acos(phi/sqrt(phi+2.))
        real,parameter :: beta = atan(2*phi*phi)
        real,parameter :: gamma = PI/2.-atan(0.5)
        real,parameter :: lambda(0:1) =(/sqrt(sqrt(5.)/3.),1./) 

        public :: initializePolyLists, getArbitaryPolyhedra
        private :: getArbitaryIcosahedron,getArbitaryDodecahedron
        integer,protected :: dodecahedronList(2,3,3,10),icosahedronList(2,5,3,6)

contains
        subroutine initializePolyLists()

                dodecahedronList = reshape((/15,  6, 11,  7, 14,  8,  4,  6,  9,  7, 11,  8, 10, 15,  3, 11,  5,&
                                                14, 16,  6, 12,  8, 10,  9,  5,  6, 10,  8,  7,  9, 11, 16,  4, 12,&
                                                1, 10, 17,  6, 13,  9, 11, 10,  1,  6, 11,  9,  8, 10,  7, 17,  5,&
                                                13,  2, 11, 18,  6, 14, 10, 12, 11,  2,  6,  7, 10,  9, 11,  8, 18,&
                                                1, 14,  3, 12, 19,  6, 13,  7, 10, 11,  3,  6, 10,  7,  8, 11,  9,&
                                                19,  4, 13,  2, 10,  8,  7,  7,  8, 10,  4,  5,  7,  3,  8,  6,  4,&
                                                11,  8,  9,  7,  0, 10,  9,  8,  8,  9, 11,  5,  1,  8,  4,  9,  6,&
                                                5,  7,  9, 10,  8,  0, 11, 12,  1,  5,  9,  9, 10,  6,  1,  2,  9,&
                                                5, 10,  0, 12,  8,  5, 11,  9, 13,  2,  6, 10,  5, 11,  6,  2,  3,&
                                                10,  1, 11,  0, 13,  9,  6,  7,  5,  6,  7, 14,  3,  7, 11,  2,  7,&
                                                6,  3,  4, 11,  8,  6,  0, 14, 10,  7/),(/2, 3, 3, 10/))

                icosahedronList = reshape((/10, 10, 11, 11,  7, 12,  8, 13,  9, 14,  5, 10,  6, 11,  7, 12,  8,&
                                                13,  9, 14, 15, 10, 16, 11, 17,  7, 18,  8, 19,  9,  3, 10,  4, 14,&
                                                11, 15,  6,  7,  8, 19,  1, 10,  3, 14, 18, 15, 12,  7, 16, 19, 11,&
                                                3, 13,  4,  8, 11,  2,  6,  6,  8,  5, 10,  4, 11,  9, 15,  7, 16,&
                                                6,  8,  4, 10,  2, 11, 17, 15, 19, 16, 13,  8, 14,  5, 12,  4,  7,&
                                                9,  9,  7,  3,  6,  1, 11,  5, 12, 10, 16,  8, 17,  6,  9,  0, 11,&
                                                3, 12, 18, 16, 15, 17, 14,  9, 10,  1, 13,  5,  8, 10,  5,  8,  4,&
                                                6,  6,  5,  1, 13,  2, 12,  9, 18, 11, 17, 10,  5,  4, 13,  1, 12,&
                                                16, 18, 19, 17,  0,  6, 14,  1, 11,  2,  6,  9,  9, 11,  6,  6,  7,&
                                                18,  2, 14,  3, 13, 10, 19, 11,  6, 15, 18,  0, 14,  2, 13, 17, 19,&
                                                1,  6,  5,  7, 10,  2, 12,  3,  7, 10/),(/2, 5, 3, 6/))

        end subroutine initializePolyLists

        subroutine getArbitaryPolyhedra(k1,k2)
                real, intent(inout) :: k1(0:2,0:5),k2(0:2,0:9)
                call getArbitaryIcosahedron(k1)
                call getArbitaryDodecahedron(k2)
        end subroutine getArbitaryPolyhedra

        subroutine getArbitaryIcosahedron(k)
                real, intent(inout) :: k(0:5,0:2)
                real :: theta(0:5)
                real :: psi(0:5)
                theta(0) = 0.
                theta(1:5) = gamma
                psi(0:1) = 0
                psi(2:5) = (/2.,4.,6.,8./)*PI/5.
                k = reshape((/sin(theta)*cos(psi),sin(theta)*sin(psi),cos(theta)/),(/6,3/))
        end subroutine getArbitaryIcosahedron


        subroutine getArbitaryDodecahedron(k)
                real, intent(inout) :: k(0:9,0:2)
                real :: theta(0:9)
                real :: psi(0:9)
                theta(:4) = alpha
                theta(5:9) = beta
                psi(:4) = (/1.,3.,5.,7.,9./)*PI/5.
                psi(5:9) = psi(:4)
                k = reshape((/sin(theta)*cos(psi),sin(theta)*sin(psi),cos(theta)/),(/10,3/))
                !print ()
        end subroutine getArbitaryDodecahedron

end module Neighbours_lists
