module IO
        use HDF5
        use h5d
        USE ISO_C_BINDING
        !use h5lib
        !use h5f
        use Shell
        
        INTEGER, PARAMETER :: int_kind_16 = SELECTED_INT_KIND(18)
        INTEGER, PARAMETER :: real_kind_15 = SELECTED_REAL_KIND(15,307)
        CHARACTER(LEN=21), PARAMETER :: filename  = "spaceConfiguration.h5"
        
        public :: writeSpaceConfig
        public :: writeCoordinates
        

         

contains
        subroutine writeSpaceConfig(typeNodes)
                class(Nodes), intent(in) :: typeNodes
                
                integer :: length
                integer :: error
                CHARACTER(LEN=3) , PARAMETER     :: dataset   = "links"
                
                integer(hsize_t), DIMENSION(3)   :: arraydims
                integer                          :: rank = 3
                TYPE(hvl_t), allocatable, TARGET :: wdata(:)
                INTEGER(hsize_t), DIMENSION(1:1) :: dims = (/1/)
                INTEGER(HID_T)                   :: file, filetype,plsit_id, datatype, memtype, space, dset ! Handles
                INTEGER(HID_T), allocatable      :: dt_id(:),st_id(:)
                
                integer(8)                       ::  offset
                integer(8), allocatable          :: sz(:)
                integer(hsize_t), dimension(1:1) :: ms_dims
                integer, allocatable, target     :: denarray(:,:,:)
                CHARACTER(LEN=15)                :: strName
                
                ! allocate datatype id's, space id's and size of arrays
                length = typeNodes%NumberOfShells*8
                allocate( dt_id(0:length-1) )
                allocate( st_id(0:length-1) )
                allocate( sz(0:length-1))
                ms_dims = (/1/)
                
                ! open HDF API
                call h5open_f(error)
                
                ! create file, truncate if exists
                call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, error)
                
                !create simple dataspace
                call h5screate_simple_f(1, ms_dims, space, error)
                
                ! make a layout for a dataspace 
                offset = 0
                do i = 0,length-1    
                        arraydims = shape(typeNodes%interactingNodesFlattened(i)%vector)
                        write(strName,*)i
                        call H5Tarray_create_f(H5T_NATIVE_INTEGER, rank, arraydims, dt_id(i), error);
                        call H5Tget_size_f(dt_id(i), sz(i), error)
                        offset = offset + sz(i)
                        
                        
                enddo
                
                ! create datatype of "memtype" for arrays storing 
                call H5Tcreate_f(H5T_COMPOUND_F, offset, memtype, error);
                offset = 0
                 do i = 0,length-1
                         arraydims = shape(typeNodes%interactingNodesFlattened(i)%vector)
                         write(strName,*)i
                         call H5Tinsert_f(memtype, strName, offset,  dt_id(i), error)
                         offset = offset+sz(i)
                enddo
                ! create dataset
                CALL h5dcreate_f(file, dataset, memtype, space, dset, error)
                
                ! create memory types and write data to dataset
                offset = 0
                do i = 0,length-1
                        allocate(denarray, SOURCE = typeNodes%interactingNodesFlattened(i)%vector)
                        write(strName,*)i
                        call h5tcreate_f(H5T_COMPOUND_F, sz(i), st_id(i), error)
                        call h5tinsert_f(st_id(i), strName, offset, dt_id(i), error)
                        call h5dwrite_f(dset, st_id(i), denarray, arraydims, error)
                        deallocate(denarray)
                enddo

                CALL h5tclose_f(memtype, error)
                CALL h5sclose_f(space, error)
                CALL h5dclose_f(dset , error)
                CALL h5fclose_f(file , error)
        end subroutine writeSpaceConfig
        
        subroutine writeCoordinates(typeNodes)
                class(Nodes), intent(in) :: typeNodes
                integer                  :: length
                
        end subroutine writeCoordinates
        

end module IO 
