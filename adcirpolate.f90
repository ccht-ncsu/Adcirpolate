!
! The following include file contains all of the C preprocessor macros used in the code.
!
#include "c_preprocessor.h"

!> @author Ali Samii - 2016 - Updated 2018
!! Ali Samii - Department of ASE/EM, UT Austin
!! @brief This module is an interface between parallel Adcirc input files and ESMF library.
module adcirpolate

   use ESMF
   use MPI
   use errors_and_msgs

   !> \author Ali Samii - 2016
   !! \brief This object stores the data required for construction of a parallel or serial
   !! ESMF_Mesh from <tt>fort.14, fort.18, partmesh.txt</tt> files.
   !!
   type meshdata
      !> \details vm is an ESMF_VM object.  ESMF_VM is just an ESMF virtual machine class,
      !! which we will use to get the data about the local PE and PE count.
      type(ESMF_VM)                      :: vm
      !> \details This array contains the node coordinates of the mesh. For
      !! example, in a 2D mesh, the \c jth coordinate of the \c nth node
      !! is stored in location <tt> 2*(n-1)+j</tt> of this array.
      real(ESMF_KIND_R8), allocatable    :: NdCoords(:)
      !> \details This array contains the elevation of different nodes of the mesh
      real(ESMF_KIND_R8), allocatable    :: bathymetry(:)
      !> \details Number of nodes present in the current PE. This is different from the
      !! number of nodes owned by this PE (cf. NumOwnedNd)
      integer(ESMF_KIND_I4)              :: NumNd
      !> \details Number of nodes owned by this PE. This is different from the number of
      !! nodes present in the current PE (cf. NumNd)
      integer(ESMF_KIND_I4)              :: NumOwnedNd
      !> \details Number of elements in the current PE. This includes ghost elements and
      !! owned elements. However, we do not bother to distinguish between owned
      !! element and present element (as we did for the nodes).
      integer(ESMF_KIND_I4)              :: NumEl
      !> \details Number of nodes of each element, which is simply three in 2D ADCIRC.
      integer(ESMF_KIND_I4)              :: NumND_per_El
      !> \details Global node numbers of the nodes which are present in the current PE.
      integer(ESMF_KIND_I4), allocatable :: NdIDs(:)
      !> \details Global element numbers which are present in the current PE.
      integer(ESMF_KIND_I4), allocatable :: ElIDs(:)
      !> \details The element connectivity array, for the present elements in the current PE.
      !! The node numbers are the local numbers of the present nodes. All the element
      !! connectivities are arranged in this one-dimensional array.
      integer(ESMF_KIND_I4), allocatable :: ElConnect(:)
      !> \details The number of the PE's which own each of the nodes present this PE.
      !! This number is zero-based.
      integer(ESMF_KIND_I4), allocatable :: NdOwners(:)
      !> \details An array containing the element types, which are all triangles in our
      !! application.
      integer(ESMF_KIND_I4), allocatable :: ElTypes(:)
      !> \details This is an array, which maps the indices of the owned nodes to the indices of the present
      !! nodes. For example, assume we are on <tt>PE = 1</tt>, and we have four nodes present, and the
      !! first and third nodes belong to <tt>PE = 0</tt>. So we have:
      !! \code
      !! NumNd = 4
      !! NumOwnedNd = 2
      !! NdOwners = (/0, 1, 0, 1/)
      !! NdIDs = (/2, 3, 5, 6/)
      !! owned_to_present = (/2, 4/)    <-- Because the first node owned by this PE is actually
      !!                                    the second node present on this PE, and so on.
      !! \endcode
      integer(ESMF_KIND_I4), allocatable :: owned_to_present_nodes(:)
      !> \details The directory where the files are located
      character(len=:), allocatable      :: dir_name
      !> \details This flag tells if the meshdata has been initialized
      logical                            :: is_initialized = .false.
   end type meshdata

   !>
   !! \author Ali Samii - 2016
   !! This structure stores the data from an ADCIRC hotstart file. To know more about
   !! different members of this stucture, consult ADCIRC manual or user refernce.
   type hotdata
      real(ESMF_KIND_R8)                 :: TimeLoc
      real(ESMF_KIND_R8), allocatable    :: ETA1(:), ETA2(:), ETADisc(:), UU2(:), VV2(:), CH1(:), realNODECODE(:)
      integer(ESMF_KIND_I4), allocatable :: NNODECODE(:), NOFF(:)
      integer(ESMF_KIND_I4)              :: InputFileFmtVn, IMHS, ITHS, NP_G_IN, NE_G_IN, NP_A_IN, NE_A_IN, &
                                            IESTP, NSCOUE, IVSTP, NSCOUV, ICSTP, NSCOUC, IPSTP, IWSTP, NSCOUM, &
                                            IGEP, NSCOUGE, IGVP, NSCOUGV, IGCP, NSCOUGC, IGPP, IGWP, NSCOUGW
   end type

   !>
   !!
   !!
   type regrid_data
      real(ESMF_KIND_R8), pointer     :: src_fieldptr(:), mapped_fieldptr(:), unmapped_fieldptr(:), dst_maskptr(:)
      type(ESMF_Field)                :: src_datafield, dst_mask_field, dst_mapped_field, dst_unmapped_field
      type(ESMF_RouteHandle)          :: mapped_route_handle, unmapped_route_handle
   end type

contains

   !> \details As the name of this function suggests, this funciton creates a parallel
   !! ESMF_Mesh from meshdata object. This function should be called collectively by
   !! all PEs for the parallel mesh to be created. The function, extract_parallel_data_from_mesh()
   !! should be called prior to calling this function.
   !! \param the_data This the input meshdata object.
   !! \param out_esmf_mesh This is the ouput ESMF_Mesh object.
   subroutine create_parallel_esmf_mesh_from_meshdata(the_data, out_esmf_mesh)
      implicit none
      type(ESMF_Mesh), intent(out)                  :: out_esmf_mesh
      integer                                       :: rc
      type(meshdata), intent(in)                    :: the_data
      integer(ESMF_KIND_I4), parameter              :: dim1 = 2, spacedim = 2, NumND_per_El = 3

      if (.not. the_data%is_initialized) then
         call throw_fatal_error(__LINE__, __FILE__, &
            "The mesh is not initialized before calling "//&
            "create_parallel_esmf_mesh_from_meshdata")
      endif
      out_esmf_mesh = ESMF_MeshCreate(parametricDim=dim1, spatialDim=spacedim, &
                                      nodeIDs=the_data%NdIDs, nodeCoords=the_data%NdCoords, &
                                      nodeOwners=the_data%NdOwners, elementIDs=the_data%ElIDs, &
                                      elementTypes=the_data%ElTypes, elementConn=the_data%ElConnect, &
                                      rc=rc)
      CHECK_ERR_CODE(0, rc)
   end subroutine

   !> \details This function is similar to create_parallel_esmf_mesh_from_meshdata(), except that
   !! it creates a masked mesh. A masked mesh is used for example to exclude the interpolation onto
   !! some nodes, when using ESMF interpolation routines.
   !! \param in_meshdata This is the input meshdata object.
   !! \param mask_array This is an array of length NumNd (number of present nodes on this PE)
   !! which contains integer numbers. When we plan to exclude a group of nodes from interpolation,
   !! we use these mask values in the interpolation routines.
   !! \param out_masked_mesh This is the output masked ESMF_Mesh.
   subroutine create_masked_esmf_mesh_from_data(in_meshdata, mask_array, out_maked_esmf_mesh)
      implicit none
      type(ESMF_Mesh), intent(out)       :: out_maked_esmf_mesh
      type(meshdata), intent(in)         :: in_meshdata
      integer(ESMF_KIND_I4), intent(in)  :: mask_array(:)
      integer(ESMF_KIND_I4), parameter   :: dim1 = 2, spacedim = 2, NumND_per_El = 3
      integer(ESMF_KIND_I4)              :: rc, localPet, petCount

      call ESMF_VMGet(vm=in_meshdata%vm, localPet=localPet, petCount=petCount)
      if (localPet == 0 .and. .not. in_meshdata%is_initialized) then
         call throw_fatal_error(__LINE__, __FILE__,&
            "The mesh is not initialized before calling "//&
            "create_masked_esmf_mesh_from_data")
      endif
      out_maked_esmf_mesh = ESMF_MeshCreate(parametricDim=dim1, spatialDim=spacedim, &
                                            nodeIDs=in_meshdata%NdIDs, nodeCoords=in_meshdata%NdCoords, &
                                            nodeOwners=in_meshdata%NdOwners, elementIDs=in_meshdata%ElIDs, &
                                            elementTypes=in_meshdata%ElTypes, elementConn=in_meshdata%ElConnect, &
                                            nodeMask=mask_array, rc=rc)
      CHECK_ERR_CODE(localPet, rc)
   end subroutine create_masked_esmf_mesh_from_data

   !> @details Using the data available in <tt> fort.14, fort.18, partmesh.txt</tt> files
   !! this function extracts the scalars and arrays required for construction of a
   !! meshdata object.
   !! After calling this fucntion, one can call create_parallel_esmf_mesh_from_meshdata()
   !! or create_masked_esmf_mesh_from_data() to create an ESMF_Mesh.
   !! @param vm This is an ESMF_VM object, which will be used to obtain the \c localPE
   !! and \c peCount of the \c MPI_Communicator.
   !! @param global_fort14_dir This is the directory path (relative to the executable
   !! or an absolute path) which contains the global \c fort.14 file (not the fort.14
   !! after decomposition).
   !! @param the_data This is the output meshdata object.
   !!
   subroutine extract_parallel_data_from_mesh(vm, global_fort14_dir, the_data)
      implicit none
      type(ESMF_VM), intent(in)             :: vm
      type(meshdata), intent(inout)         :: the_data
      character(len=*), intent(in)          :: global_fort14_dir
      character(len=6)                      :: PE_ID, garbage1
      character(len=200)                    :: fort14_filename, fort18_filename, partmesh_filename
      integer(ESMF_KIND_I4)                 :: i1, j1, i_num, localPet, petCount, num_global_nodes, &
                                               garbage2, garbage3, iNum, dir_name_length, el_nd_ids(3)
      logical                               :: iOpen, iExist
      integer(ESMF_KIND_I4), allocatable    :: local_node_numbers(:), local_elem_numbers(:), node_owner(:), &
                                               global_2_local_node_map(:)
      integer(ESMF_KIND_I4), parameter      :: dim1 = 2, NumND_per_El = 3

      the_data%vm = vm
      call ESMF_VMGet(vm=vm, localPet=localPet, petCount=petCount)
      write (PE_ID, "(A,I4.4)") "PE", localPet
      !
      dir_name_length = len_trim(global_fort14_dir)
      allocate (character(len=dir_name_length)::the_data%dir_name)
      the_data%dir_name = trim(global_fort14_dir)
      !
      fort14_filename = trim(global_fort14_dir//PE_ID//"/fort.14")
      fort18_filename = trim(global_fort14_dir//PE_ID//"/fort.18")
      partmesh_filename = trim(global_fort14_dir//"/partmesh.txt")

      if (localPet == 0) then
         call show_message("Looking for global fort.14 in the directory: "//&
            the_data%dir_name//new_line("A"))
      endif
      inquire (FILE=fort14_filename, opened=iOpen, exist=iExist, number=iNum)
      if (.not. iExist) call throw_fatal_error("fort.14 file is not found. I tried to access file in: " &
                        //fort14_filename)

      open (unit=14, file=fort14_filename, form="FORMATTED", status="OLD", action="READ")
      open (unit=18, file=fort18_filename, form="FORMATTED", status="OLD", action="READ")
      open (unit=100, file=partmesh_filename, form="FORMATTED", status="OLD", action="READ")

      read (unit=14, fmt=*)
      read (unit=14, fmt=*) the_data%NumEl, the_data%NumNd
      allocate (the_data%NdIDs(the_data%NumNd))
      allocate (local_node_numbers(the_data%NumNd))
      allocate (the_data%ElIDs(the_data%NumEl))
      allocate (local_elem_numbers(the_data%NumEl))
      allocate (the_data%NdCoords(dim1*the_data%NumNd))
      allocate (the_data%bathymetry(the_data%NumNd))
      allocate (the_data%ElConnect(NumND_per_El*the_data%NumEl))
      allocate (the_data%NdOwners(the_data%NumNd))
      allocate (the_data%ElTypes(the_data%NumEl))

      read (unit=18, fmt=*)
      read (unit=18, fmt=*)
      read (unit=18, fmt=*) local_elem_numbers
      the_data%ElIDs = abs(local_elem_numbers)
      read (unit=18, fmt=*) garbage1, num_global_nodes, garbage2, garbage3
      if (garbage3 .NE. the_data%NumNd) call throw_fatal_error(__LINE__, __FILE__, &
                                             "issue with mesh partitioning.")
      read (unit=18, fmt=*) local_node_numbers
      the_data%NumOwnedND = 0
      do i1 = 1, the_data%NumNd, 1
         if (local_node_numbers(i1) > 0) then
            the_data%NumOwnedND = the_data%NumOwnedND + 1
         end if
      end do
      the_data%NdIDs = abs(local_node_numbers)
      allocate (node_owner(num_global_nodes))
      allocate (global_2_local_node_map(num_global_nodes))
      allocate (the_data%owned_to_present_nodes(the_data%NumOwnedND))
      global_2_local_node_map(:) = -1
      read (unit=100, fmt=*) node_owner

      do i1 = 1, the_data%NumNd, 1
         local_node_numbers(i1) = i1
         read (unit=14, fmt=*) j1, &
            the_data%NdCoords((i1 - 1)*dim1 + 1), &
            the_data%NdCoords((i1 - 1)*dim1 + 2), &
            the_data%bathymetry(i1)
         global_2_local_node_map(j1) = i1
      end do
      do i1 = 1, the_data%NumEl, 1
         read (unit=14, fmt=*) local_elem_numbers(i1), i_num, el_nd_ids(:)
         el_nd_ids(:) = global_2_local_node_map(el_nd_ids(:))
         do j1 = 1, 3
            the_data%ElConnect((i1 - 1)*NumND_per_El + j1) = el_nd_ids(j1)
         end do
      end do

      do i1 = 1, the_data%NumNd, 1
         the_data%NdOwners(i1) = node_owner(the_data%NdIDs(i1)) - 1
      end do

      j1 = 0
      do i1 = 1, the_data%NumNd, 1
         if (the_data%NdOwners(i1) == localPet) then
            j1 = j1 + 1
            the_data%owned_to_present_nodes(j1) = i1
         end if
      end do
      the_data%ElTypes = ESMF_MESHELEMTYPE_TRI
      the_data%is_initialized = .true.

      close (14)
      close (18)
      close (100)
   end subroutine extract_parallel_data_from_mesh

   !> \details This function writes the input meshdata object to a \c vtu file.
   !! The \c vtu file is in \c XML format. This function can be used for both parallel
   !! and serial mesh writing. If one uses this function for parallel write, the
   !! processing element with \c localPE=0 should also enter this function, otherwise
   !! the \c pvtu file will not be written. This function assumes that the \c vtu file
   !! which we want to write does not exist. If we want to add fields to the files which
   !! are created before, we have to call write_node_field_to_vtu() function. If the user
   !! wants to add more data fields to the created \c vtu file, the \c last_write parameter
   !! should be passed <tt>.false.</tt> so that the program do not close the file.
   !! By closing we mean writing the last three closing lines in the XML \c vtu files.
   !! However, if this is the last time we want to write on the same \c vtu file, we
   !! have to pass \c last_write equal to <tt>.true.</tt>
   !! \param the_data This is the input data for which we create the vtu file
   !! \param vtu_filename This is the name of the vtu file
   !! \param last_write This parameter indicates if this is the last time we want to
   !! write something to this \c vtu file.
   subroutine write_meshdata_to_vtu(the_data, vtu_filename, last_write, in_hot_data)
      implicit none
      type(meshdata), intent(in)            :: the_data
      character(len=*), intent(in)          :: vtu_filename
      integer(ESMF_KIND_I4)                 :: localPet, petCount
      logical, intent(in)                   :: last_write
      type(hotdata), optional, intent(in)   :: in_hot_data
      integer(ESMF_KIND_I4)                 :: i1, indent, offset_counter, rc, indent2
      integer(ESMF_KIND_I4), parameter      :: dim1 = 2, spacedim = 2, NumND_per_El = 3, vtk_triangle = 5
      indent = 0

      call ESMF_VMGet(vm=the_data%vm, localPet=localPet, petCount=petCount, rc=rc)
      if (rc .NE. ESMF_Success) then
         localPet = 0
         petCount = 1
      end if

      open (unit=1014, file=vtu_filename, form="FORMATTED", &
            status="UNKNOWN", action="WRITE")
      write (unit=1014, fmt="(A,A)") '<VTKFile type="UnstructuredGrid"', &
         ' version="0.1" byte_order="LittleEndian">'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         "<UnstructuredGrid>"
      indent = indent + 2
      write (unit=1014, fmt="(A,A,I0,A,I0,A)") repeat(" ", indent), &
         '<Piece NumberOfPoints="', the_data%NumNd, &
         '" NumberOfCells="', the_data%NumEl, '">'
      indent = indent + 2

      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<Points>'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
      indent = indent + 2
      do i1 = 1, the_data%NumNd, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,F0.4,' ',F0.4,' ',F0.4,' ')", advance='no') &
            repeat(" ", indent2), &
            the_data%NdCoords((i1 - 1)*dim1 + 1), &
            the_data%NdCoords((i1 - 1)*dim1 + 2), 0.0
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</Points>'

      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<Cells>'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Int32" Name="connectivity" Format="ascii">'
      indent = indent + 2
      do i1 = 1, the_data%NumEl, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,I0,' ',I0,' ',I0,' ')", advance='no') &
            repeat(" ", indent2), &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 1) - 1, &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 2) - 1, &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 3) - 1
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Int32" Name="offsets" Format="ascii">'
      indent = indent + 2
      offset_counter = 0
      do i1 = 1, the_data%NumEl, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         offset_counter = offset_counter + 3
         write (unit=1014, fmt="(A,I0)", advance='no') repeat(" ", indent2), &
            offset_counter
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Int32" Name="types" Format="ascii">'
      indent = indent + 2
      do i1 = 1, the_data%NumEl, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,I2)", advance='no') &
            repeat(" ", indent2), vtk_triangle
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</Cells>'

      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<CellData Scalars="scalars">'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Int32" Name="subdomain_id" NumberOfComponents="1" Format="ascii">'
      indent = indent + 2
      do i1 = 1, the_data%NumEl, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,I0)", advance='no') &
            repeat(" ", indent2), localPet
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'
      !
      if (present(in_hot_data)) then
          write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
             '<DataArray type="Int32" Name="NOFF" NumberOfComponents="1" Format="ascii">'
          indent = indent + 2
          do i1 = 1, the_data%NumEl, 1
             indent2 = 1
             if (i1 == 1) then; indent2 = indent; endif
             write (unit=1014, fmt="(A,I0)", advance='no') &
                repeat(" ", indent2), in_hot_data%NOFF(i1)
          end do
          write (unit=1014, fmt="(A)") ""
          indent = indent - 2
          write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
             '</DataArray>'
      end if
      !
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</CellData>'

      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<PointData Scalars="scalars">'
      indent = indent + 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Float32" Name="bathymetry" NumberOfComponents="1" Format="ascii">'
      indent = indent + 2
      do i1 = 1, the_data%NumNd, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,F0.4)", advance='no') &
            repeat(" ", indent2), &
            the_data%bathymetry(i1)
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'

      if (last_write) then
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</PointData>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</Piece>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</UnstructuredGrid>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</VTKFile>'
      end if
      close (1014)
   end subroutine

   !> \details This function writes the input array (\c field_array) and its name (\c field_name)
   !! to the vtu file (which should already exist and not closed). Refer to write_meshdata_to_vtu()
   !! to know more about opening vtu file and closing them. If the parameter \c last_write is true
   !! then we close this file and as such we should not write anything else on this file.
   subroutine write_node_field_to_vtu(field_array, field_name, vtu_filename, last_write)
      implicit none
      character(len=*), intent(in)   :: vtu_filename, field_name
      logical, intent(in)            :: last_write
      real(ESMF_KIND_R8), intent(in) :: field_array(:)
      integer(ESMF_KIND_I4)          :: i1, indent, indent2, num_recs
      open (unit=1014, file=vtu_filename, form='FORMATTED', &
            position='APPEND', status='OLD', action='WRITE')

      indent = 8
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Float32" Name="'//field_name//'" NumberOfComponents="1" Format="ascii">'
      indent = indent + 2
      num_recs = size(field_array)
      do i1 = 1, num_recs, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,F0.4)", advance='no') &
            repeat(" ", indent2), field_array(i1)
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'

      if (last_write) then
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</PointData>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</Piece>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</UnstructuredGrid>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</VTKFile>'
      end if
      close (1014)
   end subroutine

   !> \details This function writes the input array (\c field_array) and its name (\c field_name)
   !! to the vtu file (which should already exist and not closed). Refer to write_meshdata_to_vtu()
   !! to know more about opening vtu file and closing them. If the parameter \c last_write is true
   !! then we close this file and as such we should not write anything else on this file.
   subroutine write_int_node_field_to_vtu(field_array, field_name, vtu_filename, last_write)
      implicit none
      character(len=*), intent(in)      :: vtu_filename, field_name
      logical, intent(in)               :: last_write
      integer(ESMF_KIND_I4), intent(in) :: field_array(:)
      integer(ESMF_KIND_I4)             :: i1, indent, indent2, num_recs
      open (unit=1014, file=vtu_filename, form='FORMATTED', &
            position='APPEND', status='OLD', action='WRITE')

      indent = 8
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '<DataArray type="Int32" Name="'//field_name//'" NumberOfComponents="1" Format="ascii">'
      indent = indent + 2
      num_recs = size(field_array)
      do i1 = 1, num_recs, 1
         indent2 = 1
         if (i1 == 1) then; indent2 = indent; endif
         write (unit=1014, fmt="(A,I4)", advance='no') &
            repeat(" ", indent2), field_array(i1)
      end do
      write (unit=1014, fmt="(A)") ""
      indent = indent - 2
      write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
         '</DataArray>'

      if (last_write) then
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</PointData>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</Piece>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</UnstructuredGrid>'
         indent = indent - 2
         write (unit=1014, fmt="(A,A)") repeat(" ", indent), &
            '</VTKFile>'
      end if
      close (1014)
   end subroutine

   !> \details This function creates an object of type meshdata from the fort14 file given by
   !!fort14_filename. Unlike extract_parallel_data_from_mesh(), this function does not
   !! create a parallel meshdata, so it can be called by only one PE and the created meshdata
   !! object can later be used to create an ESMF_Mesh object.
   subroutine extract_global_data_from_fort14(the_data)
      implicit none
      type(meshdata), intent(inout)         :: the_data
      integer(ESMF_KIND_I4)                 :: i1, i_num, localPet, petCount, rc
      integer(ESMF_KIND_I4), parameter      :: dim1 = 2, spacedim = 2, NumND_per_El = 3
      logical                               :: iOpen, iExist

      call ESMF_VMGet(vm=the_data%vm, localPet=localPet, petCount=petCount, rc=rc)
      if (localPet == 0) then
         call show_message("looking for global fort.14 in the directory: "//the_data%dir_name)
         inquire (FILE=TRIM(the_data%dir_name)//"/fort.14", opened=iOpen, exist=iExist, number=i_num)
         if (.not. iExist) call throw_fatal_error("fort.14 file is not found."// &
            "I tried to access file in: "//TRIM(the_data%dir_name)//"/fort.14")
      endif
      open (unit=14, file=TRIM(the_data%dir_name)//"/fort.14", form='FORMATTED', status='OLD', action='READ')
      read (unit=14, fmt=*)
      read (unit=14, fmt=*) the_data%NumEl, the_data%NumNd
      allocate (the_data%NdIDs(the_data%NumNd))
      allocate (the_data%ElIDs(the_data%NumEl))
      allocate (the_data%NdCoords(dim1*the_data%NumNd))
      allocate (the_data%bathymetry(the_data%NumNd))
      allocate (the_data%ElConnect(NumND_per_El*the_data%NumEl))
      allocate (the_data%NdOwners(the_data%NumNd))
      allocate (the_data%ElTypes(the_data%NumEl))
      do i1 = 1, the_data%NumNd, 1
         read (unit=14, fmt=*) the_data%NdIDs(i1), &
            the_data%NdCoords((i1 - 1)*dim1 + 1), &
            the_data%NdCoords((i1 - 1)*dim1 + 2), &
            the_data%bathymetry(i1)
      end do
      do i1 = 1, the_data%NumEl, 1
         read (unit=14, fmt=*) the_data%ElIDs(i1), i_num, &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 1), &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 2), &
            the_data%ElConnect((i1 - 1)*NumND_per_El + 3)
      end do
      the_data%NdOwners = 0
      the_data%ElTypes = ESMF_MESHELEMTYPE_TRI
      close (14)
   end subroutine

   !> \details Given a local array in each PE i.e. \c fieldarray, we use MPI_Gather method to
   !! gather their elements into an array (\c out_fieldarray) in \c PE=root. For this
   !! process we use an ESMF_VM which is given to this function as an input. Since, MPI_Gather
   !! is collective this function should also be called collectively.
   subroutine gather_int_datafield_on_root(vm1, fieldarray, root, out_fieldarray, dir_name, partNd)
      implicit none

      type(ESMF_VM), intent(in)                         :: vm1
      integer(ESMF_KIND_I4), allocatable, intent(in)    :: fieldarray(:)
      integer(ESMF_KIND_I4), allocatable, intent(inout) :: out_fieldarray(:)
      integer(ESMF_KIND_I4), intent(in)                 :: root
      character(len=*), intent(in)                      :: dir_name
!Casey 210524
      integer(ESMF_KIND_I4), allocatable, intent(in), optional :: partNd(:)
      integer(ESMF_KIND_I4)                             :: send_count, localPet, petCount, num_total_nodes, &
                                                           i1, j1, k1, i_num, j_num1, rc, trash2, trash3
      integer(ESMF_KIND_I4), allocatable                :: recv_counts(:), gather_displs(:)
      integer(ESMF_KIND_I4), pointer                    :: temp_fieldarray(:)
      character(len=6)                                  :: PE_ID
      character(len=4)                                  :: trash1

      call ESMF_VMGet(vm=vm1, localPet=localPet, petCount=petCount, rc=rc)
      if (localPet == root) then
         if (.not. allocated(fieldarray) .or. .not. allocated(out_fieldarray)) then
            call throw_fatal_error( "The input arrays into "// &
               "gather_int_datafield_on_root should be allocated.")
         end if
      end if
      send_count = size(fieldarray)
      num_total_nodes = size(out_fieldarray)
      if (localPet == root) then
         allocate (recv_counts(petCount))
         allocate (gather_displs(petCount))
         gather_displs(1) = 0
         recv_counts = 0
!Casey 210524
#ifdef OLD
         open (unit=100, file=dir_name//"/partmesh.txt", form='FORMATTED', &
               status='OLD', action='READ')
#endif
         do i1 = 1, num_total_nodes, 1
#ifdef OLD
            read (100, *) i_num
#else
            i_num = partNd(i1) + 1
#endif
            recv_counts(i_num) = recv_counts(i_num) + 1
         end do
         do i1 = 2, petCount, 1
            gather_displs(i1) = gather_displs(i1 - 1) + recv_counts(i1 - 1)
         end do
         allocate (temp_fieldarray(num_total_nodes))
#ifdef OLD
         close (100)
#endif
      end if
      !
      call MPI_Gatherv(fieldarray, send_count, MPI_INTEGER, &
                       temp_fieldarray, recv_counts, gather_displs, MPI_INTEGER, &
                       root, MPI_COMM_WORLD, rc)
      !
      CHECK_ERR_CODE(localPet, rc)
      if (localPet == 0) call show_message("    gathered on root with code: 0"//&
                                           new_line('A'))
      if (localPet == 0) then
         do i1 = 1, petCount, 1
!Casey 210524
#ifdef OLD
            write (PE_ID, "(A,I4.4)") 'PE', i1 - 1
            open (unit=18, file=dir_name//PE_ID//"/fort.18", form='FORMATTED', &
                  status='OLD', action='READ')
            read (unit=18, fmt=*)
            read (unit=18, fmt=*) trash1, trash2, trash3, i_num
            do j1 = 1, i_num, 1
               read (unit=18, fmt=*)
            end do
            k1 = 0
            read (unit=18, fmt=*) trash1, trash2, trash3, i_num
            do j1 = 1, i_num, 1
               read (unit=18, fmt=*) j_num1
               if (j_num1 > 0) then
                  k1 = k1 + 1
                  out_fieldarray(j_num1) = temp_fieldarray(gather_displs(i1) + k1)
               end if
            end do
#else
            k1 = 0
            do j1 = 1, num_total_nodes
              if( partNd(j1).eq.(i1-1) )then
                k1 = k1 + 1
                out_fieldarray(j1) = temp_fieldarray(gather_displs(i1) + k1)
              endif
            enddo
#endif
         end do
      end if
   end subroutine

   !> \details Given a local array in each PE i.e. \c fieldarray, we use MPI_Gather method to
   !! gather their elements into an array (\c out_fieldarray) in \c PE=root. For this
   !! process we use an ESMF_VM which is given to this function as an input. Since, MPI_Gather
   !! is collective this function should also be called collectively.
   subroutine gather_datafield_on_root(vm1, fieldarray, root, out_fieldarray, dir_name, partNd)
      implicit none

      type(ESMF_VM), intent(in)                       :: vm1
      real(ESMF_KIND_R8), allocatable, intent(in)     :: fieldarray(:)
      real(ESMF_KIND_R8), allocatable, intent(inout)  :: out_fieldarray(:)
      integer(ESMF_KIND_I4), intent(in)               :: root
      character(len=*), intent(in)                    :: dir_name
!Casey 210524
      integer(ESMF_KIND_I4),allocatable,intent(in),optional :: partNd(:)
      integer(ESMF_KIND_I4)                           :: send_count, localPet, petCount, num_total_nodes, &
                                                         i1, j1, k1, i_num, j_num1, rc, trash2, trash3
      integer(ESMF_KIND_I4), allocatable              :: recv_counts(:), gather_displs(:)
      real(ESMF_KIND_R8), allocatable                 :: temp_fieldarray(:)
      character(len=6)                                :: PE_ID
      character(len=4)                                :: trash1

      call ESMF_VMGet(vm=vm1, localPet=localPet, petCount=petCount, rc=rc)
      if (localPet == root) then
         if (.not. allocated(fieldarray) .or. .not. allocated(out_fieldarray)) then
            call throw_fatal_error("The input arrays into"//&
               "gather_datafield_on_root should be allocated.")
         end if
      endif
      send_count = size(fieldarray)
      num_total_nodes = size(out_fieldarray)
      if (localPet == root) then
         allocate (recv_counts(petCount))
         allocate (gather_displs(petCount))
         gather_displs(1) = 0
         recv_counts = 0
!Casey 210524
#ifdef OLD
         open (unit=100, file=dir_name//"/partmesh.txt", form='FORMATTED', &
               status='OLD', action='READ')
#endif
         do i1 = 1, num_total_nodes, 1
#ifdef OLD
            read (100, *) i_num
#else
            i_num = partNd(i1) + 1
#endif
            recv_counts(i_num) = recv_counts(i_num) + 1
         end do
         do i1 = 2, petCount, 1
            gather_displs(i1) = gather_displs(i1 - 1) + recv_counts(i1 - 1)
         end do
         allocate (temp_fieldarray(num_total_nodes))
#ifdef OLD
         close (100)
#endif
      end if

      call MPI_Gatherv(fieldarray, send_count, MPI_DOUBLE_PRECISION, &
                       temp_fieldarray, recv_counts, gather_displs, MPI_DOUBLE_PRECISION, &
                       0, MPI_COMM_WORLD, rc)
      CHECK_ERR_CODE(localPet, rc)
      if (localPet == 0) call show_message("    gathered on root with code: 0"//&
                                           new_line('A'))
      if (localPet == 0) then
         do i1 = 1, petCount, 1
!Casey 210524
#ifdef OLD
            write (PE_ID, "(A,I4.4)") 'PE', i1 - 1
            open (unit=18, file=dir_name//PE_ID//"/fort.18", form='FORMATTED', &
                  status='OLD', action='READ')
            read (unit=18, fmt=*)
            read (unit=18, fmt=*) trash1, trash2, trash3, i_num
            do j1 = 1, i_num, 1
               read (unit=18, fmt=*)
            end do
            k1 = 0
            read (unit=18, fmt=*) trash1, trash2, trash3, i_num
            do j1 = 1, i_num, 1
               read (unit=18, fmt=*) j_num1
               if (j_num1 > 0) then
                  k1 = k1 + 1
                  out_fieldarray(j_num1) = temp_fieldarray(gather_displs(i1) + k1)
               end if
            end do
#else
            k1 = 0
            do j1 = 1, num_total_nodes
              if( partNd(j1).eq.(i1-1) )then
                k1 = k1 + 1
                out_fieldarray(j1) = temp_fieldarray(gather_displs(i1) + k1)
              endif
            enddo
#endif
         end do
      end if
   end subroutine

   !>
   !!
   !!
   subroutine gather_dst_nodal_hotdata_on_root(localized_hotdata, global_hotdata, &
                                               localized_meshdata, global_meshdata, root)
      implicit none
      real(ESMF_KIND_R8), allocatable     :: localfieldarray_ptr(:), globalfieldarray_ptr(:)
      integer(ESMF_KIND_I4), allocatable  :: int_localfieldarray_ptr(:), int_globalfieldarray_ptr(:)
      type(hotdata), intent(in)           :: localized_hotdata
      type(hotdata), intent(inout)        :: global_hotdata
      type(meshdata), intent(in)          :: localized_meshdata
      type(meshdata), intent(in)          :: global_meshdata
      integer(ESMF_KIND_I4)               :: root, localPet, petCount, rc
!Casey 210524
      integer(ESMF_KIND_I4),allocatable   :: partNd(:)
      integer(ESMF_KIND_I4) :: gvv
      integer(ESMF_KIND_I4) :: ie
      integer(ESMF_KIND_I4) :: ivv
      integer(ESMF_KIND_I4) :: owner

      call ESMF_VMGet(vm=localized_meshdata%vm, localPet=localPet, petCount=petCount, rc=rc)
      allocate (localfieldarray_ptr(localized_meshdata%NumOwnedNd))
      allocate (int_localfieldarray_ptr(localized_meshdata%NumOwnedNd))
      if (localPet == root) then
         allocate (globalfieldarray_ptr(global_meshdata%NumNd))
         allocate (int_globalfieldarray_ptr(global_meshdata%NumNd))
!Casey 210524
         allocate(partNd(global_meshdata%NumNd))
         partNd = petCount - 1
         do ie=1,global_meshdata%NumEl
           owner = floor( real(ie-1)                                   &
     &       / real(floor(real(global_meshdata%NumEl)/real(petCount))) )
           do ivv=1,3
             gvv = global_meshdata%ElConnect((ie-1)*3+ivv)
             if( partNd(gvv).gt.owner )then
               partNd(gvv) = owner
             endif
           enddo
         enddo
      endif

      if (localPet == root) call show_message("Gathering dst eta1 on root:")
      localfieldarray_ptr = localized_hotdata%ETA1
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name, &
!Casey 210524
     &                              partNd)
      if (localPet == root) global_hotdata%ETA1 = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering dst eta2 on root:")
      localfieldarray_ptr = localized_hotdata%ETA2
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name, &
!Casey 210524
     &                              partNd)
      if (localPet == root) global_hotdata%ETA2 = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering dst etaDisc on root:")
      localfieldarray_ptr = localized_hotdata%ETADisc
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name, &
!Casey 210524
     &                              partNd)
      if (localPet == root) global_hotdata%ETADisc = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering dst UU2 on root:")
      localfieldarray_ptr = localized_hotdata%UU2
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name, &
!Casey 210524
     &                              partNd)
      if (localPet == root) global_hotdata%UU2 = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering dst VV2 on root:")
      localfieldarray_ptr = localized_hotdata%VV2
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name, &
!Casey 210524
     &                              partNd)
      if (localPet == root) global_hotdata%VV2 = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering dst CH1 on root:")
      localfieldarray_ptr = localized_hotdata%CH1
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name, &
!Casey 210524
     &                              partNd)
      if (localPet == root) global_hotdata%CH1 = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering dst NODECODE on root:")
      int_localfieldarray_ptr = localized_hotdata%NNODECODE
      call gather_int_datafield_on_root(localized_meshdata%vm, int_localfieldarray_ptr, &
                                        root, int_globalfieldarray_ptr, localized_meshdata%dir_name, &
!Casey 210524
     &                                  partNd)
      if (localPet == root) global_hotdata%NNODECODE = int_globalfieldarray_ptr
   end subroutine

   !>
   !!
   !!
   subroutine gather_src_nodal_hotdata_on_root(localized_hotdata, global_hotdata, &
                                               localized_meshdata, global_meshdata, root)
      implicit none
      real(ESMF_KIND_R8), allocatable     :: localfieldarray_ptr(:), globalfieldarray_ptr(:)
      integer(ESMF_KIND_I4), allocatable  :: int_localfieldarray_ptr(:), int_globalfieldarray_ptr(:)
      type(hotdata), intent(in)           :: localized_hotdata
      type(hotdata), intent(inout)        :: global_hotdata
      type(meshdata), intent(in)          :: localized_meshdata
      type(meshdata), intent(in)          :: global_meshdata
      integer(ESMF_KIND_I4)               :: root, localPet, petCount, rc, i1

      call ESMF_VMGet(vm=localized_meshdata%vm, localPet=localPet, petCount=petCount, rc=rc)
      allocate (localfieldarray_ptr(localized_meshdata%NumOwnedNd))
      allocate (int_localfieldarray_ptr(localized_meshdata%NumOwnedNd))
      if (localPet == root) then
         allocate (globalfieldarray_ptr(global_meshdata%NumNd))
         allocate (int_globalfieldarray_ptr(global_meshdata%NumNd))
      end if

      if (localPet == root) call show_message("Gathering src eta1 on root:")
      localfieldarray_ptr(:) = localized_hotdata%ETA1(localized_meshdata%owned_to_present_nodes(:))
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name)
      if (localPet == root) global_hotdata%ETA1 = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering src eta2 on root:")
      localfieldarray_ptr(:) = localized_hotdata%ETA2(localized_meshdata%owned_to_present_nodes(:))
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name)
      if (localPet == root) global_hotdata%ETA2 = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering src etaDisc on root:")
      localfieldarray_ptr(:) = localized_hotdata%ETADisc(localized_meshdata%owned_to_present_nodes(:))
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name)
      if (localPet == root) global_hotdata%ETADisc = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering src UU2 on root:")
      localfieldarray_ptr(:) = localized_hotdata%UU2(localized_meshdata%owned_to_present_nodes(:))
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name)
      if (localPet == root) global_hotdata%UU2 = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering src VV2 on root:")
      localfieldarray_ptr(:) = localized_hotdata%VV2(localized_meshdata%owned_to_present_nodes(:))
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name)
      if (localPet == root) global_hotdata%VV2 = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering src CH1 on root:")
      localfieldarray_ptr(:) = localized_hotdata%CH1(localized_meshdata%owned_to_present_nodes(:))
      call gather_datafield_on_root(localized_meshdata%vm, localfieldarray_ptr, &
                                    root, globalfieldarray_ptr, localized_meshdata%dir_name)
      if (localPet == root) global_hotdata%CH1 = globalfieldarray_ptr
      !
      if (localPet == root) call show_message("Gathering src NODECODE on root:")
      int_localfieldarray_ptr(:) = localized_hotdata%NNODECODE(localized_meshdata%owned_to_present_nodes(:))
      call gather_int_datafield_on_root(localized_meshdata%vm, int_localfieldarray_ptr, &
                                        root, int_globalfieldarray_ptr, localized_meshdata%dir_name)
      if (localPet == root) global_hotdata%NNODECODE = int_globalfieldarray_ptr

   end subroutine

   !>
   !!
   !!
   subroutine extract_hotdata_from_parallel_binary_fort_67(the_meshdata, the_hotdata, global_fort14_dir, write_ascii)
      implicit none
      type(meshdata), intent(in)   :: the_meshdata
      type(hotdata), intent(out)   :: the_hotdata
      character(len=*), intent(in) :: global_fort14_dir
      logical, intent(in)          :: write_ascii
      integer(ESMF_KIND_I4)        :: i1, localPet, petCount, rc, ihotstp, iNum, io_stat
      character(len=6)             :: PE_ID
      character(len=200)           :: fort67_filename, fort67_ascii_filename
      logical                      :: iExist, iOpen

      call ESMF_VMGet(vm=the_meshdata%vm, localPet=localPet, petCount=petCount, rc=rc)
      call allocate_hotdata(the_hotdata, the_meshdata)

      write (PE_ID, "(A,I4.4)") "PE", localPet
      fort67_filename = trim(global_fort14_dir//PE_ID//"/fort.67")
      inquire (FILE=fort67_filename, opened=iOpen, exist=iExist, number=iNum)
      if (.not. iExist) then
         call throw_fatal_error("parallel binary fort.67 file is not found. "//&
                                "I tried to access file in: "//fort67_filename)
      endif

      open (unit=67, file=fort67_filename, action='READ', &
            access='DIRECT', recl=8, iostat=rc, status='OLD')
      ihotstp = 1
      read (unit=67, REC=ihotstp, iostat=io_stat) the_hotdata%InputFileFmtVn;
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp, iostat=io_stat) the_hotdata%IMHS;
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp, iostat=io_stat) the_hotdata%TimeLoc;
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp, iostat=io_stat) the_hotdata%ITHS;
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp, iostat=io_stat) the_hotdata%NP_G_IN;
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp, iostat=io_stat) the_hotdata%NE_G_IN;
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp, iostat=io_stat) the_hotdata%NP_A_IN;
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp, iostat=io_stat) the_hotdata%NE_A_IN;
      ihotstp = ihotstp + 1

      do i1 = 1, the_meshdata%NumNd, 1
         read (unit=67, REC=ihotstp) the_hotdata%ETA1(i1)
         ihotstp = ihotstp + 1
      end do
      do i1 = 1, the_meshdata%NumNd, 1
         read (unit=67, REC=ihotstp) the_hotdata%ETA2(i1)
         ihotstp = ihotstp + 1
      end do
      do i1 = 1, the_meshdata%NumNd, 1
         read (unit=67, REC=ihotstp) the_hotdata%ETADisc(i1)
         ihotstp = ihotstp + 1
      end do
      do i1 = 1, the_meshdata%NumNd, 1
         read (unit=67, REC=ihotstp) the_hotdata%UU2(i1)
         ihotstp = ihotstp + 1
      end do
      do i1 = 1, the_meshdata%NumNd, 1
         read (unit=67, REC=ihotstp) the_hotdata%VV2(i1)
         ihotstp = ihotstp + 1
      end do
      if (the_hotdata%IMHS .EQ. 10) then
         do i1 = 1, the_meshdata%NumNd, 1
            read (unit=67, REC=ihotstp) the_hotdata%CH1(i1)
            ihotstp = ihotstp + 1
         end do
      end if
      do i1 = 1, the_meshdata%NumNd, 1
         read (unit=67, REC=ihotstp) the_hotdata%NNODECODE(i1)
         the_hotdata%realNODECODE(i1) = dble(the_hotdata%NNODECODE(i1))
         ihotstp = ihotstp + 1
      end do
      do i1 = 1, the_meshdata%NumEl, 1
         read (unit=67, REC=ihotstp) the_hotdata%NOFF(i1)
         ihotstp = ihotstp + 1
      end do

      read (unit=67, REC=ihotstp) the_hotdata%IESTP
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%NSCOUE
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%IVSTP
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%NSCOUV
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%ICSTP
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%NSCOUC
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%IPSTP
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%IWSTP
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%NSCOUM
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%IGEP
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%NSCOUGE
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%IGVP
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%NSCOUGV
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%IGCP
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%NSCOUGC
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%IGPP
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%IGWP
      ihotstp = ihotstp + 1
      read (unit=67, REC=ihotstp) the_hotdata%NSCOUGW
      ihotstp = ihotstp + 1
      close (67)

      if (write_ascii) then
         fort67_ascii_filename = trim(global_fort14_dir//PE_ID//"/fort.67")//".txt"
         open (unit=670, file=fort67_ascii_filename, form='FORMATTED', &
               action='WRITE', iostat=rc)
         write (unit=670, fmt=*) the_hotdata%InputFileFmtVn
         write (unit=670, fmt=*) the_hotdata%IMHS
         write (unit=670, fmt=*) the_hotdata%TimeLoc
         write (unit=670, fmt=*) the_hotdata%ITHS
         write (unit=670, fmt=*) the_hotdata%NP_G_IN
         write (unit=670, fmt=*) the_hotdata%NE_G_IN
         write (unit=670, fmt=*) the_hotdata%NP_A_IN
         write (unit=670, fmt=*) the_hotdata%NE_A_IN
         write (unit=670, fmt=*) the_hotdata%ETA1
         write (unit=670, fmt=*) the_hotdata%ETA2
         write (unit=670, fmt=*) the_hotdata%ETADisc
         write (unit=670, fmt=*) the_hotdata%UU2
         write (unit=670, fmt=*) the_hotdata%VV2
         write (unit=670, fmt=*) the_hotdata%CH1
         write (unit=670, fmt=*) the_hotdata%NNODECODE
         write (unit=670, fmt=*) the_hotdata%NOFF
         write (unit=670, fmt=*) the_hotdata%IESTP, the_hotdata%NSCOUE, the_hotdata%IVSTP, &
            the_hotdata%NSCOUV, the_hotdata%ICSTP, the_hotdata%NSCOUC, the_hotdata%IPSTP, &
            the_hotdata%IWSTP, the_hotdata%NSCOUM, the_hotdata%IGEP, the_hotdata%NSCOUGE, &
            the_hotdata%IGVP, the_hotdata%NSCOUGV, the_hotdata%IGCP, the_hotdata%NSCOUGC, &
            the_hotdata%IGPP, the_hotdata%IGWP, the_hotdata%NSCOUGW
         close (670)
      end if
      return
   end subroutine

   !>
   !!
   !!
   subroutine write_serial_hotfile_to_fort_67(the_meshdata, the_hotdata, global_fort14_dir, write_ascii)
      implicit none
      type(meshdata), intent(in)   :: the_meshdata
      type(hotdata), intent(in)    :: the_hotdata
      character(len=*), intent(in) :: global_fort14_dir
      logical, intent(in)          :: write_ascii
      integer(ESMF_KIND_I4)        :: i1, rc, ihotstp
      character(len=200)           :: fort67_filename, fort67_ascii_filename

      fort67_filename = trim(global_fort14_dir//"/fort.67")
      open (unit=67, file=fort67_filename, action='WRITE', &
            access='DIRECT', recl=8, iostat=rc, status='UNKNOWN')

      ihotstp = 1
      write (unit=67, REC=ihotstp) the_hotdata%InputFileFmtVn;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IMHS;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%TimeLoc;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%ITHS;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NP_G_IN;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NE_G_IN;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NP_A_IN;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NE_A_IN;
      ihotstp = ihotstp + 1

      do i1 = 1, the_meshdata%NumNd, 1
         write (unit=67, REC=ihotstp) the_hotdata%ETA1(i1)
         ihotstp = ihotstp + 1
      end do
      do i1 = 1, the_meshdata%NumNd, 1
         write (unit=67, REC=ihotstp) the_hotdata%ETA2(i1)
         ihotstp = ihotstp + 1
      end do
      do i1 = 1, the_meshdata%NumNd, 1
         write (unit=67, REC=ihotstp) the_hotdata%ETADisc(i1)
         ihotstp = ihotstp + 1
      end do
      do i1 = 1, the_meshdata%NumNd, 1
         write (unit=67, REC=ihotstp) the_hotdata%UU2(i1)
         ihotstp = ihotstp + 1
      end do
      do i1 = 1, the_meshdata%NumNd, 1
         write (unit=67, REC=ihotstp) the_hotdata%VV2(i1)
         ihotstp = ihotstp + 1
      end do
      if (the_hotdata%IMHS .EQ. 10) then
         do i1 = 1, the_meshdata%NumNd, 1
            write (unit=67, REC=ihotstp) the_hotdata%CH1(i1)
            ihotstp = ihotstp + 1
         end do
      end if
      do i1 = 1, the_meshdata%NumNd, 1
         write (unit=67, REC=ihotstp) the_hotdata%NNODECODE(i1)
         ihotstp = ihotstp + 1
      end do
      do i1 = 1, the_meshdata%NumEl, 1
         write (unit=67, REC=ihotstp) the_hotdata%NOFF(i1)
         ihotstp = ihotstp + 1
      end do

      write (unit=67, REC=ihotstp) the_hotdata%IESTP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUE
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IVSTP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUV
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%ICSTP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUC
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IPSTP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IWSTP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUM
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IGEP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUGE
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IGVP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUGV
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IGCP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUGC
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IGPP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IGWP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUGW
      ihotstp = ihotstp + 1
      close (67)

      if (write_ascii) then
         fort67_ascii_filename = trim(global_fort14_dir//"/fort.67")//".txt"
         open (unit=670, file=fort67_ascii_filename, form='FORMATTED', &
               action='WRITE', iostat=rc)
         write (unit=670, fmt=*) the_hotdata%InputFileFmtVn
         write (unit=670, fmt=*) the_hotdata%IMHS
         write (unit=670, fmt=*) the_hotdata%TimeLoc
         write (unit=670, fmt=*) the_hotdata%ITHS
         write (unit=670, fmt=*) the_hotdata%NP_G_IN
         write (unit=670, fmt=*) the_hotdata%NE_G_IN
         write (unit=670, fmt=*) the_hotdata%NP_A_IN
         write (unit=670, fmt=*) the_hotdata%NE_A_IN
         write (unit=670, fmt=*) the_hotdata%ETA1
         write (unit=670, fmt=*) the_hotdata%ETA2
         write (unit=670, fmt=*) the_hotdata%ETADisc
         write (unit=670, fmt=*) the_hotdata%UU2
         write (unit=670, fmt=*) the_hotdata%VV2
         write (unit=670, fmt=*) the_hotdata%CH1
         write (unit=670, fmt=*) the_hotdata%NNODECODE
         write (unit=670, fmt=*) the_hotdata%NOFF
         write (unit=670, fmt=*) the_hotdata%IESTP, the_hotdata%NSCOUE, the_hotdata%IVSTP, &
            the_hotdata%NSCOUV, the_hotdata%ICSTP, the_hotdata%NSCOUC, the_hotdata%IPSTP, &
            the_hotdata%IWSTP, the_hotdata%NSCOUM, the_hotdata%IGEP, the_hotdata%NSCOUGE, &
            the_hotdata%IGVP, the_hotdata%NSCOUGV, the_hotdata%IGCP, the_hotdata%NSCOUGC, &
            the_hotdata%IGPP, the_hotdata%IGWP, the_hotdata%NSCOUGW
         close (670)
      end if
   end subroutine

   !>
   !!
   !!
   subroutine fast_write_serial_hotfile_to_fort_67(the_meshdata, the_hotdata, global_fort14_dir)
      implicit none
      type(meshdata), intent(in)   :: the_meshdata
      type(hotdata), intent(in)    :: the_hotdata
      character(len=*), intent(in) :: global_fort14_dir
      integer(ESMF_KIND_I4)        :: i1, rc, ihotstp
      character(len=200)           :: fort67_filename, fort67_ascii_filename

      fort67_filename = trim(global_fort14_dir//"/fort.67")
      open (unit=67, file=fort67_filename, action='WRITE', &
            access='DIRECT', recl=8, iostat=rc, status='UNKNOWN')

      ihotstp = 1
      write (unit=67, REC=ihotstp) the_hotdata%InputFileFmtVn;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IMHS;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%TimeLoc;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%ITHS;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NP_G_IN;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NE_G_IN;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NP_A_IN;
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NE_A_IN;
      ihotstp = ihotstp + 1

      write (unit=67, rec=ihotstp) the_hotdata%ETA1(1:the_meshdata%NumNd)
      ihotstp = ihotstp + the_meshdata%NumNd
      write (unit=67) the_hotdata%ETA2(1:the_meshdata%NumNd)
      ihotstp = ihotstp + the_meshdata%NumNd
      write (unit=67) the_hotdata%ETADisc(1:the_meshdata%NumNd)
      ihotstp = ihotstp + the_meshdata%NumNd
      write (unit=67) the_hotdata%UU2(1:the_meshdata%NumNd)
      ihotstp = ihotstp + the_meshdata%NumNd
      write (unit=67) the_hotdata%VV2(1:the_meshdata%NumNd)
      ihotstp = ihotstp + the_meshdata%NumNd
      if (the_hotdata%IMHS .EQ. 10) then
         write (unit=67) the_hotdata%CH1(1:the_meshdata%NumNd)
         ihotstp = ihotstp + the_meshdata%NumNd
      end if
      write (unit=67) the_hotdata%NNODECODE(1:the_meshdata%NumNd)
      ihotstp = ihotstp + the_meshdata%NumNd
      write (unit=67) the_hotdata%NOFF(1:the_meshdata%NumEl)
      ihotstp = ihotstp + the_meshdata%NumEl

      write (unit=67, REC=ihotstp) the_hotdata%IESTP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUE
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IVSTP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUV
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%ICSTP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUC
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IPSTP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IWSTP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUM
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IGEP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUGE
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IGVP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUGV
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IGCP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUGC
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IGPP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%IGWP
      ihotstp = ihotstp + 1
      write (unit=67, REC=ihotstp) the_hotdata%NSCOUGW
      ihotstp = ihotstp + 1
      close (67)
   end subroutine

   subroutine allocate_hotdata(the_hotdata, the_meshdata)
      implicit none
      type(hotdata), intent(inout)    :: the_hotdata
      type(meshdata), intent(in)      :: the_meshdata
      allocate (the_hotdata%ETA1(the_meshdata%NumNd))
      allocate (the_hotdata%ETA2(the_meshdata%NumNd))
      allocate (the_hotdata%ETADisc(the_meshdata%NumNd))
      allocate (the_hotdata%UU2(the_meshdata%NumNd))
      allocate (the_hotdata%VV2(the_meshdata%NumNd))
      allocate (the_hotdata%CH1(the_meshdata%NumNd))
      allocate (the_hotdata%realNODECODE(the_meshdata%NumNd))
      allocate (the_hotdata%NNODECODE(the_meshdata%NumNd))
      allocate (the_hotdata%NOFF(the_meshdata%NumEl))
   end subroutine

   subroutine regrid_datafield_of_present_nodes(the_regrid_data, src_data, dst_data, src_array_of_present_nodes, print_unmapped)
      implicit none
      type(regrid_data), intent(inout)    :: the_regrid_data
      type(meshdata)                      :: src_data, dst_data
      real(ESMF_KIND_R8), intent(in)      :: src_array_of_present_nodes(:)
      integer(ESMF_KIND_I4)               :: i1, localPet, petCount, rc
      character(len=4)                    :: rc_str
      logical, optional                   :: print_unmapped
      logical                             :: pls_print
      if (present(print_unmapped)) then
         pls_print = print_unmapped
      else
         pls_print = .false.
      endif

      call ESMF_VMGet(vm=src_data%vm, localPet=localPet, petCount=petCount, rc=rc)
      do i1 = 1, src_data%NumOwnedNd, 1
         the_regrid_data%src_fieldptr(i1) = src_array_of_present_nodes(src_data%owned_to_present_nodes(i1))
      end do
      call ESMF_FieldRegrid(srcField=the_regrid_data%src_datafield, &
                            dstField=the_regrid_data%dst_mapped_field, &
                            routeHandle=the_regrid_data%mapped_route_handle, rc=rc)
      write(rc_str, "(I4)") rc
      if (localPet == 0) call show_message("    mapped regriding done with code: "//rc_str)
      call ESMF_FieldRegrid(srcField=the_regrid_data%src_datafield, &
                            dstField=the_regrid_data%dst_unmapped_field, &
                            routeHandle=the_regrid_data%unmapped_route_handle, rc=rc)
      write(rc_str, "(I4)") rc
      if (localPet == 0) call show_message("    unmapped regriding done with code: "//&
                                           rc_str//new_line('A'))
      do i1 = 1, dst_data%NumOwnedND, 1
         if (abs(the_regrid_data%dst_maskptr(i1)) < 1.d-8) then
            the_regrid_data%mapped_fieldptr(i1) = the_regrid_data%unmapped_fieldptr(i1)
         end if
      end do
      call MPI_Barrier(MPI_COMM_WORLD, rc)
   end subroutine

   !> \details This function simply deallocates the arrays created in the \c meshdata object
   !! creation steps.
   subroutine destroy_meshdata(the_data)
      implicit none
      type(meshdata), intent(inout) :: the_data
      if(allocated(the_data%NdIDs)) deallocate (the_data%NdIDs)
      if(allocated(the_data%ElIDs)) deallocate (the_data%ElIDs)
      if(allocated(the_data%NdCoords)) deallocate (the_data%NdCoords)
      if(allocated(the_data%bathymetry)) deallocate (the_data%bathymetry)
      if(allocated(the_data%ElConnect)) deallocate (the_data%ElConnect)
      if(allocated(the_data%NdOwners)) deallocate (the_data%NdOwners)
      if(allocated(the_data%ElTypes)) deallocate (the_data%ElTypes)
   end subroutine

   subroutine destroy_hotdata(the_data)
      implicit none
      type(hotdata), intent(inout)    :: the_data
   end subroutine

   subroutine destroy_regrid_data(the_regrid_data)
      implicit none
      type(regrid_data), intent(inout)    :: the_regrid_data
      call ESMF_FieldRegridRelease(the_regrid_data%mapped_route_handle)
      call ESMF_FieldRegridRelease(the_regrid_data%unmapped_route_handle)
      call ESMF_FieldDestroy(the_regrid_data%src_datafield)
      call ESMF_FieldDestroy(the_regrid_data%dst_mask_field)
      call ESMF_FieldDestroy(the_regrid_data%dst_mapped_field)
      call ESMF_FieldDestroy(the_regrid_data%dst_unmapped_field)
   end subroutine

end module adcirpolate

