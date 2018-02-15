!> \mainpage
!! ## Introduction
!! \ref adcirc_interpolation "ADCIRC interpolation module" provides required types and
!! functions to interpolate data between two ADCIRC meshes.
!! We use this module to first create an object of type \ref adcirc_interpolation::meshdata
!! "meshdata" from ADCIRC input files, and then use the created meshdata object to construct an
!! ESMF_Mesh object. When we have two ESMF_Mesh objects, we can use them to interpolate
!! data between them. All of this process can either be done
!! sequentially or in parallel.
!! ### Sequential interpolation
!! For sequential interpolation we use the following procedure:
!!   1. Use the subroutine \ref adcirc_interpolation::extract_global_data_from_fort14()
!!      "extract_global_data_from_fort14()" to create two objects of type
!!      \ref adcirc_interpolation::meshdata "meshdata". One of these objects is used
!!      as source mesh where we want to interpolate from, and the other is used as
!!      destination mesh where we want to interpolate to.
!!   2. Create two ESMF_Mesh objects from the above \ref adcirc_interpolation::meshdata
!!      "meshdata"'s, using \ref adcirc_interpolation::create_parallel_esmf_mesh_from_meshdata
!!      "create_parallel_esmf_mesh_from_meshdata" subroutine.
!!   3. Create two \c ESMF_Field 's on the above \c ESMF_Mesh 's.
!!   4. Use ESMF library to interpolate data from the source mesh to destination mesh for
!!      those nodes that we have enough data, and extrapolate for those nodes that we do
!!      not have enough data.
!!
!! ### Parallel interpolation
!! The process is very much similar to the sequential case except in the first step we
!! use the function \ref adcirc_interpolation::extract_parallel_data_from_mesh
!! "extract_parallel_data_from_mesh" to extract mesh data from the files <tt> fort.14, fort.18,
!! partmesh.txt</tt>.
!! ## Basic Usage
!! Here, we present an example of how the module adcirc_interpolation should be used.
!! Consider the following two meshes of east coast. We have decomposed each of these meshes
!! into 4 subdomain, which are shown with different colors here. The subdomain decomposition
!! is done by adcprep. The coarse mesh is used as the source mesh and the fine mesh is our
!! destination mesh. We want to interpolate the bathymetry from the source mesh to the
!! destination mesh.
!! <img src="Images/coarse_mesh_subdomains.png" width=600em />
!! <div style="text-align:center; font-size:150%;">The decomposed coarse mesh, which is used as the source mesh.</div>
!! <img src="Images/fine_mesh_subdomains.png" width=600em />
!! <div style="text-align:center; font-size:150%;">The decomposed fine mesh, which is used as the destination mesh.</div>
!!
#include "c_preprocessor.h"

program main

   use ESMF
   use MPI
   use adcirpolate
   use wetdry

   implicit none
   real(ESMF_KIND_R8), pointer   :: global_fieldptr(:), aux_global_fieldptr(:)
   type(ESMF_VM)                 :: vm1
   type(meshdata)                :: src_data, dst_data, global_src_data, global_dst_data
   type(hotdata)                 :: src_hotdata, dst_hotdata, global_dst_hotdata, global_src_hotdata
   type(regrid_data)             :: the_regrid_data
   type(ESMF_Mesh)               :: src_mesh, dst_mesh
   integer(ESMF_KIND_I4)         :: i1, rc, localPet, petCount
   character(len=6)              :: PE_ID
   character(len=*), parameter   :: src_fort14_dir = "coarse/", dst_fort14_dir = "fine/"
   real(ESMF_KIND_R8), parameter :: h0 = 0.05

#ifdef DEBUG_MODE
   if (localPet == 0) call show_message("This is adcirpolate, in debug mode.")
#endif

   !
   ! Any program using ESMF library should start with ESMF_Initialize(...).
   ! Next we get the ESMF_VM (virtual machine) and using this VM, we obtain
   ! the ID of our local PE and total number of PE's in the communicator.
   !
   call ESMF_Initialize(vm=vm1, defaultLogFilename="test.log", &
                        logKindFlag=ESMF_LOGKIND_MULTI, rc=rc)
   CHECK_ERR_CODE(localPet, rc)
   call ESMF_VMGet(vm=vm1, localPet=localPet, petCount=petCount, rc=rc)
   write (PE_ID, "(A,I4.4)") "PE", localPet

   !
   ! Next, we create our meshdata objects for source and destination meshes,
   ! and using those meshdata objects, we create the ESMF_Mesh objects, and
   ! we write the mesh into parallel vtu outputs.
   !
   call extract_parallel_data_from_mesh(vm1, src_fort14_dir, src_data)
   call write_meshdata_to_vtu(src_data, PE_ID//"_src_mesh.vtu", .true.)
   if (localPet == 0) call show_message("Creating parallel ESMF mesh from ADCIRC source mesh"//new_line("A"))
   call create_parallel_esmf_mesh_from_meshdata(src_data, src_mesh)

   call extract_parallel_data_from_mesh(vm1, dst_fort14_dir, dst_data)
   call write_meshdata_to_vtu(dst_data, PE_ID//"_dst_mesh.vtu", .true.)
   if (localPet == 0) call show_message("Creating parallel ESMF mesh from ADCIRC destination mesh"//new_line("A"))
   call create_parallel_esmf_mesh_from_meshdata(dst_data, dst_mesh)

   !
   ! Now, let us read data from fort.67. We also allocate the hotdata structure for
   ! destination mesh and fields.
   !
   if (localPet == 0) call show_message("Reading parallel hotdata from source directory."//new_line("A"))
   call extract_hotdata_from_parallel_binary_fort_67(src_data, src_hotdata, &
                                                     src_fort14_dir, .true.)
   call allocate_hotdata(dst_hotdata, dst_data)

   !
   ! After this point, we plan to overcome an important issue. The issue is
   ! if a point in the destination mesh is outside of the source mesh, we cannot
   ! use ESMF bilinear interpolation to transform our datafields. Hence, we first
   ! create a mask in the destination mesh, with its values equal to zero on the
   ! nodes outside of the source mesh domain, and one on the nodes which are inside
   ! the source mesh. Afterwards, we use bilinear interpolation for mask=1, and
   ! nearest node interpolation for mask=0. Thus, we need four ESMF_Fields:
   !   1- An ESMF_Field on the source mesh, which is used for the mask creation
   !      and also datafield interpolation.
   !   2- An ESMF_Field on the destination mesh, which will be only used for mask
   !      creation.
   !   3- An ESMF_Field on the destination mesh, which will be used for interpolating
   !      data on the destination points with mask=1.
   !   4- An ESMF_Field on the destination mesh, which will be used for interpolating
   !      data on the destination points with mask=0.
   !
   if (localPet == 0) call show_message("Creating ESMF fields:")
   the_regrid_data%src_datafield = ESMF_FieldCreate(mesh=src_mesh, &
                                                    typekind=ESMF_TYPEKIND_R8, rc=rc)
   CHECK_ERR_CODE(localPet, rc)
   if (localPet == 0) call show_message("source data field is created.")

   the_regrid_data%dst_mask_field = ESMF_FieldCreate(mesh=dst_mesh, &
                                                     typekind=ESMF_TYPEKIND_R8, rc=rc)
   if (rc == 0 .AND. localPet == 0) call show_message("destination mask field is created.")
   call check_error(__LINE__ - 2, __FILE__, rc)

   the_regrid_data%dst_mapped_field = ESMF_FieldCreate(mesh=dst_mesh, &
                                                       typekind=ESMF_TYPEKIND_R8, rc=rc)
   if (rc == 0 .AND. localPet == 0) call show_message("destination mapped data field is created.")
   call check_error(__LINE__ - 2, __FILE__, rc)

   the_regrid_data%dst_unmapped_field = ESMF_FieldCreate(mesh=dst_mesh, &
                                                         typekind=ESMF_TYPEKIND_R8, rc=rc)
   if (rc == 0 .AND. localPet == 0) call show_message("destination unmapped data field is created."//new_line("A"))
   call check_error(__LINE__ - 2, __FILE__, rc)

   !
   ! This is the preferred procedure in using ESMF to get a pointer to the
   ! ESMF_Field data array, and use that pointer for creating the mask, or
   ! assigning the data to field. ESMF_FieldGet also allocates the farrayptr
   !
   call ESMF_FieldGet(the_regrid_data%src_datafield, &
                      farrayPtr=the_regrid_data%src_fieldptr, rc=rc)
   call check_error(__LINE__ - 1, __FILE__, rc)

   call ESMF_FieldGet(the_regrid_data%dst_mask_field, &
                      farrayPtr=the_regrid_data%dst_maskptr, rc=rc)
   call check_error(__LINE__ - 1, __FILE__, rc)

   call ESMF_FieldGet(the_regrid_data%dst_mapped_field, &
                      farrayPtr=the_regrid_data%mapped_fieldptr, rc=rc)
   call check_error(__LINE__ - 1, __FILE__, rc)

   call ESMF_FieldGet(the_regrid_data%dst_unmapped_field, &
                      farrayPtr=the_regrid_data%unmapped_fieldptr, rc=rc)
   call check_error(__LINE__ - 1, __FILE__, rc)

   !
   ! At this section, we construct our interpolation operator (A matrix which maps
   ! source values to destination values). Here, no actual interpolation will happen,
   ! only the interpolation matrices will be constructed. We construct one matrix for
   ! nodal points with mask=1, and one for those points with mask=0.
   !
   if (localPet == 0) call show_message("Creating ESMF regriding operators:")
   call ESMF_FieldRegridStore(srcField=the_regrid_data%src_datafield, &
                              dstField=the_regrid_data%dst_mask_field, &
                              unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                              routeHandle=the_regrid_data%mapped_route_handle, &
                              regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=rc)
   if (localPet == 0 .AND. rc == 0) call show_message("mapped regriding operator is created.")
   call check_error(__LINE__, __FILE__, rc)
   call ESMF_FieldRegridStore(srcField=the_regrid_data%src_datafield, &
                              dstField=the_regrid_data%dst_unmapped_field, &
                              unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                              routeHandle=the_regrid_data%unmapped_route_handle, &
                              regridmethod=ESMF_REGRIDMETHOD_NEAREST_STOD, rc=rc)
   if (localPet == 0 .AND. rc == 0) call show_message("unmapped regriding operator is created."//new_line("A"))
   call check_error(__LINE__, __FILE__, rc)

   !
   ! This is the place that we create our mask on the destination mesh. By mask,
   ! we mean an array with length equal to number of nodes, whose values are equal
   ! to 1 at mapped nodes and 0 on unmapped nodes.
   !
   the_regrid_data%src_fieldptr = 1.d0
   call ESMF_FieldRegrid(srcField=the_regrid_data%src_datafield, &
                         dstField=the_regrid_data%dst_mask_field, &
                         routehandle=the_regrid_data%mapped_route_handle, rc=rc)

   !
   ! Now we map the nodal values of ETA1, ETA2, ETADisc, UU2, VV2, CH1
   ! from the source mesh to the destination mesh.
   !
   if (localPet == 0) call show_message("Regriding eta1:")
   call regrid_datafield_of_present_nodes(the_regrid_data, src_data, dst_data, src_hotdata%ETA1)
   dst_hotdata%ETA1 = the_regrid_data%mapped_fieldptr

   if (localPet == 0) call show_message("Regriding eta2:")
   call regrid_datafield_of_present_nodes(the_regrid_data, src_data, dst_data, src_hotdata%ETA2)
   dst_hotdata%ETA2 = the_regrid_data%mapped_fieldptr

   if (localPet == 0) call show_message("Regriding etaDisc:")
   call regrid_datafield_of_present_nodes(the_regrid_data, src_data, dst_data, src_hotdata%ETADisc)
   dst_hotdata%ETADisc = the_regrid_data%mapped_fieldptr

   if (localPet == 0) call show_message("Regriding UU2:")
   call regrid_datafield_of_present_nodes(the_regrid_data, src_data, dst_data, src_hotdata%UU2)
   dst_hotdata%UU2 = the_regrid_data%mapped_fieldptr

   if (localPet == 0) call show_message("Regriding VV2:")
   call regrid_datafield_of_present_nodes(the_regrid_data, src_data, dst_data, src_hotdata%VV2)
   dst_hotdata%VV2 = the_regrid_data%mapped_fieldptr

   if (localPet == 0) call show_message("Regriding CH1:")
   call regrid_datafield_of_present_nodes(the_regrid_data, src_data, dst_data, src_hotdata%CH1)
   dst_hotdata%CH1 = the_regrid_data%mapped_fieldptr

   if (localPet == 0) call show_message("Regriding NODECODE:")
   call regrid_datafield_of_present_nodes(the_regrid_data, src_data, dst_data, src_hotdata%realNODECODE)
   dst_hotdata%realNODECODE = the_regrid_data%mapped_fieldptr
   dst_hotdata%NNODECODE = nint(dst_hotdata%realNODECODE)

   if (localPet == 0) then
      call extract_global_data_from_fort14(src_data%dir_name//"/fort.14", global_src_data)
      call extract_global_data_from_fort14(dst_data%dir_name//"/fort.14", global_dst_data)
      call allocate_hotdata(global_src_hotdata, global_src_data)
      call allocate_hotdata(global_dst_hotdata, global_dst_data)
   end if

   call gather_src_nodal_hotdata_on_root(src_hotdata, global_src_hotdata, src_data, global_src_data, 0)
   call gather_dst_nodal_hotdata_on_root(dst_hotdata, global_dst_hotdata, dst_data, global_dst_data, 0)

   if (localPet == 0) then
      call show_message("Computing the wet and dry nodes and elements."//new_line("A"))
      call compute_wet_dry(global_dst_data, global_dst_hotdata, h0)
   end if

   !
   ! Finally, we want to visualize our results. This is not required in actual usage.
   ! We only do this for our presentation. So we write two meshes in the PE=0, and
   ! gather the interpolated field in PE=0. Then we plot these into vtu output.
   !
   if (localPet == 0) then
      !
      call write_meshdata_to_vtu(global_src_data, &
                                 src_data%dir_name//"/global_mesh.vtu", .false.)
      call write_int_node_field_to_vtu(global_src_hotdata%NNODECODE, &
                                       "NODECODE", src_data%dir_name//"/global_mesh.vtu", .false.)
      call write_node_field_to_vtu(global_src_hotdata%ETA1, &
                                   "ETA1", src_data%dir_name//"/global_mesh.vtu", .false.)
      call write_node_field_to_vtu(global_src_hotdata%ETA2, &
                                   "ETA2", src_data%dir_name//"/global_mesh.vtu", .false.)
      call write_node_field_to_vtu(global_src_hotdata%UU2, &
                                   "UU2", src_data%dir_name//"/global_mesh.vtu", .false.)
      call write_node_field_to_vtu(global_src_hotdata%VV2, &
                                   "VV2", src_data%dir_name//"/global_mesh.vtu", .true.)
      !
      call write_meshdata_to_vtu(global_dst_data, &
                                 dst_data%dir_name//"/global_mesh.vtu", .false., global_dst_hotdata)
      call write_int_node_field_to_vtu(global_dst_hotdata%NNODECODE, &
                                       "NODECODE", dst_data%dir_name//"/global_mesh.vtu", .false.)
      call write_node_field_to_vtu(global_dst_hotdata%ETA1, &
                                   "ETA1", dst_data%dir_name//"/global_mesh.vtu", .false.)
      call write_node_field_to_vtu(global_dst_hotdata%ETA2, &
                                   "ETA2", dst_data%dir_name//"/global_mesh.vtu", .false.)
      call write_node_field_to_vtu(global_dst_hotdata%UU2, &
                                   "UU2", dst_data%dir_name//"/global_mesh.vtu", .false.)
      call write_node_field_to_vtu(global_dst_hotdata%VV2, &
                                   "VV2", dst_data%dir_name//"/global_mesh.vtu", .true.)
   end if

   !
   ! Now we write the hotstart mesh for the destination mesh.
   !
   if (localPet == 0) then
      global_dst_hotdata%InputFileFmtVn = src_hotdata%InputFileFmtVn
      global_dst_hotdata%IMHS = src_hotdata%IMHS
      global_dst_hotdata%TimeLoc = src_hotdata%TimeLoc
      global_dst_hotdata%ITHS = src_hotdata%ITHS
      global_dst_hotdata%NP_G_IN = global_dst_data%NumNd
      global_dst_hotdata%NE_G_IN = global_dst_data%NumEl
      global_dst_hotdata%NP_A_IN = global_dst_data%NumNd
      global_dst_hotdata%NE_A_IN = global_dst_data%NumEl
      !
      ! By setting NNODECODE = 1, we consider that every node is a wet node,
      ! and we let ADCIRC correct this assumption on those nodes that do not
      ! satisfy the condition for a wet node. However, in Jan 2018, Jennifer
      ! found that ADCIRC will not correct this assumption in flood plains.
      ! So, we decided to interpolate NNODECODE along with other variables
      ! (such as ETA, UU, VV, ...).
      !
      ! global_dst_hotdata%NNODECODE = 1
      !
      global_dst_hotdata%NOFF = 1
      global_dst_hotdata%IESTP = 0
      global_dst_hotdata%NSCOUE = 0
      global_dst_hotdata%IVSTP = 0
      global_dst_hotdata%NSCOUV = 0
      global_dst_hotdata%ICSTP = 0
      global_dst_hotdata%NSCOUC = 0
      global_dst_hotdata%IPSTP = 0
      global_dst_hotdata%IWSTP = 0
      global_dst_hotdata%NSCOUM = 0
      global_dst_hotdata%IGEP = 0
      global_dst_hotdata%NSCOUGE = 0
      global_dst_hotdata%IGVP = 0
      global_dst_hotdata%NSCOUGV = 0
      global_dst_hotdata%IGCP = 0
      global_dst_hotdata%NSCOUGC = 0
      global_dst_hotdata%IGPP = 0
      global_dst_hotdata%IGWP = 0
      global_dst_hotdata%NSCOUGW = 0

      call show_message("Writing the global fort.67 in the destination mesh.")
      call write_serial_hotfile_to_fort_67(global_dst_data, global_dst_hotdata, &
                                           dst_fort14_dir, .false.)
   end if

   !
   ! Finally, we have to release the memory.
   !
   if (localPet == 0) then
      call destroy_meshdata(global_dst_data)
      call destroy_meshdata(global_src_data)
   end if
   call destroy_regrid_data(the_regrid_data)
   call ESMF_MeshDestroy(dst_mesh)
   call ESMF_MeshDestroy(src_mesh)
   call destroy_meshdata(src_data)
   call destroy_meshdata(dst_data)

   call ESMF_Finalize()

end program main
