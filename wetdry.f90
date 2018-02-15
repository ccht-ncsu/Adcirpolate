!> @author Ali Samii - 2018
!! Ali Samii - Department of ASE/EM, UT Austin
!! @brief This module is an interface between parallel Adcirc input files and ESMF library.
module wetdry
   use adcirpolate
   use ESMF

   integer, parameter                     :: NumND_per_El = 3
   integer(ESMF_KIND_I4)                  :: nc_change
   integer(ESMF_KIND_I4), allocatable     :: nibcnt(:), noffold(:)

contains

   !>
   !!
   !!
   subroutine update_wet_dry_nodes_D1(in_mesh_data, in_hot_data, h0)
      implicit none
      type(meshdata), intent(in)               :: in_mesh_data
      type(hotdata), intent(inout)             :: in_hot_data
      real(ESMF_KIND_R8), intent(in)           :: h0
      real(ESMF_KIND_R8)                       :: htot
      integer(ESMF_KIND_I4)                    :: i1, j1
      !
      ! if htot < h0, the node is dry
      !
      do i1 = 1, in_mesh_data%NumNd
         if (in_hot_data%NNODECODE(i1) == 1) then
            htot = in_mesh_data%bathymetry(i1) + in_hot_data%ETA2(i1)
            if (htot <= h0) then
               in_hot_data%NNODECODE(i1) = 0
            end if
         end if
      end do
   end subroutine

   !>
   !!
   !!
   subroutine update_wet_dry_nodes_W12(in_mesh_data, in_hot_data, h0)
      implicit none
      type(meshdata), intent(in)               :: in_mesh_data
      type(hotdata), intent(inout)             :: in_hot_data
      real(ESMF_KIND_R8),intent(in)            :: h0
      real(ESMF_KIND_R8)                       :: htotN1, htotN2, htotN3, ETAN1, ETAN2, ETAN3, HOFF
      integer(ESMF_KIND_I4)                    :: I, nm1, nm2, nm3, nm123, nctot
      !
      HOFF = 1.2d0*h0
      !
      do I = 1, in_mesh_data%NumEl
         nm1 = in_mesh_data%ElConnect((I - 1)*NumND_per_El + 1)
         nm2 = in_mesh_data%ElConnect((I - 1)*NumND_per_El + 2)
         nm3 = in_mesh_data%ElConnect((I - 1)*NumND_per_El + 3)
         ! WET...
         ! WET...Nodal Wetting Criteria W1: This depends on changes that occurred in D1
         ! WET...
         nctot = in_hot_data%NNODECODE(nm1) + in_hot_data%NNODECODE(nm2) + in_hot_data%NNODECODE(nm3)
         if (nctot .EQ. 2) then
            ETAN1 = in_hot_data%ETA2(nm1)
            ETAN2 = in_hot_data%ETA2(nm2)
            ETAN3 = in_hot_data%ETA2(nm3)
            HTOTN1 = in_mesh_data%bathymetry(nm1) + ETAN1
            HTOTN2 = in_mesh_data%bathymetry(nm2) + ETAN2
            HTOTN3 = in_mesh_data%bathymetry(nm3) + ETAN3
            if ((in_hot_data%NNODECODE(nm1) .EQ. 1) .AND. (in_hot_data%NNODECODE(nm2) .EQ. 1)) then
               if ((HTOTN1 .GE. HOFF) .AND. (HTOTN2 .GE. HOFF)) then
                  nm123 = nm1
                  if (in_hot_data%ETA2(nm2) .GT. in_hot_data%ETA2(nm1)) nm123 = nm2
               end if
            else if ((in_hot_data%NNODECODE(nm2) .EQ. 1) .AND. (in_hot_data%NNODECODE(nm3) .EQ. 1)) then
            else if ((in_hot_data%NNODECODE(nm3) .EQ. 1) .AND. (in_hot_data%NNODECODE(nm1) .EQ. 1)) then
            end if
         end if
         ! WET...
         ! WET...Nodal Wetting Criteria W2a
         ! WET...
!         NBnctot = NIBNODECODE(nm1) + NIBNODECODE(nm2) + NIBNODECODE(nm3)
!         NIBCNT(nm1) = NIBCNT(nm1) + NBnctot
!         NIBCNT(nm2) = NIBCNT(nm2) + NBnctot
!         NIBCNT(nm3) = NIBCNT(nm3) + NBnctot

      end do
   end subroutine

   !>
   !!
   !!
   subroutine init_wet_dry_elements(in_mesh_data, in_hot_data, h0)
      implicit none
      type(meshdata), intent(in)               :: in_mesh_data
      type(hotdata), intent(inout)             :: in_hot_data
      real(ESMF_KIND_R8), intent(in)           :: h0
      integer(ESMF_KIND_I4)                    :: i1, j1, el_connect(NumND_per_El), el_node_code(NumND_per_El), sum_node_code
      !
      allocate (nibcnt(in_mesh_data%NumNd))
      allocate (noffold(in_mesh_data%NumEl))
      !
      ! We set an element dry, if two of its nodes are dry.
      !
      nibcnt(:) = 0
      in_hot_data%NOFF(:) = 1
      do i1 = 1, in_mesh_data%NumEl
         do j1 = 1, NumND_per_El
            el_connect(j1) = in_mesh_data%ElConnect((i1 - 1)*NumND_per_El + j1)
            el_node_code(j1) = in_hot_data%NNODECODE(el_connect(j1))
         end do
         sum_node_code = sum(el_node_code)
         if (sum_node_code < 3) then
            in_hot_data%NOFF(i1) = 0
         end if
      end do
      !
      noffold(:) = in_hot_data%NOFF(:)
   end subroutine

   !>
   !!
   !!
   subroutine compute_wet_dry(in_mesh_data, in_hot_data, h0)
      implicit none
      type(meshdata), intent(in)               :: in_mesh_data
      type(hotdata), intent(inout)             :: in_hot_data
      real(ESMF_KIND_R8), intent(in)           :: h0
      real(ESMF_KIND_R8)                       :: habsmin, hoff
      integer(ESMF_KIND_I4)                    :: i1, j1, el_connect(NumND_per_El), el_node_code(NumND_per_El)
      !
      habsmin = 0.8d0*h0
      hoff = 1.2d0*h0
      !
      call update_wet_dry_nodes_D1(in_mesh_data, in_hot_data, h0)
      !
      ! We are not calling W12 wetting and drying.
      !
      !call update_wet_dry_nodes_W12(in_mesh_data, in_hot_data, h0)
      call init_wet_dry_elements(in_mesh_data, in_hot_data, h0)
   end subroutine

end module wetdry

