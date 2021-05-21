      MODULE splitter

      use ESMF

      implicit none

      CONTAINS

      SUBROUTINE decompose_global_mesh(vm, global_data, the_data)

      use adcirpolate, only: meshdata

      implicit none

      type(ESMF_VM),  intent(in)    :: vm
      type(meshdata), intent(in)    :: global_data
      type(meshdata), intent(inout) :: the_data

      integer(ESMF_KIND_I4) :: localPet
      integer(ESMF_KIND_I4) :: petCount

      integer(ESMF_KIND_I4) :: counter
      integer(ESMF_KIND_I4) :: ele_finish
      integer(ESMF_KIND_I4) :: ele_start
      integer(ESMF_KIND_I4) :: first_ele
      integer(ESMF_KIND_I4) :: ge
      integer(ESMF_KIND_I4) :: gv
      integer(ESMF_KIND_I4) :: gvv
      integer(ESMF_KIND_I4) :: ie
      integer(ESMF_KIND_I4) :: iv
      integer(ESMF_KIND_I4) :: ivv
      integer(ESMF_KIND_I4) :: lv
      integer(ESMF_KIND_I4) :: owner

      integer(ESMF_KIND_I4), allocatable :: partmesh(:)
      integer(ESMF_KIND_I4), allocatable :: vert_found(:)
      integer(ESMF_KIND_I4), allocatable :: vert_g2l(:)

      integer(ESMF_KIND_I4), parameter :: dim1 = 2
      integer(ESMF_KIND_I4), parameter :: NumND_per_El = 3

      the_data%vm = vm
      call ESMF_VMGet(vm=vm, localPet=localPet, petCount=petCount)

      the_data%dir_name = global_data%dir_name

      ! Determine the number of elements on this partition.
      ! Instead of doing a fancy decomposition (with METIS, etc.),
      ! let's assume each core gets NE/NP elements.
      the_data%NumEl = floor(real(global_data%NumEl)/real(petCount))
      ele_start  =  localPet    * the_data%NumEl + 1
      ele_finish = (localPet+1) * the_data%NumEl
      if( (localPet+1).eq.petCount )then
        ele_finish = global_data%NumEl
        the_data%NumEl = ele_finish - ele_start + 1
      endif

      ! Populate the list of global elements present on this partition.
      allocate (the_data%ElIDs(the_data%NumEl))
      do ie=1,the_data%NumEl
        the_data%ElIDs(ie) = (ele_start-1) + ie
      enddo

      ! Determine which global vertices are present on this partition.
      allocate (vert_found(global_data%NumNd))
      vert_found = 0
      do ie=ele_start,ele_finish
        do iv=1,NumNd_per_El
          vert_found(global_data%ElConnect((ie-1)*NumND_per_El+iv)) = 1
        enddo
      enddo
      the_data%NumNd = SUM(vert_found)
      allocate (the_data%NdIDs(the_data%NumNd))
      allocate (vert_g2l(global_data%NumNd))
      counter = 0
      do iv=1,global_data%NumNd
        if(vert_found(iv))then
          counter = counter + 1
          the_data%NdIDs(counter) = iv
          vert_g2l(iv) = counter
        endif
      enddo

      ! Populate the coordinates of the vertices present on this partition.
      allocate (the_data%NdCoords(dim1*the_data%NumNd))
      do iv=1,the_data%NumNd
        the_data%NdCoords((iv-1)*dim1+1)                               &
     &          = global_data%NdCoords((the_data%NdIDs(iv)-1)*dim1+1)
        the_data%NdCoords((iv-1)*dim1+2)                               &
     &          = global_data%NdCoords((the_data%NdIDs(iv)-1)*dim1+2)
      enddo

      ! Populate the element-to-vertex connectivity.
      ! This list must contain the local vertex numbers for each
      ! local element, so we have to develop this from global information.
      allocate (the_data%ElTypes(the_data%NumEl))
      the_data%ElTypes = ESMF_MESHELEMTYPE_TRI
      allocate (the_data%ElConnect(NumNd_per_El*the_data%NumEl))
      do ie=1,the_data%NumEl
        ge = the_data%ElIDs(ie)
        do ivv=1,NumNd_per_El
          gv = global_data%ElConnect((ge-1)*NumNd_per_El+ivv)
          lv = vert_g2l(gv)
          the_data%ElConnect((ie-1)*NumNd_per_El+ivv) = lv
        enddo
      enddo
#ifdef SLOW 
      allocate (the_data%ElConnect(NumND_per_el*the_data%NumEl))
      do ie=1,the_data%NumEl
        ge    = the_data%ElIDs(ie)
        gv(1) = global_data%ElConnect((ge-1)*NumNd_per_El+1)
        gv(2) = global_data%ElConnect((ge-1)*NumNd_per_El+2)
        gv(3) = global_data%ElConnect((ge-1)*NumNd_per_El+3)
        do iv=1,the_data%NumNd
          if( the_data%NdIDs(iv).eq.gv(1) )then
            the_data%ElConnect((ie-1)*NumND_per_el+1) = iv
          endif
          if( the_data%NdIDs(iv).eq.gv(2) )then
            the_data%ElConnect((ie-1)*NumND_per_el+2) = iv
          endif
          if( the_data%NdIDs(iv).eq.gv(3) )then
            the_data%ElConnect((ie-1)*NumND_per_el+3) = iv
          endif
        enddo
      enddo
#endif

      ! Populate the list of vertex owners.
      ! This is tricky because we don't have information from a full
      ! decomposition (e.g. fort.18 files). All of the elements in our
      ! ElIDs list will be owned by this partition, but many of the
      ! vertices in our NdIDs list will be shared among partitions.
      ! Assume that a vertex is owned by the partition that also owns
      ! the first element attached to that vertex. For example, if
      ! vertex #10 is shared by elements #4, #8, #15, and #22, then it
      ! will be owned by the partition that also owns element #4.
      allocate (partmesh(global_data%NumNd))
      partmesh = petCount - 1
      do ie=1,global_data%NumEl
        owner = floor( real(ie-1)                                      &
     &        / real(floor(real(global_data%NumEl)/real(petCount))) )
        do ivv=1,NumNd_per_El
          gvv = global_data%ElConnect((ie-1)*NumNd_per_El+ivv)
          if( partmesh(gvv).gt.owner )then
            partmesh(gvv) = owner
          endif
        enddo
      enddo
      allocate (the_data%NdOwners(the_data%NumNd))
      do iv=1,the_data%NumNd
        the_data%NdOwners(iv) = partmesh(the_data%NdIDs(iv))
      enddo
#ifdef SLOW
      allocate (the_data%NdOwners(the_data%NumNd))
      the_data%NdOwners = -1
      do iv=1,the_data%NumNd
        first_ele = global_data%NumEl
        do ie=1,global_data%NumEl
          do ivv=1,NumNd_per_El
            gvv = global_data%ElConnect((ie-1)*NumNd_per_El+ivv)
            if( gvv.eq.the_data%NdIDs(iv) .and. ie.lt.first_ele )then
              first_ele = ie
            endif
          enddo
        enddo
        the_data%NdOwners(iv) = floor( real(first_ele-1)               &
       &      / real(floor(real(global_data%NumEl)/real(petCount))) )
        if( the_data%NdOwners(iv).gt.(petCount-1) )then
          the_data%NdOwners(iv) = petCount - 1
        endif
      enddo
#endif

      ! Using the vertex owners, set up the information needed
      ! to collect information onto the root processor, so it can write
      ! the global hot-start file at the end of the process.
      the_data%NumOwnedNd = 0
      do iv=1,the_data%NumNd
        if( the_data%NdOwners(iv).eq.localPet )then
          the_data%NumOwnedNd = the_data%NumOwnedNd + 1
        endif
      enddo
      allocate (the_data%owned_to_present_nodes(the_data%NumOwnedNd))
      counter = 0
      do iv=1,the_data%NumNd
        if( the_data%NdOwners(iv).eq.localPet )then
          counter = counter + 1
          the_data%owned_to_present_nodes(counter) = iv
        endif
      enddo

      ! And we're finished?
      the_data%is_initialized = .true.

      ENDSUBROUTINE

      ENDMODULE
