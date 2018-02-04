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
!                  if (ETA2(nm2) .GT. ETA2(nm1)) nm123 = nm2
!                  DELDIST = SQRT((y(nm3) - y(nm123))**2.D0 &
!                                 + (X(nm3) - X(nm123))**2.D0)
!                  DELETA = ETA2(nm123) - ETA2(nm3)
!                  ! jgf50.60.18: Prevent numerical problems if DELETA is negative
!                  if (DELETA .lt. 0.d0) DELETA = 0.d0
!                  H1 = ETA2(nm123) + DP(nm123)
!                   !  RJW merged from Casey 071219: Added the following logic for 3D friction.
!                   !  RJW modified the following for 3D friction
!                  if (C2DDI) then
!                   ! sb46.28sb02
!                   !<<                     Convert Manning's N to Cd, if necessary.
!                     if (LoadManningsN) then
!                        FRIC(nm123) = g*ManningsN(nm123)**2.d0 &
!                                      /((DP(nm123) + ifNLFA*ETA2(nm123)) &
!                                        **(1.d0/3.d0))
!                        if (FRIC(nm123) .LT. BFCdLLimit) then
!                           FRIC(nm123) = BFCdLLimit
!                        end if
!                     end if
!                   !>>
!                     TKWET = FRIC(nm123)*(ifLINBF + (VELMIN/H1)* &
!                                          (ifNLBF + ifHYBF* &
!                                           (1.D0 + (HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
!
!                     if (TKWET .LT. 0.0001d0) TKWET = 0.0001d0
!                     VEL = G*(DELETA/DELDIST)/TKWET
!
!                  else if (C3D) then
                    !!C solve for the depth averaged velocity,U, from the relation :
                    !!C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
                    !!C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
                    !!C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
!                     if (LoadZ0B_var) then
!                        Z0B1 = Z0B_var(nm123)
!                     else if (LoadManningsN) then
!                        Z0B1 = (DP(nm123) + ifNLFA*ETA2(nm123))*exp(-(1.0D0 + &
!                                                                      ((0.41D0*(DP(nm123) + ifNLFA*ETA2(nm123))**(1.0D0/6.0D0))/ &
!                                                                       (ManningsN(nm123)*sqrt(g)))))
!                     else
!                        Z0B1 = Z0B
!                     end if
!                     VEL = sqrt(g*H1*(DELETA/DELDIST)) &
!                           *((H1 + Z0B1)*LOG((H1 + Z0B1)/Z0B1) - H1)/(H1*0.41D0)
!                  end if
!
!                  if (VEL .GT. VELMIN) then
!                   !    ....         third node met criteria and is also wet
!                     NNODECODE(nm3) = 1
!                   !  RJW merged 08/26/20008 Casey 071219: Added the following logic to obtain the correct friction.
!                     if (C2DDI) then
!                   !                           TK(nm123)=FRIC(nm123)*(ifLINBF+(VEL/H1)*
!                        TK(nm123) = FRIC(nm123)*(ifLINBF + (VEL/H1)* &
!                                                 (ifNLBF + ifHYBF* &
!                                                  (1.D0 + (HBREAK/H1)**FTHETA)** &
!                                                  (FGAMMA/FTHETA)))
!                     else if (C3D) then
!                        if (ISLIP .EQ. 0) then
!                           DUDS = (Q(nm123, 2) - Q(nm123, 1)) &
!                                  /(SIGMA(2) - SIGMA(1))
!                           BSX1(nm123) = EVTOT(1)*REAL(DUDS)
!                           BSY1(nm123) = EVTOT(1)*AIMAG(DUDS)
!                           BSX1(nm3) = EVTOT(1)*REAL(DUDS)
!                           BSY1(nm3) = EVTOT(1)*AIMAG(DUDS)
!                        end if
!                        if (ISLIP .NE. 0) then
!                           if (ISLIP .EQ. 1) then
!                              KSLIP = KP
!                           end if
!                           if (ISLIP .EQ. 2) then
!                              KSLIP = (1.D0/ &
!                                       ((1.D0/0.41D0)* &
!                                        LOG((ABS(((SIGMA(2) - SIGMA(1))/(A - B))* &
!                                                 (DP(nm123) + ifNLFA*ETA2(nm123))) &
!                                             + Z0B1) &
!                                            /(Z0B1))))**2.D0 &
!                                      *ABS(Q(nm123, 1))
!                              if (KSLIP .GT. 1.D0*ABS(Q(nm123, 1))) &
!                                 KSLIP = 1.D0*ABS(Q(nm123, 1))
!                              if (KSLIP .LT. 0.0025D0*ABS(Q(nm123, 1))) &
!                                 KSLIP = 0.0025D0*ABS(Q(nm123, 1))
!                           end if
!                           BSX1(nm123) = KSLIP*REAL(Q(nm123, 1))
!                           BSY1(nm123) = KSLIP*AIMAG(Q(nm123, 1))
!                           BSX1(nm3) = KSLIP*REAL(Q(nm123, 1))
!                           BSY1(nm3) = KSLIP*AIMAG(Q(nm123, 1))
!                        end if
!                     end if
!                  end if
               end if
            else if ((in_hot_data%NNODECODE(nm2) .EQ. 1) .AND. (in_hot_data%NNODECODE(nm3) .EQ. 1)) then
!               if ((HTOTN2 .GE. HOFF) .AND. (HTOTN3 .GE. HOFF)) then
!                  nm123 = nm2
!                  if (ETA2(nm3) .GT. ETA2(nm2)) nm123 = nm3
!                  DELDIST = SQRT((Y(nm1) - Y(nm123))**2.D0 &
!                                 + (X(nm1) - X(nm123))**2.D0)
!                  DELETA = ETA2(nm123) - ETA2(nm1)
!                  ! jgf50.60.18: Prevent numerical problems if DELETA is negative
!                  if (DELETA .lt. 0.d0) DELETA = 0.d0
!                  H1 = ETA2(nm123) + DP(nm123)
!                       !  RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
!                  if (C2DDI) then
!                       ! sb46.28sb02
!                       !C<<                     Convert Manning's N to Cd, if necessary.
!                     if (LoadManningsN) then
!                        FRIC(nm123) = g*ManningsN(nm123)**2.d0 &
!                                      /((DP(nm123) + ifNLFA*ETA2(nm123)) &
!                                        **(1.d0/3.d0))
!                        if (FRIC(nm123) .LT. BFCdLLimit) then
!                           FRIC(nm123) = BFCdLLimit
!                        end if
!                     end if
!                       !C>>
!                     TKWET = FRIC(nm123)*(ifLINBF + (VELMIN/H1)* &
!                                          (ifNLBF + ifHYBF* &
!                                           (1.D0 + (HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
!                     if (TKWET .LT. 0.0001d0) TKWET = 0.0001d0
!                     VEL = G*(DELETA/DELDIST)/TKWET
!
!                  else if (C3D) then
                    !! solve for the depth averaged velocity,U, from the relation :
                    !!        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
                    !!          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
                    !! where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
!                     if (LoadZ0B_var) then
!                        Z0B1 = Z0B_var(nm123)
!                     else if (LoadManningsN) then
!                        Z0B1 = (DP(nm123) + ifNLFA*ETA2(nm123))*exp(-(1.0D0 + &
!                                                                      ((0.41D0*(DP(nm123) + ifNLFA*ETA2(nm123))**(1.0D0/6.0D0))/ &
!                                                                       (ManningsN(nm123)*sqrt(g)))))
!                     else
!                        Z0B1 = Z0B
!                     end if
!                     VEL = sqrt(g*H1*(DELETA/DELDIST)) &
!                           *((H1 + Z0B1)*LOG((H1 + Z0B1)/Z0B1) - H1)/(H1*0.41D0)
!                  end if
!
!                  if (VEL .GT. VELMIN) then
!                     NNODECODE(nm1) = 1
!                   !  RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
!                     if (C2DDI) then
!                        TK(nm123) = FRIC(nm123)*(ifLINBF + (VEL/H1)* &
!                                                 (ifNLBF + ifHYBF* &
!                                                  (1.D0 + (HBREAK/H1)**FTHETA)** &
!                                                  (FGAMMA/FTHETA)))
!                     else if (C3D) then
!                        if (ISLIP .EQ. 0) then
!                           DUDS = (Q(nm123, 2) - Q(nm123, 1)) &
!                                  /(SIGMA(2) - SIGMA(1))
!                           BSX1(nm123) = EVTOT(1)*REAL(DUDS)
!                           BSY1(nm123) = EVTOT(1)*AIMAG(DUDS)
!                           BSX1(nm1) = EVTOT(1)*REAL(DUDS)
!                           BSY1(nm1) = EVTOT(1)*AIMAG(DUDS)
!                        end if
!                        if (ISLIP .NE. 0) then
!                           if (ISLIP .EQ. 1) then
!                              KSLIP = KP
!                           end if
!                           if (ISLIP .EQ. 2) then
!                              KSLIP = (1.D0/ &
!                                       ((1.D0/0.41D0)* &
!                                        LOG((ABS(((SIGMA(2) - SIGMA(1))/(A - B))* &
!                                                 (DP(nm123) + ifNLFA*ETA2(nm123))) &
!                                             + Z0B1) &
!                                            /(Z0B1))))**2.D0 &
!                                      *ABS(Q(nm123, 1))
!                              if (KSLIP .GT. 1.D0*ABS(Q(nm123, 1))) &
!                                 KSLIP = 1.D0*ABS(Q(nm123, 1))
!                              if (KSLIP .LT. 0.0025D0*ABS(Q(nm123, 1))) &
!                                 KSLIP = 0.0025D0*ABS(Q(nm123, 1))
!                           end if
!                           BSX1(nm123) = KSLIP*REAL(Q(nm123, 1))
!                           BSY1(nm123) = KSLIP*AIMAG(Q(nm123, 1))
!                           BSX1(nm1) = KSLIP*REAL(Q(nm123, 1))
!                           BSY1(nm1) = KSLIP*AIMAG(Q(nm123, 1))
!                        end if
!                     end if
!
!                  end if
!               end if
!               ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            else if ((in_hot_data%NNODECODE(nm3) .EQ. 1) .AND. (in_hot_data%NNODECODE(nm1) .EQ. 1)) then
!               if ((HTOTN3 .GE. HOFF) .AND. (HTOTN1 .GE. HOFF)) then
!                  nm123 = nm3
!                  if (ETA2(nm1) .GT. ETA2(nm3)) nm123 = nm1
!                  DELDIST = SQRT((Y(nm2) - Y(nm123))**2.D0 &
!                                 + (X(nm2) - X(nm123))**2.D0)
!                  DELETA = ETA2(nm123) - ETA2(nm2)
!                  ! jgf50.60.18: Prevent numerical problems if DELETA is negative
!                  if (DELETA .lt. 0.d0) DELETA = 0.d0
!                  H1 = ETA2(nm123) + DP(nm123)
!                   !  RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
!                  if (C2DDI) then
!                   ! sb46.28sb02
!                   !C<<                     Convert Manning's N to Cd, if necessary.
!                     if (LoadManningsN) then
!                        FRIC(nm123) = g*ManningsN(nm123)**2.d0 &
!                                      /((DP(nm123) + ifNLFA*ETA2(nm123)) &
!                                        **(1.d0/3.d0))
!                        if (FRIC(nm123) .LT. BFCdLLimit) then
!                           FRIC(nm123) = BFCdLLimit
!                        end if
!                     end if
!                   !C>>
!                     TKWET = FRIC(nm123)*(ifLINBF + (VELMIN/H1)* &
!                                          (ifNLBF + ifHYBF* &
!                                           (1.D0 + (HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
!                     if (TKWET .LT. 0.0001d0) TKWET = 0.0001d0
!                     VEL = G*(DELETA/DELDIST)/TKWET
!                  else if (C3D) then
!                   ! solve for the depth averaged velocity,U, from the relation :
!                   !        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
!                   !          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
!                   ! where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
!                     if (LoadZ0B_var) then
!                        Z0B1 = Z0B_var(nm123)
!                     else if (LoadManningsN) then
!                        Z0B1 = (DP(nm123) + ifNLFA*ETA2(nm123))*exp(-(1.0D0 + &
!                                                                      ((0.41D0*(DP(nm123) + ifNLFA*ETA2(nm123))**(1.0D0/6.0D0))/ &
!                                                                       (ManningsN(nm123)*sqrt(g)))))
!                     else
!                        Z0B1 = Z0B
!                     end if
!                     VEL = sqrt(g*H1*(DELETA/DELDIST)) &
!                           *((H1 + Z0B1)*LOG((H1 + Z0B1)/Z0B1) - H1)/(H1*0.41D0)
!                  end if
!
!                  if (VEL .GT. VELMIN) then
!                     NNODECODE(nm2) = 1
!                     ! RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
!                     if (C2DDI) then
!                        TK(nm123) = FRIC(nm123)*(ifLINBF + (VEL/H1)* &
!                                                 (ifNLBF + ifHYBF* &
!                                                  (1.D0 + (HBREAK/H1)**FTHETA)** &
!                                                  (FGAMMA/FTHETA)))
!                     else if (C3D) then
!                        if (ISLIP .EQ. 0) then
!                           DUDS = (Q(nm123, 2) - Q(nm123, 1)) &
!                                  /(SIGMA(2) - SIGMA(1))
!                           BSX1(nm123) = EVTOT(1)*REAL(DUDS)
!                           BSY1(nm123) = EVTOT(1)*AIMAG(DUDS)
!                           BSX1(nm2) = EVTOT(1)*REAL(DUDS)
!                           BSY1(nm2) = EVTOT(1)*AIMAG(DUDS)
!                        end if
!                        if (ISLIP .NE. 0) then
!                           if (ISLIP .EQ. 1) then
!                              KSLIP = KP
!                           end if
!                           if (ISLIP .EQ. 2) then
!                              KSLIP = (1.D0/ &
!                                       ((1.D0/0.41D0)* &
!                                        LOG((ABS(((SIGMA(2) - SIGMA(1))/(A - B))* &
!                                                 (DP(nm123) + ifNLFA*ETA2(nm123))) &
!                                             + Z0B1) &
!                                            /(Z0B1))))**2.D0 &
!                                      *ABS(Q(nm123, 1))
!                              if (KSLIP .GT. 1.D0*ABS(Q(nm123, 1))) &
!                                 KSLIP = 1.D0*ABS(Q(nm123, 1))
!                              if (KSLIP .LT. 0.0025D0*ABS(Q(nm123, 1))) &
!                                 KSLIP = 0.0025D0*ABS(Q(nm123, 1))
!                           end if
!                           BSX1(nm123) = KSLIP*REAL(Q(nm123, 1))
!                           BSY1(nm123) = KSLIP*AIMAG(Q(nm123, 1))
!                           BSX1(nm2) = KSLIP*REAL(Q(nm123, 1))
!                           BSY1(nm2) = KSLIP*AIMAG(Q(nm123, 1))
!                        end if
!                     end if
!
!                  end if
!               end if
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
      integer(ESMF_KIND_I4)                    :: i1, j1, el_connect(NumND_per_El), el_node_code(NumND_per_El)
      !
      allocate (nibcnt(in_mesh_data%NumNd))
      allocate (noffold(in_mesh_data%NumEl))
      !
      ! We set an element wet, if all of its nodes are wet
      !
      nibcnt(:) = 0
      in_hot_data%NOFF(:) = 1
      do i1 = 1, in_mesh_data%NumEl
         do j1 = 1, NumND_per_El
            el_connect(j1) = in_mesh_data%ElConnect((i1 - 1)*NumND_per_El + j1)
            el_node_code(j1) = in_hot_data%NNODECODE(el_connect(j1))
         end do
         if (el_node_code(1) /= 1 .OR. el_node_code(2) /= 1 .OR. el_node_code(3) /= 1) then
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
      call update_wet_dry_nodes_W12(in_mesh_data, in_hot_data, h0)
      call init_wet_dry_elements(in_mesh_data, in_hot_data, h0)
   end subroutine

end module wetdry

