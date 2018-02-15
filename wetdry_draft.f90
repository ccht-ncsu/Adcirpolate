      module wetdry

         real(sz) :: habsmin
         real(sz) :: hoff
         real(sz) :: h0 ! read h0 from fort.15

         integer :: mnp ! number of nodes
         integer :: mne ! number of cells

         integer, allocatable ::    nibcnt(:)
         integer, allocatable ::    noffold(:)

         subroutine initializeWettingAndDrying(h0)
            use sizes, only: mne, mnp
            use global, only: h0, nodecode, nnodecode
            implicit none

            allocate (nibcnt(mnp))
            allocate (noffold(mne))
            nnodecode = 1
            nodecode = 1
            noffold(:) = 1
            habsmin = 0.8d0*h0
            hoff = 1.2d0*h0

!----------------------------------------------------------------------
         end subroutine initializeWettingAndDrying
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!                   S U B R O U T I N E
!     C O M P U T E   W E T T I N G   A N D   D R Y I N G
!----------------------------------------------------------------------
!     Determines which nodes should be wet and which should be dry.
!----------------------------------------------------------------------
         subroutine computeWettingAndDrying(it)
            use sizes, only: sz, mne, mnp
            use global, only: noff, nodecode, nnodecode, eta2, tk, nolifa, &
               bsx1, bsy1, btime_end, C2DDI, C3D, g, h0, ifnlfa, nddt, &
               nibnodecode, ilump, ncchange &


            use mesh, only: ne, np, dp, mju, totalArea, nm, x, y, areas
            use nodalattributes, only: BFCdLLimit, fgamma, ftheta, fric, &
               manningsn, hbreak, ifhybf, ifnlbf, iflinbf, loadManningsN, &
               loadZ0B_var, z0b_var
            use global_3dvs, only: a, b, islip, kp, z0b, sigma, evtot, q
            use subdomain, only: subdomainOn, enforceBN, enforceWDcb, enforceWDob
            implicit none

            integer, intent(in) :: it ! time step number
            complex(sz) :: duds !jgf48.50 declare size SZ instead of plain COMPLEX
            integer :: nc1, nc2, nc3
            integer :: nm1, nm2, nm3
            integer :: nm123
            integer :: ncele
            integer :: nctot
            real(sz) :: kslip
            real(sz) :: vel
            real(sz) :: z0b1
            real(sz) :: areaEle
            real(sz) :: tkWet
            real(sz) :: etaN1, etaN2, etaN3
            real(sz) :: hTotN1, hTotN2, hTotN3
            real(sz) :: deldist, deleta
            real(sz) :: htot
            real(sz) :: h1
            integer :: nbnctot
            integer :: i
            integer :: ie

! WET...
! WET...WET/DRY - INITIALIZATIONS FOR WET/DRY LOOP
! WET...
            DO I = 1, NP
               NIBCNT(I) = 0
            ENDDO
            DO I = 1, NE
               NOFFOLD(I) = NOFF(I)
               NOFF(I) = 1
            ENDDO
! WET...
! WET...WET/DRY - PART 1 - NODAL DRYING CRITERIA D1
! WET....Drying Criteria D1: this depends on NODECODE and updates NODECODE
! WET...
            DO I = 1, NP
               IF (NODECODE(I) .EQ. 1) THEN
                  HTOT = DP(I) + ETA2(I)
                  IF (HTOT .LE. H0) THEN
                     IF (HTOT .LT. HABSMIN) ETA2(I) = HABSMIN - DP(I)
                     NNODECODE(I) = 0
                     NODECODE(I) = 0
                     NCCHANGE = NCCHANGE + 1 !NCCHANGE=0 set near beginning of GWCE
!                  ENDIF
                  ENDIF
               ENDIF
            ENDDO
! WET...
! WET...END WET/DRY SECTION - PART 1
! WET...

! WET...
! WET...WET/DRY SECTION PART 2 - NODAL WETTING LOOPS W1 AND W2
! WET...
            DO I = 1, NE
               NM1 = NM(I, 1)
               NM2 = NM(I, 2)
               NM3 = NM(I, 3)

! WET...
! WET...Nodal Wetting Criteria W1: This depends on changes that occurred in D1
! WET...
               NCTOT = NODECODE(NM1) + NODECODE(NM2) + NODECODE(NM3)
               IF (NCTOT .EQ. 2) THEN
                  ETAN1 = ETA2(NM1)
                  ETAN2 = ETA2(NM2)
                  ETAN3 = ETA2(NM3)
                  HTOTN1 = DP(NM1) + ETA2(NM1)
                  HTOTN2 = DP(NM2) + ETA2(NM2)
                  HTOTN3 = DP(NM3) + ETA2(NM3)
                  IF ((NODECODE(NM1) .EQ. 1) .AND. (NODECODE(NM2) .EQ. 1)) THEN
                     IF ((HTOTN1 .GE. HOFF) .AND. (HTOTN2 .GE. HOFF)) THEN
                        NM123 = NM1
                        IF (ETA2(NM2) .GT. ETA2(NM1)) NM123 = NM2
                        DELDIST = SQRT((y(NM3) - y(NM123))**2.D0 &
                                       + (X(NM3) - X(NM123))**2.D0)
                        DELETA = ETA2(NM123) - ETA2(NM3)
                        ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                        IF (DELETA .lt. 0.d0) DELETA = 0.d0
                        H1 = ETA2(NM123) + DP(NM123)
!  RJW merged from Casey 071219: Added the following logic for 3D friction.
!  RJW modified the following for 3D friction
                        IF (C2DDI) THEN
! sb46.28sb02
!<<                     Convert Manning's N to Cd, if necessary.
                           IF (LoadManningsN) THEN
                              FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                            /((DP(NM123) + IFNLFA*ETA2(NM123)) &
                                              **(1.d0/3.d0))
                              IF (FRIC(NM123) .LT. BFCdLLimit) THEN
                                 FRIC(NM123) = BFCdLLimit
                              ENDIF
                           ENDIF
!>>
                           TKWET = FRIC(NM123)*(IFLINBF + (VELMIN/H1)* &
                                                (IFNLBF + IFHYBF* &
                                                 (1.D0 + (HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))

                           IF (TKWET .LT. 0.0001d0) TKWET = 0.0001d0
                           VEL = G*(DELETA/DELDIST)/TKWET

                        ELSEIF (C3D) THEN
!C solve for the depth averaged velocity,U, from the relation :
!C        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
!C          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
!C where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                           IF (LoadZ0B_var) THEN
                              Z0B1 = Z0B_var(NM123)
                           ELSEIF (LoadManningsN) THEN
                              Z0B1 = (DP(NM123) + IFNLFA*ETA2(NM123))*exp(-(1.0D0 + &
                                                                        ((0.41D0*(DP(NM123) + IFNLFA*ETA2(NM123))**(1.0D0/6.0D0))/ &
                                                                             (ManningsN(NM123)*sqrt(g)))))
                           ELSE
                              Z0B1 = Z0B
                           ENDIF
                           VEL = sqrt(g*H1*(DELETA/DELDIST)) &
                                 *((H1 + Z0B1)*LOG((H1 + Z0B1)/Z0B1) - H1)/(H1*0.41D0)
                        ENDIF

                        IF (VEL .GT. VELMIN) THEN
!    ....         third node met criteria and is also wet
                           NNODECODE(NM3) = 1
!  RJW merged 08/26/20008 Casey 071219: Added the following logic to obtain the correct friction.
                           IF (C2DDI) THEN
!                           TK(NM123)=FRIC(NM123)*(IFLINBF+(VEL/H1)*
                              TK(NM123) = FRIC(NM123)*(IFLINBF + (VEL/H1)* &
                                                       (IFNLBF + IFHYBF* &
                                                        (1.D0 + (HBREAK/H1)**FTHETA)** &
                                                        (FGAMMA/FTHETA)))
                           ELSEIF (C3D) THEN
                              IF (ISLIP .EQ. 0) THEN
                                 DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                        /(SIGMA(2) - SIGMA(1))
                                 BSX1(NM123) = EVTOT(1)*REAL(DUDS)
                                 BSY1(NM123) = EVTOT(1)*AIMAG(DUDS)
                                 BSX1(NM3) = EVTOT(1)*REAL(DUDS)
                                 BSY1(NM3) = EVTOT(1)*AIMAG(DUDS)
                              ENDIF
                              IF (ISLIP .NE. 0) THEN
                                 IF (ISLIP .EQ. 1) THEN
                                    KSLIP = KP
                                 ENDIF
                                 IF (ISLIP .EQ. 2) THEN
                                    KSLIP = (1.D0/ &
                                             ((1.D0/0.41D0)* &
                                              LOG((ABS(((SIGMA(2) - SIGMA(1))/(A - B))* &
                                                       (DP(NM123) + IFNLFA*ETA2(NM123))) &
                                                   + Z0B1) &
                                                  /(Z0B1))))**2.D0 &
                                            *ABS(Q(NM123, 1))
                                    IF (KSLIP .GT. 1.D0*ABS(Q(NM123, 1))) &
                                       KSLIP = 1.D0*ABS(Q(NM123, 1))
                                    IF (KSLIP .LT. 0.0025D0*ABS(Q(NM123, 1))) &
                                       KSLIP = 0.0025D0*ABS(Q(NM123, 1))
                                 ENDIF
                                 BSX1(NM123) = KSLIP*REAL(Q(NM123, 1))
                                 BSY1(NM123) = KSLIP*AIMAG(Q(NM123, 1))
                                 BSX1(NM3) = KSLIP*REAL(Q(NM123, 1))
                                 BSY1(NM3) = KSLIP*AIMAG(Q(NM123, 1))
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ELSEIF ((NODECODE(NM2) .EQ. 1) .AND. (NODECODE(NM3) .EQ. 1)) &
                     THEN
                     IF ((HTOTN2 .GE. HOFF) .AND. (HTOTN3 .GE. HOFF)) THEN
                        NM123 = NM2
                        IF (ETA2(NM3) .GT. ETA2(NM2)) NM123 = NM3
                        DELDIST = SQRT((Y(NM1) - Y(NM123))**2.D0 &
                                       + (X(NM1) - X(NM123))**2.D0)
                        DELETA = ETA2(NM123) - ETA2(NM1)
                        ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                        IF (DELETA .lt. 0.d0) DELETA = 0.d0
                        H1 = ETA2(NM123) + DP(NM123)
!  RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                        IF (C2DDI) THEN
! sb46.28sb02
!C<<                     Convert Manning's N to Cd, if necessary.
                           IF (LoadManningsN) THEN
                              FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                            /((DP(NM123) + IFNLFA*ETA2(NM123)) &
                                              **(1.d0/3.d0))
                              IF (FRIC(NM123) .LT. BFCdLLimit) THEN
                                 FRIC(NM123) = BFCdLLimit
                              ENDIF
                           ENDIF
!C>>
                           TKWET = FRIC(NM123)*(IFLINBF + (VELMIN/H1)* &
                                                (IFNLBF + IFHYBF* &
                                                 (1.D0 + (HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
                           IF (TKWET .LT. 0.0001d0) TKWET = 0.0001d0
                           VEL = G*(DELETA/DELDIST)/TKWET

                        ELSEIF (C3D) THEN
! solve for the depth averaged velocity,U, from the relation :
!        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
!          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
! where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                           IF (LoadZ0B_var) THEN
                              Z0B1 = Z0B_var(NM123)
                           ELSEIF (LoadManningsN) THEN
                              Z0B1 = (DP(NM123) + IFNLFA*ETA2(NM123))*exp(-(1.0D0 + &
                                                                        ((0.41D0*(DP(NM123) + IFNLFA*ETA2(NM123))**(1.0D0/6.0D0))/ &
                                                                             (ManningsN(NM123)*sqrt(g)))))
                           ELSE
                              Z0B1 = Z0B
                           ENDIF
                           VEL = sqrt(g*H1*(DELETA/DELDIST)) &
                                 *((H1 + Z0B1)*LOG((H1 + Z0B1)/Z0B1) - H1)/(H1*0.41D0)
                        ENDIF

                        IF (VEL .GT. VELMIN) THEN
                           NNODECODE(NM1) = 1
!  RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                           IF (C2DDI) THEN
                              TK(NM123) = FRIC(NM123)*(IFLINBF + (VEL/H1)* &
                                                       (IFNLBF + IFHYBF* &
                                                        (1.D0 + (HBREAK/H1)**FTHETA)** &
                                                        (FGAMMA/FTHETA)))
                           ELSEIF (C3D) THEN
                              IF (ISLIP .EQ. 0) THEN
                                 DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                        /(SIGMA(2) - SIGMA(1))
                                 BSX1(NM123) = EVTOT(1)*REAL(DUDS)
                                 BSY1(NM123) = EVTOT(1)*AIMAG(DUDS)
                                 BSX1(NM1) = EVTOT(1)*REAL(DUDS)
                                 BSY1(NM1) = EVTOT(1)*AIMAG(DUDS)
                              ENDIF
                              IF (ISLIP .NE. 0) THEN
                                 IF (ISLIP .EQ. 1) THEN
                                    KSLIP = KP
                                 ENDIF
                                 IF (ISLIP .EQ. 2) THEN
                                    KSLIP = (1.D0/ &
                                             ((1.D0/0.41D0)* &
                                              LOG((ABS(((SIGMA(2) - SIGMA(1))/(A - B))* &
                                                       (DP(NM123) + IFNLFA*ETA2(NM123))) &
                                                   + Z0B1) &
                                                  /(Z0B1))))**2.D0 &
                                            *ABS(Q(NM123, 1))
                                    IF (KSLIP .GT. 1.D0*ABS(Q(NM123, 1))) &
                                       KSLIP = 1.D0*ABS(Q(NM123, 1))
                                    IF (KSLIP .LT. 0.0025D0*ABS(Q(NM123, 1))) &
                                       KSLIP = 0.0025D0*ABS(Q(NM123, 1))
                                 ENDIF
                                 BSX1(NM123) = KSLIP*REAL(Q(NM123, 1))
                                 BSY1(NM123) = KSLIP*AIMAG(Q(NM123, 1))
                                 BSX1(NM1) = KSLIP*REAL(Q(NM123, 1))
                                 BSY1(NM1) = KSLIP*AIMAG(Q(NM123, 1))
                              ENDIF
                           ENDIF

                        ENDIF
                     ENDIF
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ELSEIF ((NODECODE(NM3) .EQ. 1) .AND. (NODECODE(NM1) .EQ. 1)) &
                     THEN
                     IF ((HTOTN3 .GE. HOFF) .AND. (HTOTN1 .GE. HOFF)) THEN
                        NM123 = NM3
                        IF (ETA2(NM1) .GT. ETA2(NM3)) NM123 = NM1
                        DELDIST = SQRT((Y(NM2) - Y(NM123))**2.D0 &
                                       + (X(NM2) - X(NM123))**2.D0)
                        DELETA = ETA2(NM123) - ETA2(NM2)
                        ! jgf50.60.18: Prevent numerical problems if DELETA is negative
                        IF (DELETA .lt. 0.d0) DELETA = 0.d0
                        H1 = ETA2(NM123) + DP(NM123)
!  RJW merged 08/26/2008 Casey 071219: Added the following logic for 3D friction.
                        IF (C2DDI) THEN
! sb46.28sb02
!C<<                     Convert Manning's N to Cd, if necessary.
                           IF (LoadManningsN) THEN
                              FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                            /((DP(NM123) + IFNLFA*ETA2(NM123)) &
                                              **(1.d0/3.d0))
                              IF (FRIC(NM123) .LT. BFCdLLimit) THEN
                                 FRIC(NM123) = BFCdLLimit
                              ENDIF
                           ENDIF
!C>>
                           TKWET = FRIC(NM123)*(IFLINBF + (VELMIN/H1)* &
                                                (IFNLBF + IFHYBF* &
                                                 (1.D0 + (HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
                           IF (TKWET .LT. 0.0001d0) TKWET = 0.0001d0
                           VEL = G*(DELETA/DELDIST)/TKWET
                        ELSEIF (C3D) THEN
! solve for the depth averaged velocity,U, from the relation :
!        tau=rho*g*(h+eta)*(deta/dx)=rho*Cd*|U|*U
!          U=sqrt(g*(h+eta)*(deta/dx)/Cd )
! where:  Cd=kappa^2/(ln(z+zo)/z0)^2 is the depth integrated drag coefficient
                           IF (LoadZ0B_var) THEN
                              Z0B1 = Z0B_var(NM123)
                           ELSEIF (LoadManningsN) THEN
                              Z0B1 = (DP(NM123) + IFNLFA*ETA2(NM123))*exp(-(1.0D0 + &
                                                                        ((0.41D0*(DP(NM123) + IFNLFA*ETA2(NM123))**(1.0D0/6.0D0))/ &
                                                                             (ManningsN(NM123)*sqrt(g)))))
                           ELSE
                              Z0B1 = Z0B
                           ENDIF
                           VEL = sqrt(g*H1*(DELETA/DELDIST)) &
                                 *((H1 + Z0B1)*LOG((H1 + Z0B1)/Z0B1) - H1)/(H1*0.41D0)
                        ENDIF

                        IF (VEL .GT. VELMIN) THEN
                           NNODECODE(NM2) = 1
!  RJW merged 08/26/2008 Casey 071219: Added the following logic to obtain the correct friction.
                           IF (C2DDI) THEN
                              TK(NM123) = FRIC(NM123)*(IFLINBF + (VEL/H1)* &
                                                       (IFNLBF + IFHYBF* &
                                                        (1.D0 + (HBREAK/H1)**FTHETA)** &
                                                        (FGAMMA/FTHETA)))
                           ELSEIF (C3D) THEN
                              IF (ISLIP .EQ. 0) THEN
                                 DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                        /(SIGMA(2) - SIGMA(1))
                                 BSX1(NM123) = EVTOT(1)*REAL(DUDS)
                                 BSY1(NM123) = EVTOT(1)*AIMAG(DUDS)
                                 BSX1(NM2) = EVTOT(1)*REAL(DUDS)
                                 BSY1(NM2) = EVTOT(1)*AIMAG(DUDS)
                              ENDIF
                              IF (ISLIP .NE. 0) THEN
                                 IF (ISLIP .EQ. 1) THEN
                                    KSLIP = KP
                                 ENDIF
                                 IF (ISLIP .EQ. 2) THEN
                                    KSLIP = (1.D0/ &
                                             ((1.D0/0.41D0)* &
                                              LOG((ABS(((SIGMA(2) - SIGMA(1))/(A - B))* &
                                                       (DP(NM123) + IFNLFA*ETA2(NM123))) &
                                                   + Z0B1) &
                                                  /(Z0B1))))**2.D0 &
                                            *ABS(Q(NM123, 1))
                                    IF (KSLIP .GT. 1.D0*ABS(Q(NM123, 1))) &
                                       KSLIP = 1.D0*ABS(Q(NM123, 1))
                                    IF (KSLIP .LT. 0.0025D0*ABS(Q(NM123, 1))) &
                                       KSLIP = 0.0025D0*ABS(Q(NM123, 1))
                                 ENDIF
                                 BSX1(NM123) = KSLIP*REAL(Q(NM123, 1))
                                 BSY1(NM123) = KSLIP*AIMAG(Q(NM123, 1))
                                 BSX1(NM2) = KSLIP*REAL(Q(NM123, 1))
                                 BSY1(NM2) = KSLIP*AIMAG(Q(NM123, 1))
                              ENDIF
                           ENDIF

                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
! WET...
! WET...Nodal Wetting Criteria W2a
! WET...
               NBNCTOT = NIBNODECODE(NM1) + NIBNODECODE(NM2) + NIBNODECODE(NM3)
               NIBCNT(NM1) = NIBCNT(NM1) + NBNCTOT
               NIBCNT(NM2) = NIBCNT(NM2) + NBNCTOT
               NIBCNT(NM3) = NIBCNT(NM3) + NBNCTOT

            ENDDO

            if (subdomainOn .and. enforceBN .eq. 1) call enforceWDcb() ! NCSU Subdomain
            if (subdomainOn .and. enforceBN .eq. 2) call enforceWDob() ! NCSU Subdomain

! wet...
! WET... ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
! *******************************************************************************************
! tcm v50.66.01 -- This is an additional test for wetting only when time varying
!                  bathymetry is used and is only performed during the period of
!                  bathymetry evolution.
!
            IF ((NDDT .NE. 0) .AND. (IT .LE. BTIME_END + 1)) THEN
               DO I = 1, NE
                  NM1 = NM(I, 1)
                  NM2 = NM(I, 2)
                  NM3 = NM(I, 3)
                  NCTOT = NODECODE(NM1) + NODECODE(NM2) + NODECODE(NM3)
                  IF (NCTOT .lt. 3) THEN !If not wet from previous time step
                     NCTOT = NNODECODE(NM1) + NNODECODE(NM2) + NNODECODE(NM3)
                     if (NCTOT .lt. 3) then !if not alreay made wet for this time step
                        ETAN1 = ETA2(NM1)
                        ETAN2 = ETA2(NM2)
                        ETAN3 = ETA2(NM3)
                        HTOTN1 = DP(NM1) + ETA2(NM1)
                        HTOTN2 = DP(NM2) + ETA2(NM2)
                        HTOTN3 = DP(NM3) + ETA2(NM3)

!                    if all nodes have a depth greater than or equal to
!                    hoff = 1.2*H0, then make the element wet
                        IF ((HTOTN1 .GE. HOFF) .AND. (HTOTN2 .GE. HOFF) .AND. &
                            (HTOTN3 .GE. HOFF)) THEN
                           !THE ELEMENT SHOULD BE WET, SO WET THE DRY NODES
                           !  Make Node 1 Wet and set parameters
                           IF (NNODECODE(NM1) .NE. 1) THEN !node 1
                              NNODECODE(NM1) = 1
                              NM123 = NM1
                              IF (C2DDI) THEN
!C<<                           Convert Manning's N to Cd, if necessary.
                                 IF (LoadManningsN) THEN
                                    FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                                  /((DP(NM123) + IFNLFA*ETA2(NM123)) &
                                                    **(1.d0/3.d0))
                                    IF (FRIC(NM123) .LT. BFCdLLimit) THEN
                                       FRIC(NM123) = BFCdLLimit
                                    ENDIF
                                 ENDIF
                              ENDIF
                              VEL = VELMIN
                              H1 = HTOTN1
                              IF (C2DDI) THEN
                                 TK(NM123) = FRIC(NM123)*(IFLINBF + (VEL/H1)* &
                                                          (IFNLBF + IFHYBF* &
                                                           (1.D0 + (HBREAK/H1)**FTHETA)** &
                                                           (FGAMMA/FTHETA)))
                              ELSEIF (C3D) THEN
                                 IF (ISLIP .EQ. 0) THEN
                                    DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                           /(SIGMA(2) - SIGMA(1))
                                    BSX1(NM123) = EVTOT(1)*REAL(DUDS)
                                    BSY1(NM123) = EVTOT(1)*AIMAG(DUDS)
                                    BSX1(NM3) = EVTOT(1)*REAL(DUDS)
                                    BSY1(NM3) = EVTOT(1)*AIMAG(DUDS)
                                 ENDIF
                                 IF (ISLIP .NE. 0) THEN
                                    IF (ISLIP .EQ. 1) THEN
                                       KSLIP = KP
                                    ENDIF
                                    IF (ISLIP .EQ. 2) THEN
                                       KSLIP = (1.D0/ &
                                                ((1.D0/0.41D0)* &
                                                 LOG((ABS(((SIGMA(2) - SIGMA(1))/(A - B))* &
                                                          (DP(NM123) + IFNLFA*ETA2(NM123))) &
                                                      + Z0B1) &
                                                     /(Z0B1))))**2.D0 &
                                               *ABS(Q(NM123, 1))
                                       IF (KSLIP .GT. 1.D0*ABS(Q(NM123, 1))) &
                                          KSLIP = 1.D0*ABS(Q(NM123, 1))
                                       IF (KSLIP .LT. 0.0025D0*ABS(Q(NM123, 1))) &
                                          KSLIP = 0.0025D0*ABS(Q(NM123, 1))
                                    ENDIF
                                    BSX1(NM123) = KSLIP*REAL(Q(NM123, 1))
                                    BSY1(NM123) = KSLIP*AIMAG(Q(NM123, 1))
                                    BSX1(NM3) = KSLIP*REAL(Q(NM123, 1))
                                    BSY1(NM3) = KSLIP*AIMAG(Q(NM123, 1))
                                 ENDIF
                              ENDIF
                           ENDIF !end node 1

                           !  Make Node 2 Wet and set parameters
                           IF (NNODECODE(NM2) .NE. 1) THEN
                              NNODECODE(NM2) = 1
                              NM123 = NM2
                              IF (C2DDI) THEN
!C<<                        Convert Manning's N to Cd, if necessary.
                                 IF (LoadManningsN) THEN
                                    FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                                  /((DP(NM123) + IFNLFA*ETA2(NM123)) &
                                                    **(1.d0/3.d0))
                                    IF (FRIC(NM123) .LT. BFCdLLimit) THEN
                                       FRIC(NM123) = BFCdLLimit
                                    ENDIF
                                 ENDIF
                              ENDIF
                              VEL = VELMIN
                              H1 = HTOTN2
                              IF (C2DDI) THEN
                                 TK(NM123) = FRIC(NM123)*(IFLINBF + (VEL/H1)* &
                                                          (IFNLBF + IFHYBF* &
                                                           (1.D0 + (HBREAK/H1)**FTHETA)** &
                                                           (FGAMMA/FTHETA)))
                              ELSEIF (C3D) THEN
                                 IF (ISLIP .EQ. 0) THEN
                                    DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                           /(SIGMA(2) - SIGMA(1))
                                    BSX1(NM123) = EVTOT(1)*REAL(DUDS)
                                    BSY1(NM123) = EVTOT(1)*AIMAG(DUDS)
                                    BSX1(NM1) = EVTOT(1)*REAL(DUDS)
                                    BSY1(NM1) = EVTOT(1)*AIMAG(DUDS)
                                 ENDIF
                                 IF (ISLIP .NE. 0) THEN
                                    IF (ISLIP .EQ. 1) THEN
                                       KSLIP = KP
                                    ENDIF
                                    IF (ISLIP .EQ. 2) THEN
                                       KSLIP = (1.D0/ &
                                                ((1.D0/0.41D0)* &
                                                 LOG((ABS(((SIGMA(2) - SIGMA(1))/(A - B))* &
                                                          (DP(NM123) + IFNLFA*ETA2(NM123))) &
                                                      + Z0B1) &
                                                     /(Z0B1))))**2.D0 &
                                               *ABS(Q(NM123, 1))
                                       IF (KSLIP .GT. 1.D0*ABS(Q(NM123, 1))) &
                                          KSLIP = 1.D0*ABS(Q(NM123, 1))
                                       IF (KSLIP .LT. 0.0025D0*ABS(Q(NM123, 1))) &
                                          KSLIP = 0.0025D0*ABS(Q(NM123, 1))
                                    ENDIF
                                    BSX1(NM123) = KSLIP*REAL(Q(NM123, 1))
                                    BSY1(NM123) = KSLIP*AIMAG(Q(NM123, 1))
                                    BSX1(NM1) = KSLIP*REAL(Q(NM123, 1))
                                    BSY1(NM1) = KSLIP*AIMAG(Q(NM123, 1))
                                 ENDIF
                              ENDIF
                           ENDIF !node 2

                           !  Make Node 3 Wet and set parameters
                           IF (NNODECODE(NM3) .NE. 1) THEN
                              NNODECODE(NM3) = 1
                              NM123 = NM3
                              IF (C2DDI) THEN
!C<<                         Convert Manning's N to Cd, if necessary.
                                 IF (LoadManningsN) THEN
                                    FRIC(NM123) = g*ManningsN(NM123)**2.d0 &
                                                  /((DP(NM123) + IFNLFA*ETA2(NM123)) &
                                                    **(1.d0/3.d0))
                                    IF (FRIC(NM123) .LT. BFCdLLimit) THEN
                                       FRIC(NM123) = BFCdLLimit
                                    ENDIF
                                 ENDIF
                              ENDIF
                              VEL = VELMIN
                              H1 = HTOTN3
                              IF (C2DDI) THEN
                                 TK(NM123) = FRIC(NM123)*(IFLINBF + (VEL/H1)* &
                                                          (IFNLBF + IFHYBF* &
                                                           (1.D0 + (HBREAK/H1)**FTHETA)** &
                                                           (FGAMMA/FTHETA)))
                              ELSEIF (C3D) THEN
                                 IF (ISLIP .EQ. 0) THEN
                                    DUDS = (Q(NM123, 2) - Q(NM123, 1)) &
                                           /(SIGMA(2) - SIGMA(1))
                                    BSX1(NM123) = EVTOT(1)*REAL(DUDS)
                                    BSY1(NM123) = EVTOT(1)*AIMAG(DUDS)
                                    BSX1(NM2) = EVTOT(1)*REAL(DUDS)
                                    BSY1(NM2) = EVTOT(1)*AIMAG(DUDS)
                                 ENDIF
                                 IF (ISLIP .NE. 0) THEN
                                    IF (ISLIP .EQ. 1) THEN
                                       KSLIP = KP
                                    ENDIF
                                    IF (ISLIP .EQ. 2) THEN
                                       KSLIP = (1.D0/ &
                                                ((1.D0/0.41D0)* &
                                                 LOG((ABS(((SIGMA(2) - SIGMA(1))/(A - B))* &
                                                          (DP(NM123) + IFNLFA*ETA2(NM123))) &
                                                      + Z0B1) &
                                                     /(Z0B1))))**2.D0 &
                                               *ABS(Q(NM123, 1))
                                       IF (KSLIP .GT. 1.D0*ABS(Q(NM123, 1))) &
                                          KSLIP = 1.D0*ABS(Q(NM123, 1))
                                       IF (KSLIP .LT. 0.0025D0*ABS(Q(NM123, 1))) &
                                          KSLIP = 0.0025D0*ABS(Q(NM123, 1))
                                    ENDIF
                                    BSX1(NM123) = KSLIP*REAL(Q(NM123, 1))
                                    BSY1(NM123) = KSLIP*AIMAG(Q(NM123, 1))
                                    BSX1(NM2) = KSLIP*REAL(Q(NM123, 1))
                                    BSY1(NM2) = KSLIP*AIMAG(Q(NM123, 1))
                                 ENDIF
                              ENDIF
                           ENDIF !node3

                        ENDIF !ALL DEPTHS GREATER THAN HOFF
                     ENDIF !IF NNODECODE SUM LESS THAN 3
                  ENDIF ! IF NODECODE SUM LESS THAN 3
               ENDDO !LOOP OVER ELEMENTS

               if (subdomainOn .and. enforceBN .eq. 1) call enforceWDcb() ! NCSU Subdomain
               if (subdomainOn .and. enforceBN .eq. 2) call enforceWDob() ! NCSU Subdomain


            ENDIF !IT TIME VARYING BATHYMETRY AND WITHIN CHANGE TIME

! ... END OF ADDITIONAL WETTING FOR TIME VARYING BATHYMETRY
! WET..
! WET... ELEMENTAL WETTING CRITERIA WETBATHYCHANGE
! *******************************************************************************************

! WET...
! WET...Nodal Wetting Criteria W2b
! WET...Check for adjacent nodes and force nodes wet when attached
! WET...to receiving barrier nodes
! WET...
            DO I = 1, NP
               IF ((NIBCNT(I) .GT. 0) .AND. (NNODECODE(I) .EQ. 0)) THEN
                  NNODECODE(I) = 1
               ENDIF
            ENDDO

! WET...
! WET...END WET/DRY SECTION - PART 2
! WET...

! WET...
! WET...START WET/DRY SECTION  - PART 3
! WET...Elemental drying criteria DE1
! WET...This is an elemental check section designed to avoid artificial wetting of
! WET....of control sections
! WET...All elements where downhill flow originates from a barely wet node
! WET....into wet nodes are forced inactive; the only exception is receiving
! WET....overtopped barrier nodes
! WET...
            DO I = 1, NE
               NM1 = NM(I, 1)
               NM2 = NM(I, 2)
               NM3 = NM(I, 3)
               NBNCTOT = NIBCNT(NM1)*NIBCNT(NM2)*NIBCNT(NM3)
               IF (NBNCTOT .EQ. 0) THEN !No barrier/pipe receiving nodes in this elem
                  ETAN1 = ETA2(NM1)
                  ETAN2 = ETA2(NM2)
                  ETAN3 = ETA2(NM3)
                  HTOTN1 = DP(NM1) + ETA2(NM1)
                  HTOTN2 = DP(NM2) + ETA2(NM2)
                  HTOTN3 = DP(NM3) + ETA2(NM3)
                  ! jgf52.08.08: Analyst can eliminate noff from
                  ! consideration in fort.15 file.
                  if (noffActive .eqv. .true.) then
! ..ABC pattern
                     IF ((ETAN1 .GE. ETAN2) .AND. (ETAN2 .GT. ETAN3)) THEN
                        IF ((HTOTN1 .LT. HOFF) .OR. (HTOTN2 .LT. HOFF)) NOFF(I) = 0
                     ENDIF
                     IF ((ETAN2 .GE. ETAN3) .AND. (ETAN3 .GT. ETAN1)) THEN
                        IF ((HTOTN2 .LT. HOFF) .OR. (HTOTN3 .LT. HOFF)) NOFF(I) = 0
                     ENDIF
                     IF ((ETAN3 .GE. ETAN1) .AND. (ETAN1 .GT. ETAN2)) THEN
                        IF ((HTOTN3 .LT. HOFF) .OR. (HTOTN1 .LT. HOFF)) NOFF(I) = 0
                     ENDIF
! ..ACB pattern
                     IF ((ETAN1 .GE. ETAN3) .AND. (ETAN3 .GT. ETAN2)) THEN
                        IF ((HTOTN1 .LT. HOFF) .OR. (HTOTN3 .LT. HOFF)) NOFF(I) = 0
                     ENDIF
                     IF ((ETAN2 .GE. ETAN1) .AND. (ETAN1 .GT. ETAN3)) THEN
                        IF ((HTOTN2 .LT. HOFF) .OR. (HTOTN1 .LT. HOFF)) NOFF(I) = 0
                     ENDIF
                     IF ((ETAN3 .GE. ETAN2) .AND. (ETAN2 .GT. ETAN1)) THEN
                        IF ((HTOTN3 .LT. HOFF) .OR. (HTOTN2 .LT. HOFF)) NOFF(I) = 0
                     ENDIF
                  endif
               ENDIF
            ENDDO

! WET...
! WET...END WET/DRY SECTION  - PART 3
! WET...

! WET...
! WET...START WET/DRY SECTION PART 4 - NODAL DRYING LOOP D2
! WET...Update number of active elements (MJU) and the total area (TotalArea) connected
! WET...to a node. If these are zero, the node is landlocked and should be dried.
! WET...These depend on NNODECODE which varies during the time step
! WET...
            DO I = 1, NP
               MJU(I) = 0
               TotalArea(I) = 0.d0
            ENDDO
            DO IE = 1, NE
               NM1 = NM(IE, 1)
               NM2 = NM(IE, 2)
               NM3 = NM(IE, 3)
               NC1 = NNODECODE(NM1)
               NC2 = NNODECODE(NM2)
               NC3 = NNODECODE(NM3)

               NCEle = NC1*NC2*NC3*NOFF(IE)
               AreaEle = NCEle*Areas(IE)/2.d0
               MJU(NM1) = MJU(NM1) + NCEle
               MJU(NM2) = MJU(NM2) + NCEle
               MJU(NM3) = MJU(NM3) + NCEle
               TotalArea(NM1) = TotalArea(NM1) + AreaEle
               TotalArea(NM2) = TotalArea(NM2) + AreaEle
               TotalArea(NM3) = TotalArea(NM3) + AreaEle
            ENDDO

! jjwnote - looks like this is used later in momentum equations
! jjwnote - this has implications on making this into a subroutine

            DO I = 1, NP
               IF ((NNODECODE(I) .EQ. 1) .AND. (MJU(I) .EQ. 0)) THEN
                  NNODECODE(I) = 0
               ENDIF
               IF (MJU(I) .EQ. 0) MJU(I) = 1 !Because MJU is also used to solve Mom Eq. !Eliminate this?
            ENDDO

!     WET...
!     WET...END WET/DRY SECTION - PART 4
!     WET...

! jjwnote - may have to pass TotalArea and mju as well

            if (subdomainOn .and. enforceBN .eq. 1) call enforceWDcb() ! NCSU Subdomain
            if (subdomainOn .and. enforceBN .eq. 2) call enforceWDob() ! NCSU Subdomain

! WET...
! WET...WET/DRY SECTION - PART 5 - RESET NODECODE USING NNODECODE
! WET...Check to see if any wetting occurred & update NODECODE
! WET...Note, NCCHANGE=0 set near the beginning of GWCE subroutine
! WET...
            DO I = 1, NP
               IF (NNODECODE(I) .NE. NODECODE(I)) THEN
                  NODECODE(I) = NNODECODE(I)
                  NCCHANGE = NCCHANGE + 1
               ENDIF
            ENDDO
! WET...
! WET...END WET/DRY SECTION - PART 5
! WET...

! WET...
! WET...WET/DRY SECTION - PART 6
! WET...Check to see if any NOFF changed requiring the matrix to be reset
! WET...Note, NCCHANGE=0 set near the beginning of GWCE subroutine
! WET...
            DO I = 1, NE
               IF (NOFF(I) .NE. NOFFOLD(I)) NCCHANGE = NCCHANGE + 1
            ENDDO
! WET...
! WET... jgf45.06 If there has been any wetting or drying in any
! WET... of the subdomains, the NCCHANGE flag will be activated on all
! WET... of the subdomains, to prevent them from getting out of sync
! WET... with their MPI calls as some reset the GWCE and others do not.
! WET...
! WET...END WET/DRY SECTION - PART 6
! WET...

! ....
!-----------------------------------------------------------------------
         end subroutine computeWettingAndDrying
!-----------------------------------------------------------------------

      end module wetDry

