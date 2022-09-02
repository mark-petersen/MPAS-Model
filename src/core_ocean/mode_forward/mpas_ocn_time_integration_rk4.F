! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_rk4
!
!> \brief MPAS ocean RK4 Time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the RK4 time integration routine.
!
!-----------------------------------------------------------------------

module ocn_time_integration_rk4

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_tendency
   use ocn_diagnostics
   use ocn_gm

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average_coupled

   use ocn_effective_density_in_land_ice

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_time_integrator_rk4

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_rk4
!
!> \brief MPAS ocean RK4 Time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine integrates one timestep (dt) using an RK4 time integrator.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_rk4(domain, dt)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step using
   !   4th order Runge-Kutta
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:))
   !                 plus mesh meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      real (kind=RKIND), intent(in) :: dt !< Input: timestep

      integer :: iCell, iEdge, k, i, err
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: tracersTendPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: scratchPool

      integer :: rk_step

      !type (mpas_pool_type), pointer :: nextProvisPool, prevProvisPool

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

      real (kind=RKIND) :: coef
      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfEdge, vertDiffTopOfCell

      ! Dimensions
      integer, pointer :: nCells, nEdges, nVertLevels, num_tracers

      ! Config options
      logical, pointer :: config_prescribe_velocity, config_prescribe_thickness
      logical, pointer :: config_filter_btr_mode, config_use_freq_filtered_thickness
      logical, pointer :: config_use_GM
      logical, pointer :: config_use_cvmix_kpp
      logical, pointer :: config_use_tracerGroup
      logical, pointer :: config_disable_thick_all_tend
      logical, pointer :: config_disable_vel_all_tend
      logical, pointer :: config_disable_tr_all_tend
      real (kind=RKIND), pointer :: config_mom_del4
      real (kind=RKIND), pointer :: config_drying_min_cell_height
      character (len=StrKIND), pointer :: config_land_ice_flux_mode

      ! State indices
      integer, pointer :: indexTemperature
      integer, pointer :: indexSalinity

      ! Diagnostics Indices
      integer, pointer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
      integer, pointer :: indexSSHGradientZonal, indexSSHGradientMeridional

      ! Mesh array pointers
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop
      real (kind=RKIND), dimension(:), pointer :: bottomDepth

      ! Provis Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityNew, layerThicknessNew
!      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupProvis

      ! Tend Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityTend, &
                                                    layerThicknessTend
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupTend

      ! Diagnostics Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge
      real (kind=RKIND), dimension(:,:), pointer :: vertAleTransportTop
      real (kind=RKIND), dimension(:,:), pointer :: normalTransportVelocity, normalGMBolusVelocity
      real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ
      real (kind=RKIND), dimension(:,:), pointer :: velocityZonal, velocityMeridional
      real (kind=RKIND), dimension(:), pointer :: gradSSH
      real (kind=RKIND), dimension(:), pointer :: gradSSHX, gradSSHY, gradSSHZ
      real (kind=RKIND), dimension(:), pointer :: gradSSHZonal, gradSSHMeridional
      real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, sshGradient

      ! State Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityCur
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur
      real (kind=RKIND), dimension(:), pointer :: sshCur, sshNew

      real (kind=RKIND), dimension(:,:,:), pointer :: tracerGroup, tracersCur, tracersNew

      ! Diagnostics Field Pointers
      type (field1DReal), pointer :: boundaryLayerDepthField, effectiveDensityField
      type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField

      ! State/Tend Field Pointers
      type (field2DReal), pointer :: normalVelocityField, layerThicknessField
      type (field3DReal), pointer :: tracersGroupField

      ! Tracer Group Iteartion
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      ! Tidal boundary condition
      logical, pointer :: config_use_tidal_forcing
      character (len=StrKIND), pointer :: config_tidal_forcing_type
      real (kind=RKIND), dimension(:), pointer :: tidalInputMask, tidalBCValue
      real (kind=RKIND), dimension(:,:), pointer :: restingThickness
      real (kind=RKIND) :: totalDepth, currentTime
      real (kind=RKIND), pointer :: daysSinceStartOfSim
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocitySolution
      integer, dimension(:,:), pointer :: boundaryEdge

      ! Get config options
      call mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
      call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)
      call mpas_pool_get_config(domain % configs, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
      call mpas_pool_get_config(domain % configs, 'config_use_GM', config_use_GM)
      call mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
      call mpas_pool_get_config(domain % configs, 'config_disable_vel_all_tend', config_disable_vel_all_tend)
      call mpas_pool_get_config(domain % configs, 'config_disable_thick_all_tend', config_disable_thick_all_tend)
      call mpas_pool_get_config(domain % configs, 'config_disable_tr_all_tend', config_disable_tr_all_tend)

      !
      ! Initialize time_levs(2) with state at current time
      ! Initialize first RK state
      ! Couple tracers time_levs(2) with layerThickness in time-levels
      ! Initialize RK weights
      !
      !block => domain % blocklist
      !do while (associated(block))
      !   call mpas_pool_get_subpool(block % structs, 'state', statePool)
      !   call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      !   call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      !   call mpas_pool_create_pool(statePool)

      !   call mpas_pool_clone_pool(statePool, statePool, 1)
      !   call mpas_pool_add_subpool(block % structs, 'state', statePool)

      !   call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      !   call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

      !   call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
      !   call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
      !   call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
      !   call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

      !   call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      !   call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      !   !$omp parallel
      !   !$omp do schedule(runtime) private(k)
      !   do iEdge = 1, nEdges
      !      do k = 1, maxLevelEdgeTop(iEdge)
      !         normalVelocityNew(k, iEdge) = normalVelocityCur(k, iEdge)
      !      end do
      !   end do
      !   !$omp end do

      !   !$omp do schedule(runtime) private(k)
      !   do iCell = 1, nCells
      !      do k = 1, maxLevelCell(iCell)
      !         layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell)
      !      end do
      !   end do
      !   !$omp end do
      !   !$omp end parallel

      !   block => block % next
      !end do

      !block => domain % blocklist
      !do while(associated(block))
      !   if (associated(block % prev)) then
      !      call mpas_pool_get_subpool(block % prev % structs, 'state', prevProvisPool)
      !   else
      !      nullify(prevProvisPool)
      !   end if

      !   if (associated(block % next)) then
      !      call mpas_pool_get_subpool(block % next % structs, 'state', nextProvisPool)
      !   else
      !      nullify(nextProvisPool)
      !   end if

      !   call mpas_pool_get_subpool(block % structs, 'state', statePool)

      !   if (associated(prevProvisPool) .and. associated(nextProvisPool)) then
      !      call mpas_pool_link_pools(statePool, prevProvisPool, nextProvisPool)
      !   else if (associated(prevProvisPool)) then
      !      call mpas_pool_link_pools(statePool, prevProvisPool)
      !   else if (associated(nextProvisPool)) then
      !      call mpas_pool_link_pools(statePool, nextPool=nextProvisPool)
      !   else
      !      call mpas_pool_link_pools(statePool)
      !   end if

      !   call mpas_pool_link_parinfo(block, statePool)

      !   block => block % next
      !end do

      ! Fourth-order Runge-Kutta, solving dy/dt = f(t,y) is typically written as follows
      ! where h = delta t is the large time step.  Here f(t,y) is the right hand side,
      ! called the tendencies in the code below.
      ! k_1 = h f(t_n        , y_n)
      ! k_2 = h f(t_n + 1/2 h, y_n + 1/2 k_1)
      ! k_3 = h f(t_n + 1/2 h, y_n + 1/2 k_2)
      ! k_4 = h f(t_n +     h, y_n +     k_3)
      ! y_{n+1} = y_n + 1/6 k_1 + 1/3 k_2 + 1/3 k_3 + 1/6 k_4

      ! in index notation:
      ! k_{j+1} = h f(t_n + a_j h, y_n + a_j k_j)
      ! y_{n+1} = y_n + sum ( b_j k_j )

      ! The coefficients of k_j are b_j = (1/6, 1/3, 1/3, 1/6) and are
      ! initialized here as delta t * b_j:

! mrp forward backward:
      rk_weights(1) = dt
      !rk_weights(1) = dt/6.0_RKIND
      rk_weights(2) = dt/3.0_RKIND
      rk_weights(3) = dt/3.0_RKIND
      rk_weights(4) = dt/6.0_RKIND

      ! The a_j coefficients of h in the computation of k_j are typically written (0, 1/2, 1/2, 1).
      ! However, in the algorithm below we pre-compute the state for the tendency one iteration early.
      ! That is, on j=1 (rk_step=1, below) we pre-compute y_n + 1/2 k_1 and save it in state.
      ! Then we compute 1/6 k_1 and add it to state % time_levs(2).
      ! That is why the coefficients of h are one index early in the following, i.e.
      ! a = (1/2, 1/2, 1)

      rk_substep_weights(1) = dt ! mrp/2.0_RKIND
      rk_substep_weights(2) = dt/2.0_RKIND
      rk_substep_weights(3) = dt
      rk_substep_weights(4) = dt ! a_4 only used for ALE step, otherwise it is skipped.

      call mpas_timer_start("RK4-main loop")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN RK loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do rk_step = 1,1 !mrp , 4

        call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

        call mpas_timer_start("RK4-diagnostic halo update")
        call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
        call mpas_timer_stop("RK4-diagnostic halo update")

! mrp forward backward:
! compute tendencies
        ! Compute tendencies for velocity, thickness, and tracers.
        ! In RK4 notation, we are computing the right hand side f(t,y),
        ! which is the same as k_j / h.
        call mpas_timer_start("RK4 vel/thick tendency computations")
        block => domain % blocklist
        do while (associated(block))
           call ocn_time_integrator_rk4_compute_vel_tends( block, dt, rk_substep_weights(rk_step), err )
           call ocn_time_integrator_rk4_compute_thick_tends( block, dt, rk_substep_weights(rk_step), err )
           block => block % next
        end do
        call mpas_timer_stop("RK4 vel/thick tendency computations")

! mrp forward backward:
! update halo
        ! Update halos for prognostic variables.
        call mpas_timer_start("RK4 vel/thick prognostic halo update")
        call mpas_dmpar_field_halo_exch(domain, 'tendNormalVelocity')
        call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')
        call mpas_timer_stop("RK4 vel/thick prognostic halo update")

! mrp forward backward:
! remove this. It is an extra diagnostics in the middle
        !call mpas_timer_start("RK4-update diagnostic variables")
        !if (rk_step < 4) then
        !   block => domain % blocklist
        !   do while (associated(block))
        !      call ocn_time_integrator_rk4_diagnostic_update(block, dt, rk_substep_weights(rk_step), err)
        !      block => block % next
        !   end do
        !end if
        !call mpas_timer_stop("RK4-update diagnostic variables")

! mrp forward backward:
! add RHS
        ! Accumulate update.
        ! In RK4 notation, we are computing b_j k_j and adding it to an accumulating sum so that we have
        !    y_{n+1} = y_n + sum ( b_j k_j )
        ! after the fourth iteration.

        call mpas_timer_start("RK4-accumulate update")

        block => domain % blocklist
        do while (associated(block))
           call ocn_time_integrator_rk4_accumulate_update(block, rk_weights(rk_step), err)

           block => block % next
        end do

        call mpas_timer_stop("RK4-accumulate update")

      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END RK loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_stop("RK4-main loop")

      !
      !  A little clean up at the end: rescale tracer fields and compute diagnostics for new state
      !
      call mpas_timer_start("RK4-cleanup phase")

! mrp forward backward:
! diagnostics
      block => domain % blocklist
      do while(associated(block))
        call ocn_time_integrator_rk4_cleanup(block, dt, err)

        block => block % next
      end do

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)

         call mpas_pool_get_array(diagnosticsPool, 'daysSinceStartOfSim',daysSinceStartOfSim)
         ! daysSinceStartOfSim is at old time
         currentTime = daysSinceStartOfSim*86400.0_RKIND + dt
         call ocn_diagnostic_solve(dt, currentTime, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)

         ! ------------------------------------------------------------------
         ! Accumulating various parameterizations of the transport velocity
         ! ------------------------------------------------------------------
         !!$omp parallel
         !!$omp do schedule(runtime) 
         !do iEdge = 1, nEdges
         !   normalTransportVelocity(:, iEdge) = normalVelocityNew(:, iEdge)
         !end do
         !!$omp end do
         !!$omp end parallel

         block => block % next
      end do

      call mpas_timer_stop("RK4-cleanup phase")

      !block => domain % blocklist
      !do while(associated(block))
      !   call mpas_pool_get_subpool(block % structs, 'state', statePool)

      !   call mpas_pool_destroy_pool(statePool)

      !   call mpas_pool_remove_subpool(block % structs, 'state')
      !   block => block % next
      !end do

   end subroutine ocn_time_integrator_rk4!}}}

   subroutine ocn_time_integrator_rk4_compute_vel_tends(block, dt, &
     rkSubstepWeight, err)!{{{

      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), intent(in) :: rkSubstepWeight
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
      !real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity

      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

      call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 1, dt)

   end subroutine ocn_time_integrator_rk4_compute_vel_tends!}}}

   subroutine ocn_time_integrator_rk4_compute_thick_tends(block, dt, rkSubstepWeight, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), intent(in) :: rkSubstepWeight
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
!      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity
      real (kind=RKIND), dimension(:, :), pointer ::  normalVelocityNew

      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

      call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool, verticalMeshPool)

   end subroutine ocn_time_integrator_rk4_compute_thick_tends!}}}

   subroutine ocn_time_integrator_rk4_compute_tends(block, dt, rkWeight, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), intent(in) :: rkWeight
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: statePool, diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool

      real (kind=RKIND), dimension(:), pointer :: sshCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, normalVelocityCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity

      logical, pointer :: config_filter_btr_mode

      err = 0

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

      call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 1, dt)

      call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool, verticalMeshPool)

   end subroutine ocn_time_integrator_rk4_compute_tends!}}}

!   subroutine ocn_time_integrator_rk4_diagnostic_update(block, dt, rkSubstepWeight, err)!{{{
!      type (block_type), intent(in) :: block
!      real (kind=RKIND), intent(in) :: dt
!      real (kind=RKIND), intent(in) :: rkSubstepWeight
!      integer, intent(out) :: err
!
!      logical, pointer :: config_prescribe_velocity, config_prescribe_thickness, config_use_GM
!
!      integer, pointer :: nCells, nEdges, nVertLevels
!      integer :: iCell, iEdge, k
!
!      type (mpas_pool_type), pointer :: statePool, tendPool, meshPool, scratchPool
!      type (mpas_pool_type), pointer :: diagnosticsPool, statePool, forcingPool
!      type (mpas_pool_type), pointer :: tracersPool, tracersTendPool, tracersPool
!
!      real (kind=RKIND), dimension(:, :), pointer :: normalVelocityCur, normalVelocityNew, normalVelocityTend
!      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur, layerThicknessNew, layerThicknessTend
!      integer, dimension(:,:), pointer :: boundaryEdge
!      real (kind=RKIND) :: totalDepth, currentTime
!      real (kind=RKIND), pointer :: daysSinceStartOfSim
!
!      real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroupCur, tracersGroupProvis, tracersGroupTend
!
!      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop
!
!      logical, pointer :: config_use_tracerGroup
!      type (mpas_pool_iterator_type) :: groupItr
!      character (len=StrKIND) :: modifiedGroupName
!      character (len=StrKIND) :: configName
!
!      err = 0
!
!      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
!      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
!      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
!
!      call mpas_pool_get_subpool(block % structs, 'state', statePool)
!      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
!      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
!      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
!      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
!      call mpas_pool_get_subpool(block % structs, 'state', statePool)
!      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
!
!      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
!      call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
!
!      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
!
!      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
!      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
!
!      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 1)
!      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 1)
!
!      call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
!      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
!
!      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
!      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
!
!      !$omp parallel
!      !$omp do schedule(runtime) private(k)
!      do iEdge = 1, nEdges
!         do k = 1, maxLevelEdgeTop(iEdge)
!            normalVelocityNew(k, iEdge) = normalVelocityCur(k, iEdge) + rkSubstepWeight &
!                                           * normalVelocityTend(k, iEdge)
!         end do
!      end do
!      !$omp end do
!
!      !$omp do schedule(runtime) private(k)
!      do iCell = 1, nCells
!         do k = 1, maxLevelCell(iCell)
!            layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell) + rkSubstepWeight &
!                                           * layerThicknessTend(k, iCell)
!         end do
!      end do
!      !$omp end do
!      !$omp end parallel
!
!      call mpas_pool_get_array(diagnosticsPool, 'daysSinceStartOfSim',daysSinceStartOfSim)
!      ! daysSinceStartOfSim is at old time
!      ! rkSubstepWeight is either dt or dt/2, depending on RK4 stage
!      currentTime = daysSinceStartOfSim*86400.0_RKIND + rkSubstepWeight
!      call ocn_diagnostic_solve(dt, currentTime, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 1)
!
!   end subroutine ocn_time_integrator_rk4_diagnostic_update!}}}

   subroutine ocn_time_integrator_rk4_accumulate_update(block, rkWeight, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: rkWeight
      integer, intent(out) :: err

      integer, pointer :: nCells, nEdges
      integer :: iCell, iEdge, k

      type (mpas_pool_type), pointer :: statePool, tendPool, meshPool, diagnosticsPool

      real (kind=RKIND), dimension(:, :), pointer :: normalVelocityNew, normalVelocityTend
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessNew, layerThicknessTend
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityCur
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur

      integer, dimension(:), pointer :: maxLevelCell

      logical, pointer :: config_use_tracerGroup
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName

      err = 0

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

      call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)


      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
         do k = 1, maxLevelCell(iCell)
            layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell) + rkWeight * layerThicknessTend(k, iCell)
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime) 
      do iEdge = 1, nEdges
         normalVelocityNew(:, iEdge) = normalVelocityCur(:, iEdge) + rkWeight * normalVelocityTend(:, iEdge)
      end do
      !$omp end do
      !$omp end parallel

   end subroutine ocn_time_integrator_rk4_accumulate_update!}}}

   subroutine ocn_time_integrator_rk4_cleanup(block, dt, err)!{{{
      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      integer, pointer :: nCells, nEdges, indexTemperature, indexSalinity
      integer :: iCell, iEdge, k

      type (mpas_pool_type), pointer :: statePool, meshPool, forcingPool
      type (mpas_pool_type), pointer :: diagnosticsPool, scratchPool
      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessNew, normalVelocityNew
      !real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity, normalGMBolusVelocity
      !real (kind=RKIND), dimension(:, :, :), pointer :: tracersGroupNew

      integer, dimension(:), pointer :: maxLevelCell

      logical, pointer :: config_use_tracerGroup
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName
      real (kind=RKIND) :: totalDepth, currentTime
      real (kind=RKIND), pointer :: daysSinceStartOfSim

      logical, pointer :: config_use_GM

      err = 0

      call mpas_pool_get_config(block % configs, 'config_use_GM', config_use_GM)

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)

      !call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
      !call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)

      !call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      !call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      !call mpas_pool_begin_iteration(tracersPool)
      !do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
      !   if ( groupItr % memberType == MPAS_POOL_FIELD ) then
      !      call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)
      !      if ( associated(tracersGroupNew) ) then
      !         !$omp parallel
      !         !$omp do schedule(runtime) private(k)
      !         do iCell = 1, nCells
      !           do k = 1, maxLevelCell(iCell)
      !             tracersGroupNew(:, k, iCell) = tracersGroupNew(:, k, iCell) / layerThicknessNew(k, iCell)
      !           end do
      !         end do
      !         !$omp end do
      !         !$omp end parallel
      !      end if
      !   end if
      !end do

      call mpas_pool_get_array(diagnosticsPool, 'daysSinceStartOfSim',daysSinceStartOfSim)
      ! daysSinceStartOfSim is at old time
      currentTime = daysSinceStartOfSim*86400.0_RKIND + dt
      call ocn_diagnostic_solve(dt, currentTime, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)

      !call ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, forcingPool, scratchPool, err, 2)

   end subroutine ocn_time_integrator_rk4_cleanup!}}}

end module ocn_time_integration_rk4

! vim: foldmethod=marker
