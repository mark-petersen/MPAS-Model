! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_verification
!
!> \brief MPAS ocean initialize case -- verification
!> \author Mark Petersen
!> \date   October 2020
!> \details
!>  This module contains the routines for initializing the
!>  verification initial condition
!>
!
!-----------------------------------------------------------------------

module ocn_init_verification

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar

   use ocn_constants
   use ocn_config
   use ocn_init_vertical_grids
   use ocn_init_cell_markers

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

   public :: ocn_init_setup_verification, &
             ocn_init_validate_verification

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_verification
!
!> \brief   Setup for this initial condition
!> \author  Mark Petersen
!> \date    October 2020
!> \details
!>  This routine sets up the initial conditions for this case.
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_verification(domain, iErr)!{{{

    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: diagnosticsPool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: verticalMeshPool

    ! local variables
    integer :: iCell, iEdge, k, idx
    real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
    real (kind=RKIND) :: localVar1, localVar2
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations
    real (kind=RKIND) :: f0, beta, maxDepth, phaseSpeed, Rinv, etaHat, kY, T0, S0, v, currentTime
    logical :: flatBottom

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: maxLevelCell
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, xEdge, yEdge, yVertex, dcEdge, angleEdge
    real (kind=rkind), dimension(:), pointer :: fCell, fEdge, fVertex
    real (kind=rkind), dimension(:), pointer :: vertCoordMovementWeights, bottomDepth
    real (kind=RKIND), dimension(:), pointer :: refBottomDepth, refZMid, refLayerThickness
    real (kind=RKIND), dimension(:), pointer :: ssh
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness, normalVelocity
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
    real (kind=RKIND), pointer :: yPlanarDomainWidth, xEdgeGlobalMin

    iErr = 0

    if(config_init_configuration .ne. trim('verification')) return

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

     if ( on_a_sphere ) call mpas_log_write('The verification configuration can only be applied ' &
         // 'to a planar mesh. Exiting...', MPAS_LOG_CRIT)

    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

    !--------------------------------------------------------------------
    ! Use this section to make boundaries non-periodic
    !--------------------------------------------------------------------

    ! Initalize min/max values to large positive and negative values
    yMin = 1.0E10_RKIND
    yMax = -1.0E10_RKIND
    xMin = 1.0E10_RKIND
    xMax = -1.0E10_RKIND
    dcEdgeMin = 1.0E10_RKIND

!    ! Determine local min and max values.
!    block_ptr => domain % blocklist
!    do while(associated(block_ptr))
!       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
!       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

!
!       call mpas_pool_get_array(meshPool, 'xCell', xCell)
!       call mpas_pool_get_array(meshPool, 'yCell', yCell)
!       call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
!       call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
!
!       xCellMin = minval(xCell(1:nCellsSolve))
!       xCellMax = maxval(xCell(1:nCellsSolve))
!       yCellMin = minval(yCell(1:nCellsSolve))
!       yCellMax = maxval(yCell(1:nCellsSolve))
!       xEdgeMin = minval(xEdge(1:nEdgesSolve))
!       xEdgeMax = maxval(xEdge(1:nEdgesSolve))
!       yEdgeMin = minval(yEdge(1:nEdgesSolve))
!       yEdgeMax = maxval(yEdge(1:nEdgesSolve))
!
!       block_ptr => block_ptr % next
!    end do
!
!    ! Determine global min and max values.
!    call mpas_dmpar_min_real(domain % dminfo, xCellMin, xCellGlobalMin)
!    call mpas_dmpar_max_real(domain % dminfo, xCellMax, xCellGlobalMax)
!    call mpas_dmpar_min_real(domain % dminfo, yCellMin, yCellGlobalMin)
!    call mpas_dmpar_max_real(domain % dminfo, yCellMax, yCellGlobalMax)
!    call mpas_dmpar_min_real(domain % dminfo, xEdgeMin, xEdgeGlobalMin)
!    call mpas_dmpar_max_real(domain % dminfo, xEdgeMax, xEdgeGlobalMax)
!    call mpas_dmpar_min_real(domain % dminfo, yEdgeMin, yEdgeGlobalMin)
!    call mpas_dmpar_max_real(domain % dminfo, yEdgeMax, yEdgeGlobalMax)

    !--------------------------------------------------------------------
    ! Use this section to set initial values
    !--------------------------------------------------------------------

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
       call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

       call mpas_pool_get_array(meshPool, 'xCell', xCell)
       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
       call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
       call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
       call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
       call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
       call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
       call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'fCell', fCell)
       call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
       call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

       call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
       call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
       call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)

       call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
       call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
       call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)

       call mpas_pool_get_array(diagnosticsPool, 'yPlanarDomainWidth', yPlanarDomainWidth)
       call mpas_pool_get_array(diagnosticsPool, 'xEdgeGlobalMin', xEdgeGlobalMin)

       maxDepth = 3000.0_RKIND ! H [m]

       ! Set refBottomDepth and refZMid
       do k = 1, nVertLevels
          refBottomDepth(k) = maxDepth * interfaceLocations(k+1)
          refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * maxDepth
       end do
       ! Compute refLayerThickness and refZMid
       call ocn_compute_layerThickness_zMid_from_bottomDepth(refLayerThickness,refZMid, &
               refBottomDepth,refBottomDepth(nVertLevels), nVertLevels,nVertLevels,iErr)

       ! Set vertCoordMovementWeights
       vertCoordMovementWeights(:) = 1.0_RKIND

       select case (config_verification_case)
       case (1) ! coastal Kelvin wave
          f0 = 1.0e-4_RKIND ! coriolis parameter [1/s]
          beta = 0.0 ! beta parameter
          flatBottom = .true.
          phaseSpeed = sqrt(gravity*maxDepth) ! c=sqrt(gH)
          Rinv = f0 / phaseSpeed ! R = sqrt(gH)/f0
          etaHat = 0.001
          kY = 2.0_RKIND*pii/yPlanarDomainWidth
          T0 = 20.0_RKIND ! constant temperature
          S0 = 35.0_RKIND ! constant salinity
          currentTime = 0.0_RKIND ! start at time zero

          do iCell = 1, nCellsSolve
              ssh(iCell) = - maxDepth &
                 * etaHat*sin(Ky*(yCell(iCell) + phaseSpeed*currentTime)) &
                 * exp(-(xCell(iCell) - xEdgeGlobalMin) * Rinv)
          end do
          do iEdge = 1, nEdgesSolve
              ! u = 0.0
              v = phaseSpeed &
                 * etaHat*sin(Ky*(yEdge(iEdge) + phaseSpeed*currentTime)) &
                 * exp(-(xEdge(iEdge) - xEdgeGlobalMin) * Rinv) 
             ! angleEdge=0: edge normal points in positive x
             ! angleEdge=pi/2: edge normal points in positive y
             ! normalVelocity = -u*sin(angleEdge) + v*cos(angleEdge)
             !normalVelocity(:, iEdge) = u*cos(angleEdge(iEdge)) + v*sin(angleEdge(iEdge))
             normalVelocity(:, iEdge) = v*sin(angleEdge(iEdge))
          end do
       end select

       if (flatBottom) then
          bottomDepth(:) = maxDepth
          maxLevelCell(:) = nVertLevels
       endif

       do iCell = 1, nCellsSolve
          do k = 1, maxLevelCell(iCell)
             layerThickness(k, iCell) = refLayerThickness(k) + ssh(iCell)/maxLevelCell(iCell)
             restingThickness(k, iCell) = refLayerThickness(k)
             activeTracers(index_temperature, k, iCell) = T0
             activeTracers(index_salinity, k, iCell) = S0
          end do
       end do

       fCell(:) = f0 + beta*yCell(:)
       fEdge(:) = f0 + beta*yEdge(:)
       fVertex(:) = f0 + beta*yVertex(:)

       block_ptr => block_ptr % next
    end do

    deallocate(interfaceLocations)
    !--------------------------------------------------------------------

  end subroutine ocn_init_setup_verification!}}}

!***********************************************************************
!
!  routine interp_vert_velocity_to_zlevel
!
!> \brief   Interpolate the vertical velcity to a z level
!> \author  Phillip Wolfram
!> \date    05/08/2014
!> \details
!>  This routine interpolates the vertical velocity to a particular
!>  z-level.
!
!-----------------------------------------------------------------------
!real (kind=RKIND) function CoastalKelvinWaveFunctionalForm( etaHat, kY, y) {{{
!
!     !-----------------------------------------------------------------
!     ! input variables
!     !-----------------------------------------------------------------
!     real (kind=RKIND), intent(in) :: etaHat ! coefficient
!     real (kind=RKIND), intent(in) :: kY ! coefficient
!     real (kind=RKIND), intent(in) :: y ! coefficient
!
!     !-----------------------------------------------------------------
!     ! output variables
!     !-----------------------------------------------------------------
!     !real (kind=RKIND), intent(out) :: interp_vert_velocity_to_zlevel
!
!     !-----------------------------------------------------------------
!     ! local variables
!     !-----------------------------------------------------------------
!     real (kind=RKIND) :: alpha
!
!     eta = etaHat*sin(kY*y)

!***********************************************************************
!
!  routine ocn_init_validate_verification
!
!> \brief   Validation for this initial condition
!> \author  Mark Petersen
!> \date    October 2020
!> \details
!>  This routine validates the configuration options for this case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_verification(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_verification_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('verification')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_verification_vert_levels', config_verification_vert_levels)

      if(config_vert_levels <= 0 .and. config_verification_vert_levels > 0) then
         config_vert_levels = config_verification_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for verification. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_verification!}}}


!***********************************************************************

end module ocn_init_verification

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
