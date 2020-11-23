! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_verification_cases
!
!> \brief MPAS ocean diagnostics driver
!> \author Mark Petersen
!> \date   23 September 2011
!> \details
!>  This module contains the routines for computing
!>  verification cases, and other quantities such as vertAleTransportTop.
!
!-----------------------------------------------------------------------

module ocn_verification_cases

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_stream_manager
   use mpas_io_units

   use ocn_constants
   use ocn_config

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

   public :: ocn_verification_solution, &
             ocn_verification_cases_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_verification_solution
!
!> \brief   Computes verification cases
!> \author  Mark Petersen
!> \date    23 September 2011
!> \details
!>  This routine computes the verification cases for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_verification_solution(currentTime, statePool, meshPool, diagnosticsPool, tracersPool)

      real (kind=RKIND), intent(in) :: currentTime ! current time sec
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(inout) :: diagnosticsPool  !< Input: diagnostic fields derived from State
      type (mpas_pool_type), intent(in) :: tracersPool !< Input: tracer fields

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j
      integer :: boundaryMask, velMask, err, nCells, nEdges, nVertices
      integer, pointer  :: nVertLevels, vertexDegree
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray, nVerticesArray

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
        maxLevelVertexBot
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, &
        verticesOnEdge, edgesOnEdge, edgesOnVertex,boundaryCell, kiteIndexOnCell, &
        verticesOnCell, edgeSignOnVertex, edgeSignOnCell, edgesOnCell, edgeMask

      integer :: timeLevel

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, xEdge, yEdge, angleEdge
      real (kind=RKIND), dimension(:), pointer :: sshSolution, velocityMeridionalSolution, velocityZonalSolution
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocitySolution
      real (kind=RKIND) :: f0, beta, maxDepth, phaseSpeed, Rinv, etaHat, kY, v
      real (kind=RKIND) :: xMin,xMax,yMin,yMax,yMinGlobal,yMaxGlobal,xMinGlobal,xMaxGlobal
      real (kind=RKIND), pointer :: yPlanarDomainWidth, xEdgeGlobalMin

      if (config_verification_enable) then
       call mpas_pool_get_array(diagnosticsPool, 'sshSolution', sshSolution)
       call mpas_pool_get_array(diagnosticsPool, 'normalVelocitySolution', normalVelocitySolution)
       call mpas_pool_get_array(diagnosticsPool, 'velocityMeridionalSolution', velocityMeridionalSolution)
       call mpas_pool_get_array(diagnosticsPool, 'velocityZonalSolution', velocityZonalSolution)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
       call mpas_pool_get_array(meshPool, 'xCell', xCell)
       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
       call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
       call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
       call mpas_pool_get_array(diagnosticsPool, 'yPlanarDomainWidth', yPlanarDomainWidth)
       call mpas_pool_get_array(diagnosticsPool, 'xEdgeGlobalMin', xEdgeGlobalMin)

       nEdges = nEdgesArray( size(nEdgesArray) )
       nCells = nCellsArray( size(nCellsArray) )
       maxDepth = 3000.0_RKIND ! H [m]
       select case (config_verification_case)
       case (1) ! coastal Kelvin wave
          f0 = 1.0e-4_RKIND ! coriolis parameter [1/s]
          beta = 0.0 ! beta parameter
          phaseSpeed = sqrt(gravity*maxDepth) ! c=sqrt(gH)
          Rinv = f0 / phaseSpeed ! R = sqrt(gH)/f0
          etaHat = 0.001
          kY = 2.0_RKIND*pii/yPlanarDomainWidth

          do iCell = 1, nCells
              sshSolution(iCell) = - maxDepth &
                 * etaHat*sin(Ky*(yCell(iCell) + phaseSpeed*currentTime)) &
                 * exp(-(xCell(iCell) - xEdgeGlobalMin) * Rinv)
              velocityMeridionalSolution(iCell) = phaseSpeed &
                 * etaHat*sin(Ky*(yCell(iCell) + phaseSpeed*currentTime)) &
                 * exp(-(xCell(iCell) - xEdgeGlobalMin) * Rinv)
          end do
          do iEdge = 1, nEdges
              ! u = 0.0
              v = phaseSpeed &
                 * etaHat*sin(Ky*(yEdge(iEdge) + phaseSpeed*currentTime)) &
                 * exp(-(xEdge(iEdge) - xEdgeGlobalMin) * Rinv)
             ! angleEdge=0: edge normal points in positive x
             ! angleEdge=pi/2: edge normal points in positive y
             ! normalVelocity = -u*sin(angleEdge) + v*cos(angleEdge)
             !normalVelocity(:, iEdge) = u*cos(angleEdge(iEdge)) + v*sin(angleEdge(iEdge))
             normalVelocitySolution(:, iEdge) = v*sin(angleEdge(iEdge))
          end do
       end select
      endif

   end subroutine ocn_verification_solution!}}}

!***********************************************************************
!
!  routine ocn_verification_cases_init
!
!> \brief   Initializes flags used within diagnostics routines.
!> \author  Mark Petersen
!> \date    4 November 2011
!> \details
!>  This routine initializes flags related to quantities computed within
!>  other diagnostics routines.
!
!-----------------------------------------------------------------------
   subroutine ocn_verification_cases_init(err)!{{{
      integer, intent(out) :: err !< Output: Error flag

      err = 0

    end subroutine ocn_verification_cases_init!}}}


!***********************************************************************

end module ocn_verification_cases

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
