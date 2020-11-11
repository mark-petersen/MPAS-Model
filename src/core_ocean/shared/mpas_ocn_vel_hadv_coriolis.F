! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hadv_coriolis
!
!> \brief MPAS ocean horizontal momentum mixing driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routine for computing
!>  tendencies from the horizontal advection and coriolis force.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_hadv_coriolis

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
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

   public :: ocn_vel_hadv_coriolis_tend, &
             ocn_vel_hadv_coriolis_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: gradKEOn ! Turn on gradient of kinetic energy term
   integer :: relVorOn ! Turn on relative vorticity term
   integer :: planVorOn ! Turn on planetary vorticity term

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hadv_coriolis_tend
!
!> \brief   Computes tendency term for horizontal advection and coriolis force
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details

!>  This routine computes the horizontal advection and coriolis and
!>  advection tendencies for momentum based on current state.

!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hadv_coriolis_tend(meshPool, normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, &
                                    layerThicknessEdge, normalVelocity, kineticEnergyCell, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalizedRelativeVorticityEdge, &!< Input: relative vorticity over thickness, on an edge
         normalizedPlanetaryVorticityEdge,&!< Input: planetary vorticity over thickness, on an edge
         layerThicknessEdge,&!< Input: Thickness on edge
         normalVelocity,&    !< Input: Horizontal velocity
         kineticEnergyCell   !< Input: Kinetic Energy

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnEdge
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnEdge, edgeMask
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge
      real (kind=RKIND), dimension(:), pointer :: dcEdge

      integer :: j, k
      integer :: cell1, cell2, iEdge, eoe, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nEdgesArray
      real (kind=RKIND) :: workVorticity, invLength, edgeWeight, r_tmp
      real (kind=RKIND), dimension(:), allocatable :: qArr

      err = 0

      if (relVorOn.eq.0.and.gradKEOn.eq.0) return

      call mpas_timer_start("coriolis")

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      nEdges = nEdgesArray( 1 )

      allocate( qArr(nVertLevels) )

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(cell1, cell2, invLength, k, qArr, j, eoe, edgeWeight, workVorticity)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         invLength = 1.0_RKIND / dcEdge(iEdge)

         do k = 1, maxLevelEdgeTop(iEdge)
            qArr(k) = 0.0_RKIND
         end do

         do j = 1, nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(j, iEdge)
            edgeWeight = weightsOnEdge(j, iEdge)

            do k = 1, maxLevelEdgeTop(iEdge)
               workVorticity = 0.5_RKIND &
                  * (  relVorOn*normalizedRelativeVorticityEdge(k, iEdge) + planVorOn * normalizedPlanetaryVorticityEdge(k, iEdge) &
                     + relVorOn*normalizedRelativeVorticityEdge(k, eoe)   + planVorOn * normalizedPlanetaryVorticityEdge(k, eoe))
               qArr(k) = qArr(k) + edgeWeight * normalVelocity(k, eoe) * workVorticity * layerThicknessEdge(k, eoe)
            end do

         end do

         do k = 1, maxLevelEdgeTop(iEdge)
            tend(k, iEdge) = tend(k, iEdge) &
               + edgeMask(k, iEdge) &
                 * ( qArr(k) &
                     - gradKEOn*(  kineticEnergyCell(k, cell2) &
                                 - kineticEnergyCell(k, cell1) ) * invLength )
         end do

      end do
      !$omp end do
      !$omp end parallel

      deallocate( qArr )

      call mpas_timer_stop("coriolis")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hadv_coriolis_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hadv_coriolis_init
!
!> \brief   Initializes ocean momentum horizontal advection and
!>  coriolis tendencies
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hadv_coriolis_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! Output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      gradKEOn = 1
      relVorOn = 1
      planVorOn = 1

      ! In (eta+f)/h, do not add f/h
      if ( config_disable_vel_coriolis ) planVorOn = 0

      ! In (eta+f)/h, do not add eta/h. Also do not add grad(KE)
      if ( config_disable_vel_hadv ) then
         relVorOn = 0
         gradKEOn = 0
      end if

      if ( trim( config_time_integrator ) == 'split_explicit' &
        .or. trim( config_time_integrator ) == 'unsplit_explicit') then
         ! In (eta+f)/h, do not add f/h
         ! For split explicit, Coriolis tendency 
         ! is added separately to the momentum tendencies.
         planVorOn = 0
      end if


   !--------------------------------------------------------------------

   end subroutine ocn_vel_hadv_coriolis_init!}}}

!***********************************************************************

end module ocn_vel_hadv_coriolis

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
