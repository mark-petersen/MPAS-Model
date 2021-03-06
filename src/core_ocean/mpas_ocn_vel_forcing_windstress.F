! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_windstress
!
!> \brief MPAS ocean wind stress
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routine for computing 
!>  tendencies from wind stress.  
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_windstress

   use mpas_grid_types
   use mpas_configure

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

   public :: ocn_vel_forcing_windstress_tend, &
             ocn_vel_forcing_windstress_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: windStressOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_windstress_tend
!
!> \brief   Computes tendency term from wind stress
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    15 September 2011
!> \details 
!>  This routine computes the wind stress tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_windstress_tend(grid, normalVelocityForcing, layerThicknessEdge, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocityForcing    !< Input: wind stress

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge     !< Input: thickness at edge

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

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

      integer :: iEdge, nEdgesSolve, k
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: edgeMask

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.windStressOn) return

      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      edgeMask => grid % edgeMask % array

      do iEdge=1,nEdgesSolve
        ! efficiency note: it would be nice to avoid this
        ! if within a do.  This could be done with
        ! k =  max(maxLevelEdgeTop(iEdge),1)
        ! and then tend_u(1,iEdge) is just not used for land cells.

        do k = 1,min(maxLevelEdgeTop(iEdge),1)
           ! forcing in top layer only
           tend(k,iEdge) =  tend(k,iEdge) + edgeMask(k, iEdge) * (normalVelocityForcing(k,iEdge) / config_density0 / layerThicknessEdge(k,iEdge))
        enddo

      enddo


   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_windstress_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_windstress_init
!
!> \brief   Initializes ocean wind stress forcing
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine initializes quantities related to wind stress 
!>  in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_windstress_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      windStressOn = .true.

      if(config_disable_vel_windstress) windStressOn = .false.

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_windstress_init!}}}

!***********************************************************************

end module ocn_vel_forcing_windstress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
