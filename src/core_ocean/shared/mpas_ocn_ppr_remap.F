! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_equation_of_state
!
!> \brief MPAS ocean Piecewise Polynomial Reconstruction driver
!> \author Mark Petersen
!> \date   July 2019
!> \details
!>  This module contains the main driver routines for calling
!>  the Piecewise Polynomial Reconstruction and related quantities.
!
!-----------------------------------------------------------------------

module ocn_equation_of_state

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

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

   public :: ocn_equation_of_state_density, &
             ocn_equation_of_state_init, &
             ocn_freezing_temperature, &
             ocn_freezing_temperature_salinity_deriv

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ppr_1d
!
!> \brief   Computes density using an ocean Piecewise Polynomial Reconstruction
!> \author  Mark Petersen
!> \date    July 2019
!> \details
!>  This routine calls the Piecewise Polynomial Reconstruction to compute a density
!
!-----------------------------------------------------------------------

include "../ppr/ppr_1d.f90"

!***********************************************************************

end module ocn_equation_of_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
