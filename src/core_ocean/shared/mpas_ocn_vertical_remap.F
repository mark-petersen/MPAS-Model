! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vertical_remap
!
!> \brief MPAS ocean Piecewise Polynomial Reconstruction driver
!> \author Mark Petersen
!> \date   July 2019
!> \details
!>  This module contains the main driver routines for calling
!>  the Piecewise Polynomial Reconstruction and related quantities.
!
!-----------------------------------------------------------------------

module ocn_vertical_remap

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

   use ppr_1d

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

   public :: ocn_test_ppr_1d
!             ocn_freezing_temperature_salinity_deriv

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

subroutine ocn_test_ppr_1d()


        integer, parameter :: npos = 31 ! no. edge (old grid)
        integer, parameter :: ntmp = 23 ! no. edge (new grid)
        integer, parameter :: nvar = 1  ! no. variables to remap
        integer, parameter :: ndof = 1  ! no. FV DoF per cell
        integer :: ipos

    !------------------------------ position of cell edges !
        real*8  :: xpos(npos),xtmp(ntmp)
        real*8  :: xmid

    !-------------------------------- finite-volume arrays !

    !   Arrays represent a "block" of finite-volume tracers
    !   to remap. The 1st dim. is the no. of DoF per cell,
    !   NDOF=1 is a standard finite-volume scheme where the
    !   data is specified as cell means. NDOF>1 is reserved
    !   for future use with DG-style schemes. NVAR is the
    !   number of tracers to remap. Processing tracers in a
    !   batch is typically more efficient than one-by-one.
    !   The last dim. is the no. cells (layers) in the grid.

        real*8  :: init(ndof,nvar,npos-1)
        real*8  :: ftmp(ndof,nvar,ntmp-1)
        real*8  :: fdat(ndof,nvar,npos-1)

    !------------------------------ method data-structures !
        type(rmap_work) :: work
        type(rmap_opts) :: opts
        type(rcon_ends) :: bc_l(nvar)
        type(rcon_ends) :: bc_r(nvar)

    !------------------------------ define a simple domain !

        call linspace(0.d0,1.d0,npos,xpos)
        call linspace(0.d0,1.d0,ntmp,xtmp)

    !------------------------------ setup some simple data !

        do ipos = +1, npos-1

            xmid = xpos(ipos+0) * 0.5d+0 &
    &            + xpos(ipos+1) * 0.5d+0

            init(1,1,ipos) = xmid ** 2

        end do

    !------------------------------ specify method options !

        opts%edge_meth = p3e_method     ! 3rd-order edge interp.
        opts%cell_meth = ppm_method     ! PPM method in cells
        opts%cell_lims = null_limit     ! no slope limiter

    !------------------------------ set BC.'s at endpoints !

        bc_l%bcopt = bcon_loose         ! "loose" = extrapolate
        bc_r%bcopt = bcon_loose

    !------------------------------ init. method workspace !

        call work%init(npos,nvar,opts)

    !------------------------------ re-map back-and-forth: !

        fdat = init

        do ipos = +1, +1000

    !------------------------------ re-map from dat-to-tmp !

        call rmap1d(npos,ntmp,nvar,ndof, &
        &           xpos,xtmp,fdat,ftmp, &
        &           bc_l,bc_r,work,opts)

    !------------------------------ re-map from tmp-to-dat !

        call rmap1d(ntmp,npos,nvar,ndof, &
        &           xtmp,xpos,ftmp,fdat, &
        &           bc_l,bc_r,work,opts)

        end do

    !------------------------------ clear method workspace !

        call work%free()

    !------------------------------ dump results to stdout !

        print*,"Cell data: [INIT] [RMAP] "

        do ipos = +1, npos-1

            print *, init(1,1,ipos) &
        &          , fdat(1,1,ipos)

        end do

        print*,"Conservation defect := " &
        &     , sum(init) - sum(fdat)

end subroutine ocn_test_ppr_1d

!***********************************************************************

end module ocn_vertical_remap

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
