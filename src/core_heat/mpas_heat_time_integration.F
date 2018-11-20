! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module heat_time_integration

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_log

   use heat_constants

   contains

   subroutine heat_timestep(domain, dt, timeStamp)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool

      character (len=StrKIND), pointer :: xtime
      character (len=StrKIND), pointer :: config_time_integration

      call mpas_pool_get_config(domain % configs, 'config_time_integration', config_time_integration)

      if (trim(config_time_integration) == 'forward_Euler') then
         call heat_forward_Euler(domain, dt)
      else
         call mpas_log_write('Unknown time integration option '//trim(config_time_integration), MPAS_LOG_ERR)
         call mpas_log_write('Currently, only ''forward_Euler'' is supported.', MPAS_LOG_CRIT)
      end if

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)

         call mpas_pool_get_array(statePool, 'xtime', xtime, 2)
         xtime = timeStamp 
         block => block % next
      end do

   end subroutine heat_timestep


   subroutine heat_forward_euler(domain, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step using 
   !   4th order Runge-Kutta
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      integer :: iCell, k
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: tendPool

      integer :: rk_step

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels

      real (kind=RKIND), dimension(:,:,:), pointer :: tracersOld, tracersNew, tracersProvis, tracersTend

      integer, pointer :: config_test_case
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc4


      call mpas_pool_get_config(domain % configs, 'config_test_case', config_test_case)
      call mpas_pool_get_config(domain % configs, 'config_h_mom_eddy_visc4', config_h_mom_eddy_visc4)

! --- compute tendencies
       block => domain % blocklist
       do while (associated(block))
          call mpas_pool_get_subpool(block % structs, 'state', statePool)
          call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
          call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

          call heat_compute_tracer_tend(tendPool, statePool, meshPool, 1)

          block => block % next
       end do

! --- update halos for prognostic variables
       call mpas_dmpar_field_halo_exch(domain, 'tend_tracers')

       block => domain % blocklist
       do while (associated(block))
          call mpas_pool_get_subpool(block % structs, 'state', statePool)
          call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
          call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

          call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)
          call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)
          call mpas_pool_get_array(tendPool, 'tracers', tracersTend)

          do iCell = 1, nCells
             do k = 1, nVertLevels
                tracersNew(:,k,iCell) = tracersNew(:,k,iCell) + dt * tracersTend(:,k,iCell)
             end do
          end do
          block => block % next
       end do


   end subroutine heat_forward_Euler


   subroutine heat_compute_tracer_tend(tendPool, statePool, meshPool, timeLevelIn)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed scalar tendencies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (mpas_pool_type), intent(inout) :: tendPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: meshPool
      integer, intent(in), optional :: timeLevelIn

      integer :: iCell, iEdge, k, iTracer, cell1, cell2, i
      real (kind=RKIND) :: flux, tracer_edge, r
      real (kind=RKIND) :: invAreaCell1, invAreaCell2, tracer_turb_flux
      integer, dimension(:,:), pointer :: boundaryEdge
      real (kind=RKIND), dimension(:,:), allocatable :: boundaryMask
      real (kind=RKIND), dimension(:,:,:), allocatable:: delsq_tracer
      
      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers, tracer_tend
      integer, dimension(:,:), pointer :: cellsOnEdge, boundaryCell
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND) :: coef_3rd_order
      real (kind=RKIND), dimension(:,:), pointer :: u, h_edge

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, nTracers
      integer :: timeLevel

      real (kind=RKIND), pointer :: config_h_tracer_eddy_diff2, config_h_tracer_eddy_diff4

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      endif

      call mpas_pool_get_config(heatConfigs, 'config_h_tracer_eddy_diff2', config_h_tracer_eddy_diff2)
      call mpas_pool_get_config(heatConfigs, 'config_h_tracer_eddy_diff4', config_h_tracer_eddy_diff4)

      call mpas_pool_get_array(statePool, 'u', u, timeLevel)
      call mpas_pool_get_array(statePool, 'h_edge', h_edge, timeLevel)
      call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)

      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'deriv_two', deriv_two)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)
      call mpas_pool_get_array(meshPool, 'boundaryEdge', boundaryEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

      call mpas_pool_get_array(tendPool, 'tracers', tracer_tend)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
      call mpas_pool_get_dimension(meshPool, 'nTracers', nTracers)

      tracer_tend(:,:,:) = 0.0

      ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
      !
      if ( config_h_tracer_eddy_diff2 > 0.0 ) then

         !
         ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
         !
         allocate(boundaryMask(nVertLevels, nEdges+1))
         boundaryMask = 1.0
         where(boundaryEdge.eq.1) boundaryMask=0.0

         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invAreaCell1 = 1.0/areaCell(cell1)
            invAreaCell2 = 1.0/areaCell(cell2)

            do k = 1, nVertLevels
              do iTracer = 1,  nTracers
                 ! \kappa_2 \nabla \phi on edge
                 tracer_turb_flux = config_h_tracer_eddy_diff2 &
                    *( tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge)

                 ! div(h \kappa_2 \nabla \phi) at cell center
                 flux = dvEdge(iEdge) * h_edge(k,iEdge) * tracer_turb_flux * boundaryMask(k, iEdge)
                 tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) + flux * invAreaCell1
                 tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) - flux * invAreaCell2
              end do
            end do

         end do

        deallocate(boundaryMask)

      end if

      !
      ! tracer tendency: del4 horizontal tracer diffusion, &
      !    div(h \kappa_4 \nabla [div(h \nabla \phi)])
      !
      if ( config_h_tracer_eddy_diff4 > 0.0 ) then

         !
         ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
         !
         allocate(boundaryMask(nVertLevels, nEdges+1))
         boundaryMask = 1.0
         where(boundaryEdge.eq.1) boundaryMask=0.0

         allocate(delsq_tracer(nTracers, nVertLevels, nCells+1))

         delsq_tracer(:,:,:) = 0.

         ! first del2: div(h \nabla \phi) at cell center
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k = 1, nVertLevels
              do iTracer = 1, nTracers
                 delsq_tracer(iTracer,k,cell1) = delsq_tracer(iTracer,k,cell1) &
                    + dvEdge(iEdge) * h_edge(k,iEdge) * (tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge) * boundaryMask(k,iEdge)
                 delsq_tracer(iTracer,k,cell2) = delsq_tracer(iTracer,k,cell2) &
                    - dvEdge(iEdge) * h_edge(k,iEdge) * (tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge) * boundaryMask(k,iEdge)
              end do
            end do

         end do

         do iCell = 1, nCells
            r = 1.0 / areaCell(iCell)
            do k = 1, nVertLevels
            do iTracer = 1, nTracers
               delsq_tracer(iTracer,k,iCell) = delsq_tracer(iTracer,k,iCell) * r
            end do
            end do
         end do

         ! second del2: div(h \nabla [delsq_tracer]) at cell center
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invAreaCell1 = 1.0 / areaCell(cell1)
            invAreaCell2 = 1.0 / areaCell(cell2)

            do k = 1, nVertLevels
            do iTracer = 1, nTracers
               tracer_turb_flux = config_h_tracer_eddy_diff4 * (delsq_tracer(iTracer,k,cell2) - delsq_tracer(iTracer,k,cell1)) / dcEdge(iEdge)
               flux = dvEdge(iEdge) * tracer_turb_flux
               tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) - flux * invAreaCell1 * boundaryMask(k,iEdge)
               tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) + flux * invAreaCell2 * boundaryMask(k,iEdge)
            end do
            enddo

         end do

         deallocate(delsq_tracer)
         deallocate(boundaryMask)

      end if

   end subroutine heat_compute_tracer_tend

end module heat_time_integration
