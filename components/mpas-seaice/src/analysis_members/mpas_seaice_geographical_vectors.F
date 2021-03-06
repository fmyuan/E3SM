!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_geographical_vectors
!
!> \brief MPAS sea ice analysis mode member: geographical_vectors
!> \authorAdrian K. Turner
!> \date  14th October 2015
!> \details
!>  MPAS sea ice analysis mode member: geographical_vectors
!>
!-----------------------------------------------------------------------

module seaice_geographical_vectors

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

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

   public :: seaice_bootstrap_geographical_vectors, &
             seaice_init_geographical_vectors, &
             seaice_precompute_geographical_vectors, &
             seaice_compute_geographical_vectors, &
             seaice_restart_geographical_vectors, &
             seaice_finalize_geographical_vectors

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_geographical_vectors
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    10th November 2015
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_geographical_vectors(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

    end subroutine seaice_bootstrap_geographical_vectors!}}}

!***********************************************************************
!
!  routine seaice_init_geographical_vectors
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th October 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_geographical_vectors(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

    end subroutine seaice_init_geographical_vectors!}}}

!***********************************************************************
!
!  routine seaice_precompute_geographical_vectors
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th October 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_geographical_vectors(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine seaice_precompute_geographical_vectors!}}}

!***********************************************************************
!
!  routine seaice_compute_geographical_vectors
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th October 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_geographical_vectors(domain, instance, timeLevel, err)!{{{

     use seaice_mesh, only: &
          seaice_latlon_vector_rotation_backward

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           velocitySolverPool, &
           geographicalVectorsAMPool

      ! mesh quantities
      real(kind=RKIND), dimension(:), pointer :: &
           latVertex, &
           lonVertex, &
           xVertex, &
           yVertex, &
           zVertex

      ! other quantities
      real(kind=RKIND), pointer :: &
           sphereRadius

      logical, pointer :: &
           rotateCartesianGrid

      ! non-geographical arrays
      real(kind=RKIND), dimension(:), pointer :: &
           uVelocity, &
           vVelocity, &
           stressDivergenceU, &
           stressDivergenceV, &
           airStressVertexU, &
           airStressVertexV, &
           oceanStressU, &
           oceanStressV, &
           surfaceTiltForceU, &
           surfaceTiltForceV, &
           uOceanVelocityVertex, &
           vOceanVelocityVertex

      ! geographical arrays
      real(kind=RKIND), dimension(:), pointer :: &
           uVelocityGeo, &
           vVelocityGeo, &
           stressDivergenceUGeo, &
           stressDivergenceVGeo, &
           airStressVertexUGeo, &
           airStressVertexVGeo, &
           oceanStressUGeo, &
           oceanStressVGeo, &
           surfaceTiltForceUGeo, &
           surfaceTiltForceVGeo, &
           uOceanVelocityVertexGeo, &
           vOceanVelocityVertexGeo

      integer, pointer :: &
           nVerticesSolve

      integer :: &
           iVertex

      err = 0

      dminfo = domain % dminfo

      ! get configs
      call MPAS_pool_get_config(domain % configs, "config_rotate_cartesian_grid", rotateCartesianGrid)

      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
         call MPAS_pool_get_subpool(block % structs, "geographicalVectorsAM", geographicalVectorsAMPool)

         ! mesh arrays
         call MPAS_pool_get_config(meshPool, "sphere_radius", sphereRadius)
         call MPAS_pool_get_array(meshPool, "latVertex", latVertex)
         call MPAS_pool_get_array(meshPool, "lonVertex", lonVertex)
         call MPAS_pool_get_array(meshPool, "xVertex", xVertex)
         call MPAS_pool_get_array(meshPool, "yVertex", yVertex)
         call MPAS_pool_get_array(meshPool, "zVertex", zVertex)

         call MPAS_pool_get_dimension(meshPool, "nVerticesSolve", nVerticesSolve)

         ! get non-geographical arrays
         call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
         call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
         call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
         call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
         call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
         call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)
         call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
         call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
         call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
         call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
         call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
         call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)

         ! get geographical arrays
         call MPAS_pool_get_array(geographicalVectorsAMPool, "uVelocityGeo", uVelocityGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "vVelocityGeo", vVelocityGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "stressDivergenceUGeo", stressDivergenceUGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "stressDivergenceVGeo", stressDivergenceVGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "airStressVertexUGeo", airStressVertexUGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "airStressVertexVGeo", airStressVertexVGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "oceanStressUGeo", oceanStressUGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "oceanStressVGeo", oceanStressVGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "surfaceTiltForceUGeo", surfaceTiltForceUGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "surfaceTiltForceVGeo", surfaceTiltForceVGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "uOceanVelocityVertexGeo", uOceanVelocityVertexGeo)
         call MPAS_pool_get_array(geographicalVectorsAMPool, "vOceanVelocityVertexGeo", vOceanVelocityVertexGeo)

         ! perform rotations
         do iVertex = 1, nVerticesSolve

            ! ice velocity
            call seaice_latlon_vector_rotation_backward(&
                 uVelocityGeo(iVertex), &
                 vVelocityGeo(iVertex), &
                 uVelocity(iVertex), &
                 vVelocity(iVertex), &
                 latVertex(iVertex), &
                 lonVertex(iVertex), &
                 xVertex(iVertex), &
                 yVertex(iVertex), &
                 zVertex(iVertex), &
                 sphereRadius, &
                 rotateCartesianGrid)

            ! stress divergence
            call seaice_latlon_vector_rotation_backward(&
                 stressDivergenceUGeo(iVertex), &
                 stressDivergenceVGeo(iVertex), &
                 stressDivergenceU(iVertex), &
                 stressDivergenceV(iVertex), &
                 latVertex(iVertex), &
                 lonVertex(iVertex), &
                 xVertex(iVertex), &
                 yVertex(iVertex), &
                 zVertex(iVertex), &
                 sphereRadius, &
                 rotateCartesianGrid)

            ! air stress
            call seaice_latlon_vector_rotation_backward(&
                 airStressVertexUGeo(iVertex), &
                 airStressVertexVGeo(iVertex), &
                 airStressVertexU(iVertex), &
                 airStressVertexV(iVertex), &
                 latVertex(iVertex), &
                 lonVertex(iVertex), &
                 xVertex(iVertex), &
                 yVertex(iVertex), &
                 zVertex(iVertex), &
                 sphereRadius, &
                 rotateCartesianGrid)

            ! ocean stress
            call seaice_latlon_vector_rotation_backward(&
                 oceanStressUGeo(iVertex), &
                 oceanStressVGeo(iVertex), &
                 oceanStressU(iVertex), &
                 oceanStressV(iVertex), &
                 latVertex(iVertex), &
                 lonVertex(iVertex), &
                 xVertex(iVertex), &
                 yVertex(iVertex), &
                 zVertex(iVertex), &
                 sphereRadius, &
                 rotateCartesianGrid)

            ! surface tilt
            call seaice_latlon_vector_rotation_backward(&
                 surfaceTiltForceUGeo(iVertex), &
                 surfaceTiltForceVGeo(iVertex), &
                 surfaceTiltForceU(iVertex), &
                 surfaceTiltForceV(iVertex), &
                 latVertex(iVertex), &
                 lonVertex(iVertex), &
                 xVertex(iVertex), &
                 yVertex(iVertex), &
                 zVertex(iVertex), &
                 sphereRadius, &
                 rotateCartesianGrid)

            ! ocean velocity
            call seaice_latlon_vector_rotation_backward(&
                 uOceanVelocityVertexGeo(iVertex), &
                 vOceanVelocityVertexGeo(iVertex), &
                 uOceanVelocityVertex(iVertex), &
                 vOceanVelocityVertex(iVertex), &
                 latVertex(iVertex), &
                 lonVertex(iVertex), &
                 xVertex(iVertex), &
                 yVertex(iVertex), &
                 zVertex(iVertex), &
                 sphereRadius, &
                 rotateCartesianGrid)

         enddo ! iVertex

         block => block % next
      end do

   end subroutine seaice_compute_geographical_vectors!}}}

!***********************************************************************
!
!  routine seaice_restart_geographical_vectors
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th October 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_geographical_vectors(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine seaice_restart_geographical_vectors!}}}

!***********************************************************************
!
!  routine seaice_finalize_geographical_vectors
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    14th October 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_geographical_vectors(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine seaice_finalize_geographical_vectors!}}}

!-----------------------------------------------------------------------

end module seaice_geographical_vectors

! vim: foldmethod=marker
