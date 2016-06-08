!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_ALLOCATE                                           C
!                                                                      C
!  Purpose: subroutines to allocate all DEM arrays                     C
!                                                                      C
!  Author: Rahul Garg                               Date: 1-Dec-2013   C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE DES_ALLOCATE

  PUBLIC:: DES_ALLOCATE_ARRAYS, ADD_PAIR, PARTICLE_GROW, PAIR_GROW, ALLOCATE_DEM_MI

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_ALLOCATE_ARRAYS                                     C
!  Purpose: Original allocate arrays subroutines for DES               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DES_ALLOCATE_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE constant
      USE discretelement
      Use indices
      Use geometry
      Use compar
      Use physprop
      Use des_bc
      Use pic_bc
      use funits
      USE mfix_pic
      Use des_thermo
      Use des_rxns
      USE cutcell
      USE functions

      use run, only: ENERGY_EQ
      use run, only: ANY_SPECIES_EQ

      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_GARG
      use particle_filter, only: DES_INTERP_DPVM
      use particle_filter, only: DES_INTERP_GAUSS
      use particle_filter, only: FILTER_SIZE
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_WEIGHT

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: IJK
!-----------------------------------------------

      CALL INIT_ERR_MSG("DES_ALLOCATE_ARRAYS")

! For parallel processing the array size required should be either
! specified by the user or could be determined from total particles
! with some factor.
      MAX_PIP = merge(0, PARTICLES/numPEs, PARTICLES==UNDEFINED_I)
      MAX_PIP = MAX(MAX_PIP,4)

      WRITE(ERR_MSG,1000) trim(iVal(MAX_PIP))
      CALL FLUSH_ERR_MSG(HEADER = .FALSE., FOOTER = .FALSE.)

 1000 FORMAT('Initial DES Particle array size: ',A)

! DES Allocatable arrays
!-----------------------------------------------
! Dynamic particle info including another index for parallel
! processing for ghost
      ALLOCATE( PEA (MAX_PIP, 4) )
      ALLOCATE (iglobal_id(max_pip))

! R.Garg: Allocate necessary arrays for PIC mass inlet/outlet BCs
      IF(PIC_BCMI /= 0 .OR. PIC_BCMO /=0) CALL ALLOCATE_PIC_MIO

! Particle attributes
! Radius, density, mass, moment of inertia
      Allocate(  DES_RADIUS (MAX_PIP) )
      Allocate(  RO_Sol (MAX_PIP) )
      Allocate(  PVOL (MAX_PIP) )
      Allocate(  PMASS (MAX_PIP) )
      Allocate(  OMOI (MAX_PIP) )

! Old and new particle positions, velocities (translational and
! rotational)
      Allocate(  DES_POS_NEW (DIMN,MAX_PIP) )
      Allocate(  DES_VEL_NEW (DIMN,MAX_PIP) )
      Allocate(  OMEGA_NEW (DIMN,MAX_PIP) )

      IF(PARTICLE_ORIENTATION) Allocate(  ORIENTATION (DIMN,MAX_PIP) )

      IF (DO_OLD) THEN
         Allocate(  DES_POS_OLD (DIMN,MAX_PIP) )
         Allocate(  DES_VEL_OLD (DIMN,MAX_PIP) )
         Allocate(  DES_ACC_OLD (DIMN,MAX_PIP) )
         Allocate(  OMEGA_OLD (DIMN,MAX_PIP) )
         Allocate(  ROT_ACC_OLD (DIMN,MAX_PIP))
      ENDIF

! Allocating user defined array
      IF(DES_USR_VAR_SIZE > 0) &
         Allocate( DES_USR_VAR(DES_USR_VAR_SIZE,MAX_PIP) )

! Particle positions at the last call neighbor search algorithm call
      Allocate(  PPOS (DIMN,MAX_PIP) )

! Total, normal and tangetial forces
      Allocate(  FC (DIMN,MAX_PIP) )

! Torque
      Allocate(  TOW (DIMN,MAX_PIP) )


! allocate variable for des grid binning
      allocate(dg_pijk(max_pip)); dg_pijk=0
      allocate(dg_pijkprv(max_pip)); dg_pijkprv=0

! allocate variables related to ghost particles
      allocate(ighost_updated(max_pip))



      Allocate(  wall_collision_facet_id (COLLISION_ARRAY_MAX, MAX_PIP) )
      wall_collision_facet_id(:,:) = -1
      Allocate(  wall_collision_PFT (DIMN, COLLISION_ARRAY_MAX, MAX_PIP) )

! Temporary variables to store wall position, velocity and normal vector
      Allocate(  WALL_NORMAL  (NWALLS,DIMN) )

      OLD_PAIR_NUM = 0
      PAIR_NUM = 0
      PAIR_MAX = 1024

      Allocate(  PAIRS (2,PAIR_MAX) )
      Allocate(  PAIRS_OLD (2,PAIR_MAX) )
      Allocate(  PV_PAIR (PAIR_MAX) )
      Allocate(  PV_PAIR_OLD (PAIR_MAX) )
      Allocate(  PFT_PAIR (3,PAIR_MAX) )
      Allocate(  PFT_PAIR_OLD (3,PAIR_MAX) )
      Allocate(  PFN_PAIR (3,PAIR_MAX) )
      Allocate(  PFN_PAIR_OLD (3,PAIR_MAX) )

      Allocate(  force_mag_PAIR (3,PAIR_MAX) )  !  1 normal, 2 tangential, 3 cohesive      
      
! Variable that stores the particle in cell information (ID) on the
! computational fluid grid defined by imax, jmax and kmax in mfix.dat
      ALLOCATE(PIC(DIMENSION_3))
      DO IJK=1,DIMENSION_3
        NULLIFY(pic(ijk)%p)
      ENDDO

! Particles in a computational fluid cell (for volume fraction)
      Allocate(  PINC (DIMENSION_3) )

! For each particle track its i,j,k location on computational fluid grid
! defined by imax, jmax and kmax in mfix.dat and phase no.
      Allocate(  PIJK (MAX_PIP,5) )

      ALLOCATE(DRAG_AM(DIMENSION_3))
      ALLOCATE(DRAG_BM(DIMENSION_3, DIMN))
      ALLOCATE(F_gp(MAX_PIP ))
      F_gp(1:MAX_PIP)  = ZERO

! Explict drag force acting on a particle.
      Allocate(DRAG_FC (DIMN,MAX_PIP) )

! force due to gas-pressure gradient
      ALLOCATE(P_FORCE(DIMN, DIMENSION_3))

! Volume averaged solids volume in a computational fluid cell
      Allocate(  DES_U_s (DIMENSION_3, DES_MMAX) )
      Allocate(  DES_V_s (DIMENSION_3, DES_MMAX) )
      Allocate(  DES_W_s (DIMENSION_3, DES_MMAX) )

! Volume of nodes
      ALLOCATE(DES_VOL_NODE(DIMENSION_3))

      ALLOCATE(F_GDS(DIMENSION_3))
      ALLOCATE(VXF_GDS(DIMENSION_3))

      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_DPVM, DES_INTERP_GAUSS)
         ALLOCATE(FILTER_CELL(FILTER_SIZE, MAX_PIP))
         ALLOCATE(FILTER_WEIGHT(FILTER_SIZE, MAX_PIP))
      CASE(DES_INTERP_GARG)
         ALLOCATE(DES_ROPS_NODE(DIMENSION_3, DES_MMAX))
         ALLOCATE(DES_VEL_NODE(DIMENSION_3, DIMN, DES_MMAX))
      END SELECT

! Variables for hybrid model
      IF (DES_CONTINUUM_HYBRID) THEN
         ALLOCATE(SDRAG_AM(DIMENSION_3,DIMENSION_M))
         ALLOCATE(SDRAG_BM(DIMENSION_3, DIMN,DIMENSION_M))

         ALLOCATE(F_SDS(DIMENSION_3,DIMENSION_M))
         ALLOCATE(VXF_SDS(DIMENSION_3,DIMENSION_M))
      ENDIF
! Bulk density in a computational fluid cell / for communication with
! MFIX continuum
      ALLOCATE( DES_ROP_S(DIMENSION_3, DES_MMAX) )
      ALLOCATE( DES_ROP_SO(DIMENSION_3, DES_MMAX) )

! MP-PIC related
      IF(MPPIC) THEN
         Allocate(PS_FORCE_PIC(DIMENSION_3, DIMN))
         ALLOCATE(DES_STAT_WT(MAX_PIP))
         ALLOCATE(DES_VEL_MAX(DIMN))
         ALLOCATE(PS_GRAD(MAX_PIP, DIMN))
         ALLOCATE(AVGSOLVEL_P(DIMN, MAX_PIP))
         ALLOCATE(EPG_P(MAX_PIP))

         Allocate(PIC_U_S(DIMENSION_3, DES_MMAX))
         Allocate(PIC_V_S(DIMENSION_3, DES_MMAX))
         Allocate(PIC_W_S(DIMENSION_3, DES_MMAX))

         Allocate(PIC_P_s (DIMENSION_3, DES_MMAX) )
!         ALLOCATE(MPPIC_VPTAU(MAX_PIP, DIMN))
         PIC_U_s = zero
         PIC_V_s = zero
         PIC_W_s = zero
         PIC_P_s = zero
      ENDIF


! Granular temperature in a computational fluid cell
      Allocate(DES_THETA (DIMENSION_3, DES_MMAX) )

! Averaged velocity obtained by averaging over all the particles
      ALLOCATE(DES_VEL_AVG(DIMN) )

! Global Granular Energy
      ALLOCATE(GLOBAL_GRAN_ENERGY(DIMN) )
      ALLOCATE(GLOBAL_GRAN_TEMP(DIMN) )

! variable for bed height of solids phase M
      ALLOCATE(BED_HEIGHT(DES_MMAX))

! ---------------------------------------------------------------->>>
! BEGIN COHESION
      IF(USE_COHESION) THEN
! Matrix location of particle  (should be allocated in case user wishes
! to invoke routines in /cohesion subdirectory
         Allocate(  PostCohesive (MAX_PIP) )
      ENDIF
! END COHESION
! ----------------------------------------------------------------<<<

! ---------------------------------------------------------------->>>
! BEGIN Thermodynamic Allocation
      IF(ENERGY_EQ)THEN
! Particle temperature
         Allocate( DES_T_s_OLD( MAX_PIP ) )
         Allocate( DES_T_s_NEW( MAX_PIP ) )
! Specific heat
         Allocate( DES_C_PS( MAX_PIP ) )
! Species mass fractions comprising a particle. This array may not be
! needed for all thermo problems.
         Allocate( DES_X_s( MAX_PIP, DIMENSION_N_S))
! Total rate of heat transfer to individual particles.
         Allocate( Q_Source( MAX_PIP ) )
! Average solids temperature in fluid cell
         Allocate(avgDES_T_s(DIMENSION_3) )

         Allocate(DES_ENERGY_SOURCE(DIMENSION_3) )

! Allocate the history variables for Adams-Bashforth integration
         IF (INTG_ADAMS_BASHFORTH) &
            Allocate( Q_Source0( MAX_PIP ) )
      ENDIF
! End Thermodynamic Allocation
! ----------------------------------------------------------------<<<


! ---------------------------------------------------------------->>>
! BEGIN Species Allocation
      IF(ANY_SPECIES_EQ)THEN
! Rate of solids phase production for each species
         Allocate( DES_R_sp( MAX_PIP, DIMENSION_N_s) )
! Rate of solids phase consumption for each species
         Allocate( DES_R_sc( MAX_PIP, DIMENSION_N_s) )


         Allocate( DES_R_gp( DIMENSION_3, DIMENSION_N_g ) )
         Allocate( DES_R_gc( DIMENSION_3, DIMENSION_N_g ) )
         Allocate( DES_SUM_R_g( DIMENSION_3 ) )
         Allocate( DES_R_PHASE( DIMENSION_3, DIMENSION_LM+DIMENSION_M-1 ) )
         Allocate( DES_HOR_g( DIMENSION_3 ) )


! Allocate the history variables for Adams-Bashforth integration
         IF (INTG_ADAMS_BASHFORTH) THEN
! Rate of change of particle mass
            Allocate( dMdt_OLD( MAX_PIP ) )
! Rate of change of particle mass percent species
            Allocate( dXdt_OLD( MAX_PIP, DIMENSION_N_s) )
         ENDIF

! Energy generation from reaction (cal/sec)
         Allocate( Qint( MAX_PIP ) )
      ENDIF
! End Species Allocation
! ----------------------------------------------------------------<<<

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE DES_ALLOCATE_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_DEM_MIO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE ALLOCATE_DEM_MI

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE des_bc
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------

! Particle injection factor
      Allocate( PI_FACTOR (DEM_BCMI) )
! Particle injection count (injection number)
      Allocate( PI_COUNT (DEM_BCMI) )
! Particle injection time scale
      Allocate( DEM_MI_TIME (DEM_BCMI) )
! Array used for polydisperse inlets: stores the particle number
! distribution of an inlet scaled with numfrac_limit
      Allocate( DEM_BC_POLY_LAYOUT( DEM_BCMI, NUMFRAC_LIMIT ) )
! Data structure for storing BC data.
      Allocate( DEM_MI(DEM_BCMI) )

! Initializiation
! Integer arrays
      PI_FACTOR(:) = -1
      PI_COUNT(:) = -1
      DEM_BC_POLY_LAYOUT(:,:) = -1
! Double precision arrays
      DEM_MI_TIME(:) = UNDEFINED

      allocate( DEM_BCMI_IJKSTART(DEM_BCMI) )
      allocate( DEM_BCMI_IJKEND(DEM_BCMI) )

      DEM_BCMI_IJKSTART = -1
      DEM_BCMI_IJKEND   = -1

! Boundary classification
!         Allocate( PARTICLE_PLCMNT (DES_BCMI) )
! Character precision arrays
!         PARTICLE_PLCMNT(:) = UNDEFINED_C

      RETURN
      END SUBROUTINE ALLOCATE_DEM_MI

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_PIC_MIO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: R. Garg                                    Date: 11-Jun-14  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE ALLOCATE_PIC_MIO

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE pic_bc
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------

! Allocate/Initialize for inlets
      IF(PIC_BCMI /= 0)THEN

         allocate( PIC_BCMI_IJKSTART(PIC_BCMI) )
         allocate( PIC_BCMI_IJKEND  (PIC_BCMI) )
         allocate( PIC_BCMI_NORMDIR (PIC_BCMI,3) )

         ALLOCATE( PIC_BCMI_OFFSET  (PIC_BCMI,3))

         ALLOCATE( PIC_BCMI_INCL_CUTCELL(PIC_BCMI) )

         PIC_BCMI_IJKSTART = -1
         PIC_BCMI_IJKEND   = -1

      ENDIF  ! end if PIC_BCMI /= 0



      IF(PIC_BCMO > 0)THEN
         allocate( PIC_BCMO_IJKSTART(PIC_BCMO) )
         allocate( PIC_BCMO_IJKEND(PIC_BCMO) )

         PIC_BCMO_IJKSTART = -1
         PIC_BCMO_IJKEND   = -1
      ENDIF


      RETURN
      END SUBROUTINE ALLOCATE_PIC_MIO


!``````````````````````````````````````````````````````````````````````!
! Subroutine: ADD_PAIR                                                 !
!                                                                      !
! Purpose: Adds a neighbor pair to the pairs array.                    !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE add_pair(ii,jj)
      USE discretelement
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ii,jj

      pair_num = pair_num +1

! Reallocate to double the size of the arrays if needed.
      IF(PAIR_NUM > PAIR_MAX) THEN
         PAIR_MAX = PAIR_MAX*2
         CALL PAIR_GROW
      ENDIF

      pairs(1,pair_num) = ii
      pairs(2,pair_num) = jj

      RETURN
      END SUBROUTINE add_pair

!``````````````````````````````````````````````````````````````````````!
! Subroutine: PAIR_GROW                                                !
!                                                                      !
! Purpose: Grow pair arrays to pair_max. Note that pair                !
! max should be increased before calling this routine. Also, no        !
! assumption to the previous array size is made as needed for restarts.!
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE PAIR_GROW

        USE discretelement

        IMPLICIT NONE

        call integer_grow2(pairs,pair_max)
        call integer_grow2(pairs_old,pair_max)

        call logical_grow(pv_pair,pair_max)
        call logical_grow(pv_pair_old,pair_max)

        call real_grow2(pft_pair,pair_max)
        call real_grow2(pft_pair_old,pair_max)
        call real_grow2(pfn_pair,pair_max)
        call real_grow2(pfn_pair_old,pair_max)
        
        call real_grow2(force_mag_pair,pair_max)

        RETURN

      END SUBROUTINE PAIR_GROW

!``````````````````````````````````````````````````````````````````````!
! Subroutine: PARTICLE_GROW                                            !
!                                                                      !
! Purpose: Grow pair arrays to pair_max. Note that pair                !
! max should be increased before calling this routine. Also, no        !
! assumption to the previous array size is made as needed for restarts.!
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE PARTICLE_GROW(new_max_pip)

!      use discretelement, only: pea, wall_collision_facet_id, wall_collision_pft, MAX_PIP, PIP
!      use param1

        USE des_rxns
        USE des_thermo
        USE mfix_pic
        USE discretelement
        USE particle_filter
        USE run

        IMPLICIT NONE

        integer, intent(in) :: new_max_pip

        DO WHILE (MAX_PIP < new_max_pip)
           MAX_PIP = MAX_PIP*2

           call real_grow(des_radius,MAX_PIP)
           call real_grow(RO_Sol,MAX_PIP)
           call real_grow(PVOL,MAX_PIP)
           call real_grow(PMASS,MAX_PIP)
           call real_grow(OMOI,MAX_PIP)
           call real_grow2(DES_POS_NEW,MAX_PIP)
           call real_grow2(DES_VEL_NEW,MAX_PIP)
           call real_grow2(OMEGA_NEW,MAX_PIP)
           IF(PARTICLE_ORIENTATION) call real_grow2(ORIENTATION,MAX_PIP)
           call real_grow2(PPOS,MAX_PIP)
           call logical_grow2_reverse(PEA,MAX_PIP)
           call integer_grow(iglobal_id,MAX_PIP)
           call integer_grow2_reverse(pijk,MAX_PIP)
           call integer_grow(dg_pijk,MAX_PIP)
           call integer_grow(dg_pijkprv,MAX_PIP)
           call logical_grow(ighost_updated,MAX_PIP)
           call real_grow2(FC,MAX_PIP)
           call real_grow2(TOW,MAX_PIP)
           call real_grow(F_GP,MAX_PIP)
           call integer_grow2(WALL_COLLISION_FACET_ID,MAX_PIP)
           call real_grow3(WALL_COLLISION_PFT,MAX_PIP)
           call real_grow2(DRAG_FC,MAX_PIP)

           IF(FILTER_SIZE > 0) THEN
              call integer_grow2(FILTER_CELL,MAX_PIP)
              call real_grow2(FILTER_WEIGHT,MAX_PIP)
           ENDIF

           IF(MPPIC) THEN
              call real_grow(DES_STAT_WT,MAX_PIP)
              call real_grow2_reverse(PS_GRAD,MAX_PIP)
              call real_grow2(AVGSOLVEL_P,MAX_PIP)
              call real_grow(EPG_P,MAX_PIP)
           ENDIF

           IF(USE_COHESION) THEN
              call real_grow(PostCohesive,MAX_PIP)
           ENDIF

           IF (DO_OLD) THEN
              call real_grow2(DES_POS_OLD,MAX_PIP)
              call real_grow2(DES_VEL_OLD,MAX_PIP)
              call real_grow2(DES_ACC_OLD,MAX_PIP)
              call real_grow2(OMEGA_OLD,MAX_PIP)
              call real_grow2(ROT_ACC_OLD,MAX_PIP)
           ENDIF

           IF(ENERGY_EQ)THEN
              call real_grow(DES_T_s_OLD,MAX_PIP)
              call real_grow(DES_T_s_NEW,MAX_PIP)
              call real_grow(DES_C_PS,MAX_PIP)
              call real_grow2_reverse(DES_X_s,MAX_PIP)
              call real_grow(Q_Source,MAX_PIP)

              IF (INTG_ADAMS_BASHFORTH) &
                   call real_grow(Q_Source0,MAX_PIP)
           ENDIF

           IF(ANY_SPECIES_EQ)THEN
              call real_grow2_reverse( DES_R_sp, MAX_PIP )
              call real_grow2_reverse( DES_R_sc, MAX_PIP )

              IF (INTG_ADAMS_BASHFORTH) THEN
                 call real_grow( dMdt_OLD, MAX_PIP )
                 call real_grow2_reverse( dXdt_OLD, MAX_PIP )
              ENDIF

              call real_grow( Qint, MAX_PIP )
           ENDIF

           IF(DES_USR_VAR_SIZE > 0) &
              call real_grow2(DES_USR_VAR,MAX_PIP)

           CALL DES_INIT_PARTICLE_ARRAYS(MAX_PIP/2+1,MAX_PIP)

        ENDDO

      RETURN

      END SUBROUTINE PARTICLE_GROW

      SUBROUTINE INTEGER_GROW(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE

        lSIZE = size(integer_array,1)
        allocate(integer_tmp(new_size))
        integer_tmp(1:lSIZE) = integer_array(1:lSIZE)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW

      SUBROUTINE INTEGER_GROW2_reverse(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(integer_array,1)
        lSIZE2 = size(integer_array,2)
        allocate(integer_tmp(new_size,lSIZE2))
        integer_tmp(1:lSIZE,:) = integer_array(1:lSIZE,:)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW2_reverse

      SUBROUTINE INTEGER_GROW2(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(integer_array,1)
        lSIZE2 = size(integer_array,2)
        allocate(integer_tmp(lSIZE,new_size))
        integer_tmp(:,1:lSIZE2) = integer_array(:,1:lSIZE2)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW2

      SUBROUTINE LOGICAL_GROW(logical_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: logical_array
        LOGICAL, DIMENSION(:), ALLOCATABLE :: logical_tmp
        INTEGER lSIZE2

        lSIZE2 = size(logical_array,1)
        allocate(logical_tmp(new_size))
        logical_tmp(1:lSIZE2) = logical_array(1:lSIZE2)
        call move_alloc(logical_tmp,logical_array)

      END SUBROUTINE LOGICAL_GROW

      SUBROUTINE REAL_GROW(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE

        lSIZE = size(real_array,1)
        allocate(real_tmp(new_size))
        real_tmp(1:lSIZE) = real_array(1:lSIZE)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW

      SUBROUTINE REAL_GROW2(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        allocate(real_tmp(lSIZE,new_size))
        real_tmp(:,1:lSIZE2) = real_array(:,1:lSIZE2)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW2

      SUBROUTINE REAL_GROW2_reverse(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        allocate(real_tmp(new_size,lSIZE2))
        real_tmp(1:lSIZE,:) = real_array(1:lSIZE,:)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW2_REVERSE

      SUBROUTINE REAL_GROW3(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2, lSIZE3

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        lSIZE3 = size(real_array,3)
        allocate(real_tmp(lSIZE,lSIZE2,new_size))
        real_tmp(:,:,1:lSIZE3) = real_array(:,:,1:lSIZE3)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW3

      SUBROUTINE LOGICAL_GROW2_REVERSE(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        allocate(real_tmp(new_size,lSIZE2))
        real_tmp(1:lSIZE,:) = real_array(1:lSIZE,:)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE LOGICAL_GROW2_REVERSE

    END MODULE DES_ALLOCATE
