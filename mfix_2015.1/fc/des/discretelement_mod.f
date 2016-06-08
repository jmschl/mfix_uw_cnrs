!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!   Module name: DISCRETELEMENT                                        C
!   Purpose: DES mod file                                              C
!            Common Block containing DEM conditions                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE DISCRETELEMENT

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------

! Define interface - needed when passing arrays of assumed size
! This function is used to identify fluid grid index (i,j or k) of a
! particle in a given x-, y-, z- line of space
      INTERFACE
         INTEGER FUNCTION DES_GETINDEXFROMPOS(LIM1,LIM2,PART_POS,&
            GRID_POS,AXIS,AXIS_INDEX)
            INTEGER :: LIM1, LIM2
            DOUBLE PRECISION :: PART_POS
            DOUBLE PRECISION, DIMENSION(:) :: GRID_POS
            CHARACTER(LEN=1) :: AXIS,AXIS_INDEX
         END FUNCTION DES_GETINDEXFROMPOS
      END INTERFACE


! Total number of particles in simulation: read from input or generated
      INTEGER PARTICLES

! Start particle tracking quantities
!----------------------------------------------------------------->>>
! Generally for inlet/outlet related routines but also employed in
! tracking for parallelization

! Dynamic particle count elements:
! PEA(n,1) : This column identifies particle as 'existing' if true.
! It is used with the inlet/outlet to skip indices that do not represent
! particles in the system or indices that represent particles that have
! exited the system.

! PEA(n,2) : This column identifies a particle as 'new' if true.
! Particles with a classification of 'new' do not react when in contact
! with a wall or another particle, however existing particles do collide
! and interact with 'new' particles. The classification allows new
! particles to push particles already in the system out of the way when
! entering to prevent overlap.  This flag is also used when the center
! of a particle crosses a dem outlet (i.e. an exiting particle; see
! PEA(n,3)) so that the particle will maintain its present trajectory
! until it has fully exited the system

! PEA(n,3) : This column identifies a particle as 'exiting' if true.
! If a particle initiates contact with a wall surface designated as a
! des outlet, this flag is set to true. With this classification the
! location of the particle is checked to assess if the particle has
! fully exited the system.  At this point, the particle is removed
! from the list.

! PEA(n,4) : for ghost particles
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: PEA ! (PARTICLES,4)

! PARALLEL PROCESSING: explanation of variables in parallel architecture
! pip - particles in each processor (includes the ghost particles)
! max_pip - maximum allocated particles in processor

! Number of particles in the system (current)
      INTEGER PIP
! Global sum of particles (excluding ghost) in the system
      INTEGER TOT_PAR
! Maximum particles permitted in the system at once
      INTEGER MAX_PIP

! End particle tracking quantities
!-----------------------------------------------------------------<<<

! For parallel processing: global id of particles
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGLOBAL_ID
! Ghost count of particles on each processor
      INTEGER :: IGHOST_CNT
! Maximum global id, new particles global id will be assigned based
! on this value
      Integer :: imax_global_id


! If gener_part_config is true, then the particle_input.dat file does
! not need to be supplied nor does the total number of particles as
! these are determined based on the specified volume fraction (vol_frac)
! in the specified domain (des_eps_xyzstart)
      LOGICAL :: GENER_PART_CONFIG
      DOUBLE PRECISION ::  VOL_FRAC(DIM_M)
      DOUBLE PRECISION :: DES_EPS_XSTART, &
      DES_EPS_YSTART, DES_EPS_ZSTART
! volume of the IC region for computing number of particles to be seeded
      DOUBLE PRECISION, dimension(:), allocatable :: VOL_IC_REGION!(DIMENSION_IC)
! number of particles for each phase corresponding to the IC number. This will
! be real particles for DEM but parcels or computational particles for PIC model
      INTEGER, dimension(:,:), allocatable :: PART_MPHASE_BYIC!(DIMENSION_IC, DIM_M)

! Number of real particles by IC and by solid phase. Only relevant for PIC model
      double precision, dimension(:,:), allocatable :: REALPART_MPHASE_BYIC!(DIMENSION_IC, DIM_M)
! The number of particles that belong to solid phase M according to the
! vol_frac and particle diameter. this information is used when
! gener_part_config is invoked for initialization
! This will be removed soon as PART_MPHASE_BYIC will be used from now on
      INTEGER PART_MPHASE(DIM_M)

! The number of real particles that belong to solid phase M during the initialization.
! It is equal to Part_mphase for DEM but implies real number of particles for PIC model
      double precision REALPART_MPHASE(DIM_M)
! Assigns the initial particle velocity distribution based on user
! specified mean and standard deviation (regardless if already set
! within particle_input.dat)
      DOUBLE PRECISION pvel_mean, PVEL_StDev

! Output/debug controls
!----------------------------------------------------------------->>>
! Logic that controls whether to print data dem simulations (granular or
! coupled)
      LOGICAL PRINT_DES_DATA

! logic that controls if des run time messages are printed on screen or not
      LOGICAL PRINT_DES_SCREEN

! Usr specified time interval that controls frequency of writing DEM
! output and restart for pure granular flow; otherwise (when coupled)
! the frequency of writing output and restart is controlled by the
! value of spx_dt(1) and res_dt
      DOUBLE PRECISION DES_SPX_DT, DES_RES_DT

! This specifies the file type used for outputting DES data
! options are :
!    TECPLOT - data is written in Tecplot format
!    undefined - data is written in ParaView format (default)
      CHARACTER(LEN=64) :: DES_OUTPUT_TYPE

! Used sporadically to control screen dumps (for debug purposes)
      LOGICAL :: DEBUG_DES

! Single particle no. index that is followed if debugging
      INTEGER FOCUS_PARTICLE

! Output file count for .vtp type files and for tecplot files;
! for vtp output used to label .vtp file names in sequential order
! and is saved so restarts begin at the correct count
      INTEGER VTP_FINDEX, TECPLOT_FINDEX
! End Output/debug controls
!-----------------------------------------------------------------<<<



! DES - Continuum
      LOGICAL DISCRETE_ELEMENT
      LOGICAL DES_CONTINUUM_COUPLED

! DES - Invoke hybrid model where both the DEM and continuum model
! are employed to describe solids
      LOGICAL DES_CONTINUUM_HYBRID

! DES -
! With this logic the particles see the fluid but the fluid does
! not see the particles.
      LOGICAL DES_ONEWAY_COUPLED

! Only used when coupled and represents the number of times a pure
! granular flow simulation is run before the actual coupled simulation
! is started (i.e. for particle settling w/o fluid forces)
      INTEGER NFACTOR

! Drag
      LOGICAL TSUJI_DRAG

! Collision model, options are as follows
!   linear spring dashpot model (default/undefined)
!   'hertzian' model
      CHARACTER(LEN=64) :: DES_COLL_MODEL
      INTEGER DES_COLL_MODEL_ENUM
      INTEGER,PARAMETER ::  HERTZIAN=0
      INTEGER,PARAMETER ::  LSD=1

! Integration method, options are as follows
!   'euler' first-order scheme (default)
!   'adams_bashforth' second-order scheme (by T.Li)
      CHARACTER(LEN=64) :: DES_INTG_METHOD
      LOGICAL INTG_ADAMS_BASHFORTH
      LOGICAL INTG_EULER

! Value of solids time step based on particle properties
      DOUBLE PRECISION DTSOLID
! Run time value of simulation time used in dem simulation
      DOUBLE PRECISION S_TIME


! Neighbor search related quantities
!----------------------------------------------------------------->>>
! Neighbor search method, options are as follows
!   1= nsquare, 2=quadtree, 3=octree, 4=grid/cell based search
      INTEGER DES_NEIGHBOR_SEARCH

! Quantities used to determine whether neighbor search should be called
      INTEGER NEIGHBOR_SEARCH_N
      DOUBLE PRECISION NEIGHBOR_SEARCH_RAD_RATIO
      LOGICAL DO_NSEARCH

! Flag on whether to have DES_*_OLD arrays, if either Adams Bashforth or PIC is used
      LOGICAL DO_OLD

! Factor muliplied by sum of radii in grid based neighbor search and
! nsquare search method.  increases the effective radius of a particle
! for detecting particle contacts
      DOUBLE PRECISION FACTOR_RLM

! Stores number of neighbors based on neighbor search
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PAIRS
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PAIRS_OLD
      LOGICAL, DIMENSION(:), ALLOCATABLE :: PV_PAIR
      LOGICAL, DIMENSION(:), ALLOCATABLE :: PV_PAIR_OLD
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PFT_PAIR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PFT_PAIR_OLD
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PFN_PAIR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PFN_PAIR_OLD
      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: force_mag_PAIR      

      INTEGER, DIMENSION(:), ALLOCATABLE :: CELLNEIGHBOR_FACET_NUM, CELLNEIGHBOR_FACET_MAX
      INTEGER :: PAIR_NUM,OLD_PAIR_NUM,PAIR_MAX

! Quantities used for reporting: max no. neighbors and max overlap
! that exists during last solid time step of dem simulation
      DOUBLE PRECISION OVERLAP_MAX

! The number of i, j, k divisions in the grid used to perform the
! cell based neighbor search
      INTEGER :: DESGRIDSEARCH_IMAX, DESGRIDSEARCH_JMAX, &
                 DESGRIDSEARCH_KMAX

! End neighbor search related quantities
!-----------------------------------------------------------------<<<


! User specified dimension of the system (by default 2D, but if 3D system is
! desired then it must be explicitly specified)
      INTEGER, PARAMETER :: DIMN = 3

! Variable that is set to the number of walls in the system
      INTEGER NWALLS

! Position of domain boundaries generally given as
!   (0, xlength, 0, ylength, 0, zlength)
      DOUBLE PRECISION WX1, EX2, BY1, TY2, SZ1, NZ2

! X, Y, Z position of cell faces of computational fluid grid
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XE  !(0:DIMENSION_I)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: YN  !(0:DIMENSION_J)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZT  !(0:DIMENSION_K)

! Wall normal vector
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WALL_NORMAL  !(NWALLS,3)

! Gravity vector and magnitude
      DOUBLE PRECISION :: GRAV(3)
      DOUBLE PRECISION :: GRAV_MAG


! Periodic wall BC
      LOGICAL DES_PERIODIC_WALLS
      LOGICAL DES_PERIODIC_WALLS_X
      LOGICAL DES_PERIODIC_WALLS_Y
      LOGICAL DES_PERIODIC_WALLS_Z


! Lees & Edwards wall BC (lost in current DEM)
!----------------------------------------------------------------->>>
! Logic for Lees & Edwards BC (T = turn on LE BC)
      LOGICAL DES_LE_BC
! Relative velocity of LE boundaries (distance/time)
      DOUBLE PRECISION DES_LE_REL_VEL
! Shear direction
!   2D options are DUDY or DVDX
!   3D options are DUDY, DUDZ, DVDX, DVDZ, DWDX or DWDY
!   Note that all other directions are treated as periodic boundaries
      CHARACTER(LEN=4) :: DES_LE_SHEAR_DIR
! End LE BC
!-----------------------------------------------------------------<<<



! Particle-particle and Particle-wall collision model parameters
!----------------------------------------------------------------->>>
! Spring contants
      DOUBLE PRECISION KN, KN_W  !Normal
      DOUBLE PRECISION KT, KT_W, KT_FAC, KT_W_FAC  ! Tangential factors = KT/KN and KT_w/KN_w, resp.

! Damping coeffients
      DOUBLE PRECISION ETA_DES_N, ETA_N_W  !Normal
      DOUBLE PRECISION ETA_DES_T, ETA_T_W  !Tangential

! Flag to use van der Hoef et al. (2006) model for adjusting the rotation of the
! contact plane
      LOGICAL :: USE_VDH_DEM_MODEL
! Tangential damping factors, eta_t = eta_t_factor * eta_N
      DOUBLE PRECISION DES_ETAT_FAC, DES_ETAT_W_FAC

! Damping coeffients in array form
      DOUBLE PRECISION :: DES_ETAN(DIM_M, DIM_M), DES_ETAN_WALL(DIM_M)
      DOUBLE PRECISION :: DES_ETAT(DIM_M, DIM_M), DES_ETAT_WALL(DIM_M)

! Friction coeficients
      DOUBLE PRECISION MEW, MEW_W

! coeff of restituion input in one D array, solid solid
! Tangential rest. coef. are used for hertzian collision model but not linear
      DOUBLE PRECISION DES_EN_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)
      DOUBLE PRECISION DES_ET_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)

! coeff of restitution input in one D array, solid wall
      DOUBLE PRECISION DES_EN_WALL_INPUT(DIM_M)
      DOUBLE PRECISION DES_ET_WALL_INPUT(DIM_M)

! Hertzian collision model:
      DOUBLE PRECISION :: E_YOUNG(DIM_M), Ew_YOUNG
      DOUBLE PRECISION :: V_POISSON(DIM_M), Vw_POISSON
      DOUBLE PRECISION :: HERT_KN(DIM_M, DIM_M), HERT_KWN(DIM_M)
      DOUBLE PRECISION :: HERT_KT(DIM_M, DIM_M), HERT_KWT(DIM_M)
      DOUBLE PRECISION :: G_MOD(DIM_M)

! End particle-particle and particle-wall collision model parameters
!-----------------------------------------------------------------<<<


! Particle attributes: radius, density, mass, moment of inertia
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_RADIUS !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RO_Sol     !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PVOL       !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PMASS      !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OMOI       !(PARTICLES)

! Additional quantities
      DOUBLE PRECISION :: MIN_RADIUS, MAX_RADIUS


! 'solids phase' particle diameters
      DOUBLE PRECISION DES_D_P0 (DIM_M)
! 'solids phase' particle densities
      DOUBLE PRECISION DES_RO_s (DIM_M)
! number of 'solids phases'
      INTEGER DES_MMAX

! Old and new particle positions, velocities (translational and
! rotational)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_OLD  !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_NEW  !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_OLD  !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_NEW  !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_OLD    !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_NEW    !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PPOS         !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_ACC_OLD  !(PARTICLES,3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ROT_ACC_OLD  !(PARTICLES,3)

! Particle orientation
      LOGICAL :: PARTICLE_ORIENTATION = .FALSE.
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ORIENTATION  !(3,PARTICLES)
      DOUBLE PRECISION, DIMENSION(3) :: INIT_ORIENTATION = (/0.0,1.0,0.0/)

! Defining user defined allocatable array
      INTEGER :: DES_USR_VAR_SIZE
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_USR_VAR  !(PARTICLES,3)

! Total force and torque on each particle
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FC    !(3,PARTICLES)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TOW   !(3,PARTICLES)

! Dynamic information related to computational (eulerian) fluid grid
!----------------------------------------------------------------->>>
! Dynamic variable. for each ijk computational fluid cell store the
! total number of particles and the id's of the particles in that cell
      TYPE iap1
         INTEGER, DIMENSION(:), POINTER:: p
      END TYPE iap1

!     particle can collide with at most COLLISION_ARRAY_MAX facets simultaneously
      INTEGER, PARAMETER :: COLLISION_ARRAY_MAX = 8

!     -1 value indicates no collision
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: wall_collision_facet_id       ! (COLLISION_ARRAY_MAX,PARTICLES)
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: wall_collision_PFT ! (DIMN,COLLISION_ARRAY_MAX,PARTICLES)

      TYPE cnaa1
         INTEGER, DIMENSION(:), ALLOCATABLE:: p
         INTEGER, DIMENSION(:), ALLOCATABLE:: extentdir
         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: extentmin
         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: extentmax
      END TYPE cnaa1
! in order to facilitate the parallel processing the PIC is defined
! as single array IJK
      TYPE(iap1), DIMENSION(:), ALLOCATABLE:: pic  ! (DIMENSION_3)

      TYPE(cnaa1), DIMENSION(:), ALLOCATABLE :: CELLNEIGHBOR_FACET

! Store the number of particles in a computational fluid cell
      INTEGER, DIMENSION(:), ALLOCATABLE :: PINC  ! (DIMENSION_3)

! For each particle track its i, j, k & ijk location on the fluid grid
! and solids phase no.:
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PIJK ! (PARTICLES,5)=>I,J,K,IJK,M
!-----------------------------------------------------------------<<<


! note that thse variables are needed since the existing variables (i.e.
! f_gs, f_ss, etc) are also being used to store the information between
! the gas and continuous solids phases.
! drag coefficient between gas phase and discrete particle 'phases'
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: F_GDS
! drag coefficient between continuous solids phases and discrete
! particle 'phases'
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F_SDS

! the following should probably be local to the subroutine
! solve_vel_star they are only needed when invoking the non-interpolated
! version of drag wherein they are used to determine part of the source
! in the gas-phase momentum balances, in the routine to determine
! coefficients for the pressure correction equation and in the partial
! elimination algorithm
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VXF_GDS
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: VXF_SDS

! the contribution of solids-particle drag to the mth phase continuum
! solids momentum A matrix (center coefficient)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: SDRAG_AM
! the contribution of solids-particle drag to the to mth phase continuum
! solids momentum B vector
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: SDRAG_BM


! the coefficient add to gas momentum A matrix
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DRAG_AM
! the coefficient add to gas momentum B matrix
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DRAG_BM

! Explictly calculated fluid-particle drag force.
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DRAG_FC

! An intermediate array used in calculation of mean solids velocity
! by backward interpolation, i.e., when INTERP_DES_MEAN_FIELDS is true.
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::DES_VEL_NODE
                        !(DIMENSION_3,3,DES_MMAX)

! An intermediate array used in calculation of solids volume fraction
! by backward interpolation, i.e., when INTERP_DES_MEAN_FIELDS is true.
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  DES_ROPS_NODE
                        !(DIMENSION_3,DES_MMAX)

      DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: weightp

! the gas-particle drag coefficient
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: f_gp
                        !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: wtderivp
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: sstencil
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: gstencil, vstencil, pgradstencil

! stencil for interpolation of solids pressure
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE ::  psgradstencil

! stencil for interpolation of solids velocity
      DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE::  VEL_SOL_STENCIL


! quantities are set in subroutine set_interpolation_scheme
! order = order of the interpolation method, ob2l = (order+1)/2,
! ob2r = order/2
      CHARACTER(LEN=7):: scheme, interp_scheme
      INTEGER:: order, ob2l, ob2r

! END of interpolation related data
!-----------------------------------------------------------------<<<

      ! Volume of each node. Used to obtain Eulerian fields
      double precision, allocatable, dimension(:) :: des_vol_node
                        !(DIMENSION_3)

! Ratio of actual volume of each node to volume of node not corrected for
! outside the domain or being in cut-cell
      double precision, allocatable, dimension(:) :: des_vol_node_ratio
                        !(DIMENSION_3)

! Variable to track pressure force in computational fluid cell (ijk)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_FORCE
                        !(DIMENSION_3,DES_MMAX)

! Bulk density of particles in fluid cell
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_ROP_S,&
                                                       DES_ROP_SO
                        !(DIMENSION_3,DES_MMAX)

! Volume averaged solids velocity in a fluid cell
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_U_s
                        !(DIMENSION_3,DES_MMAX)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_V_s
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_W_s

! Granular temperature in a fluid cell
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_THETA
                        !(DIMENSION_3,DES_MMAX)

! Global average velocity: obtained by averaging over all the particles
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_VEL_AVG
                        !(3)

! Global granular energy & temp: obtained by averaging over all the particles
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GLOBAL_GRAN_ENERGY
                        !(3)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GLOBAL_GRAN_TEMP
                        !(3)

! Kinetic and potential energy of the system: obtained by averaging
! over all particles
! Added rotational kinetic energy (DES_ROTE)
      DOUBLE PRECISION DES_KE, DES_PE, DES_ROTE

! Logic for bed height calculations (T = turn on bed height
! calculations)
      LOGICAL DES_CALC_BEDHEIGHT
! Used to track bed height of solids phase M
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bed_height
                       !(DES_MMAX)

! MAX velocity of particles in each direction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_VEL_MAX
                        !(3)

!flag to convert the outside box to facets if stl facet representation is used
!default is false

      Logical :: des_convert_box_to_facets

      Integer, parameter :: FACET_TYPE_NORMAL = 1, FACET_TYPE_PO = 2 &
      , FACET_TYPE_MI = 3

      !make this a short integer
      !array to specify the facet type
      Integer, dimension(:), allocatable :: STL_FACET_TYPE

      Integer :: count_facet_type_normal, count_facet_type_po, count_facet_type_mi

      Logical :: minimize_des_facet_list


!-----------------------------------------------------------------<<<


! Start Cohesion
!----------------------------------------------------------------->>>
! Includes square-well type model and a van der waals type model

! Switch to turn cohesion on and off (set in mfix.dat)
      LOGICAL USE_COHESION


! Van der Waals constants (set in mfix.dat)
      LOGICAL SQUARE_WELL, VAN_DER_WAALS
      DOUBLE PRECISION HAMAKER_CONSTANT
      DOUBLE PRECISION VDW_INNER_CUTOFF ! (in cm)
      DOUBLE PRECISION VDW_OUTER_CUTOFF
      DOUBLE PRECISION WALL_HAMAKER_CONSTANT
      DOUBLE PRECISION WALL_VDW_INNER_CUTOFF
      DOUBLE PRECISION WALL_VDW_OUTER_CUTOFF
      DOUBLE PRECISION Asperities ! average radius of asperities (default zero)

! Store postcohesive
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PostCohesive
                        !(PARTICLES)
! Store cluster information array for postprocessing
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PostCluster

! Variables for van der waals cohesion calculations:
! Surface energy used to calculate cohesive force for low separation distances
! in Van der Waals model (this variable is calculated at the beginning of each
! simulation to ensure the van der Waals force is continuous at the inner cutoff)
      DOUBLE PRECISION SURFACE_ENERGY
      DOUBLE PRECISION WALL_SURFACE_ENERGY

! END Cohesion
!-----------------------------------------------------------------<<<

      LOGICAL :: DES_EXPLICITLY_COUPLED

! particle in cell related variable
      type iap2
         integer :: isize
         integer, dimension(:), pointer:: p
      end type iap2

      type(iap2), dimension(:),allocatable:: dg_pic
      integer, dimension(:),allocatable :: dg_pijk,dg_pijkprv

! variable to clean the ghost cells
      logical,dimension(:),allocatable :: ighost_updated

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: DES_CROSSPRDCT
!  Purpose: Calculate the cross product of two vectors that both have
!           either 2 or 3 elements and return the result in the first
!           argument
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      FUNCTION DES_CROSSPRDCT (XX,YY)

!-----------------------------------------------
! Modules
!-----------------------------------------------
        USE param
        USE param1
        use geometry, only: DO_K
        IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! sent vectors
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: XX, YY
! returned result: cross product of vectors
      DOUBLE PRECISION, DIMENSION(3) :: DES_CROSSPRDCT
!-----------------------------------------------

      IF(DO_K) THEN
         DES_CROSSPRDCT(1) = XX(2)*YY(3) - XX(3)*YY(2)
         DES_CROSSPRDCT(2) = XX(3)*YY(1) - XX(1)*YY(3)
      ELSE
         DES_CROSSPRDCT(1) = ZERO
         DES_CROSSPRDCT(2) = ZERO
      ENDIF

      DES_CROSSPRDCT(3) = XX(1)*YY(2) - XX(2)*YY(1)

    END FUNCTION DES_CROSSPRDCT

  END MODULE DISCRETELEMENT
