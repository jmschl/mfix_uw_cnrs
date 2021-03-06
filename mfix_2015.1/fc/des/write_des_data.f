!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: WRITE_DES_DATA                                          !
!  Purpose: Writing DES output in Paraview format                      !
!                                                                      !
!                                                                      !
!  Author: Jay Boyalakuntla                           Date: 26-Jul-06  !
!  Reviewer: Sreekanth Pannala                        Date: 31-Oct-06  !
!                                                                      !
!  Reviewer: J. Musser                                Date: 20-Apr-10  !
!  Comments: Split original subroutine into one for ParaView *.vtp     !
!  files, and a second for TECPLOT files *.dat.                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE WRITE_DES_DATA

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE discretelement
      USE run
      USE geometry
      USE physprop
      USE sendrecv
      USE des_bc

      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------

!-----------------------------------------------
! Functions
!-----------------------------------------------

!-----------------------------------------------

      IF (TRIM(DES_OUTPUT_TYPE) .EQ. 'TECPLOT') THEN
         CALL WRITE_DES_TECPLOT
      ELSE
         CALL WRITE_DES_VTP
      ENDIF

! Invoke at own risk

! Granular temperature subroutine should be called/calculated when
! writing DES data

      IF (DES_CALC_BEDHEIGHT) THEN
         CALL CALC_DES_BEDHEIGHT
         CALL WRITE_DES_BEDHEIGHT
      ENDIF

      IF (.FALSE.) CALL DES_GRANULAR_TEMPERATURE
      IF (.FALSE.) CALL WRITE_DES_THETA

      RETURN
      END SUBROUTINE WRITE_DES_DATA
!-----------------------------------------------



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: WRITE_DES_VTP                                           !
!  Purpose: Writing DES output in Paraview format                      !
!                                                                      !
!                                                                      !
!  Reviewer: Rahul Garg                               Dare: 01-Aug-07  !
!  Comments: Added one more output file containing averaged bed height !
!                                                                      !
!  Revision : For parallel runs added cumulative and parallel IO       !
!  Author   : Pradeep G.                                               !
!                                                                      !
!  NOTE: If the system starts with zero particles, ParaView may have   !
!  trouble showing the results. To view the results in the current     !
!  version of ParaView, Version 3.6.1:                                 !
!    i - load the *.vtp files                                          !
!   ii - filter with glyph (Filters > Common > Glyph)                  !
!        a - change glyph to sphere                                    !
!        b - change scale mode to scalar                               !
!        c - check the "Edit" Set Scale Factor Box                     !
!        d - change the value to 1.0                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE WRITE_DES_VTP

      use vtp
      use run, only: ENERGY_EQ
      use discretelement, only: S_TIME
      use discretelement, only: DES_POS_NEW, DES_VEL_NEW, DES_USR_VAR
      use discretelement, only: DES_USR_VAR, DES_USR_VAR_SIZE
      use des_thermo, only: DES_T_s_NEW
      use discretelement, only: DES_RADIUS, RO_Sol, IGLOBAL_ID
      use discretelement, only: USE_COHESION, PostCohesive
      use des_rxns, only: DES_X_s
      use run, only: ANY_SPECIES_EQ
      use param, only: DIMENSION_N_S
      USE mfix_pic, only: des_stat_wt, mppic

      use error_manager

      IMPLICIT NONE


      CHARACTER(len=10) :: lNoP
      CHARACTER(len=24) :: sTIMEc
      INTEGER :: N

      sTIMEc=''; WRITE(sTIMEc,"(ES24.16)") S_TIME

! This routine opens the VTP file and calculates send/recv information.
! It returns back the number of points as a string.
      CALL VTP_OPEN_FILE(lNoP)


! Standard VTP header information:
!----------------------------------------------------------------------/
      CALL VTP_WRITE_ELEMENT('<?xml version="1.0"?>')
      CALL VTP_WRITE_ELEMENT('<!-- Time ='//sTIMEc//'s -->')
      CALL VTP_WRITE_ELEMENT('<VTKFile type="PolyData" version="0.1" &
         &byte_order="LittleEndian">')
      CALL VTP_WRITE_ELEMENT('<PolyData>')

      CALL VTP_WRITE_ELEMENT('<Piece NumberOfPoints="'//lNoP//'" &
         &NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" &
         &NumberOfPolys="0">')


! Points are the particle identified by position:
!----------------------------------------------------------------------/
      CALL VTP_WRITE_ELEMENT('<Points>')
      CALL VTP_WRITE_DATA('Position', DES_POS_NEW)
      CALL VTP_WRITE_ELEMENT('</Points>')


! PointData are individual particle properties:
!----------------------------------------------------------------------/
      CALL VTP_WRITE_ELEMENT('<PointData Scalars="Id" &
         &Vectors="Velocity">')

      CALL VTP_WRITE_DATA('Id', IGLOBAL_ID)
      CALL VTP_WRITE_DATA('Diameter', 2.0d0*DES_RADIUS)
      CALL VTP_WRITE_DATA('Density', RO_Sol)
      CALL VTP_WRITE_DATA('Velocity', DES_VEL_NEW)

      IF(DES_USR_VAR_SIZE > 0) &
         CALL VTP_WRITE_DATA('User Defined Var', DES_USR_VAR)

!      IF(MPPIC) CALL VTP_WRITE_DATA('Statwt', DES_STAT_WT)
      IF(ENERGY_EQ) &
         CALL VTP_WRITE_DATA('Temperature', DES_T_s_NEW)

      IF(ANY_SPECIES_EQ) THEN
         DO N=1, DIMENSION_N_S
            CALL VTP_WRITE_DATA(trim(iVar('X_s',N)), DES_X_s(:,N))
         ENDDO
      ENDIF

      IF(USE_COHESION) &
         CALL VTP_WRITE_DATA('CohesiveForce', PostCohesive)

      CALL VTP_WRITE_ELEMENT('</PointData>')

! Open/Close the unused VTP tags.
!----------------------------------------------------------------------/
      CALL VTP_WRITE_ELEMENT('<CellData></CellData>')
      CALL VTP_WRITE_ELEMENT('<Verts></Verts>')
      CALL VTP_WRITE_ELEMENT('<Lines></Lines>')
      CALL VTP_WRITE_ELEMENT('<Strips></Strips>')
      CALL VTP_WRITE_ELEMENT('<Polys></Polys>')


! Close all the opened tags:
!----------------------------------------------------------------------/
      CALL VTP_WRITE_ELEMENT('</Piece>')
      CALL VTP_WRITE_ELEMENT('</PolyData>')
      CALL VTP_WRITE_ELEMENT('</VTKFile>')

      CALL VTP_CLOSE_FILE

! Add the new VTP file to the PVD file for time association.
      CALL ADD_VTP_TO_PVD
      
      call WRITE_force_chain      


      RETURN
      END SUBROUTINE WRITE_DES_VTP

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111      
      
      
      SUBROUTINE WRITE_force_chain

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE discretelement
      USE run
      USE geometry
      USE physprop
      USE sendrecv
      USE des_bc
      use mpi_utility
      use compar 
      use desmpi
      USE mpi_comm_des      
      use cdist
      Use des_thermo
      IMPLICIT NONE
	  
      integer :: l, k, des_unit, l_pair_num 
      CHARACTER*64 :: FNAME_VTP = ''
      double precision, dimension(:,:), allocatable :: l_fc_p1, l_fc_p2
      double precision, dimension(:), allocatable :: l_force_n, l_force_t, l_force_c      
      integer llocalcnt,lglocnt,lgathercnts(0:numpes-1),lproc,ltotvar,lcount
      real,dimension(:,:), allocatable :: ltemp_array
      integer :: wdimn
	  
	wdimn = merge(2,3,no_k)
	ltotvar = merge(7,9,no_k)
              
        l_pair_num = 0 
	do l=1, pair_num
	      if((use_cohesion.and.force_mag_pair(3,l).gt.zero) .or. &
	        (.not.use_cohesion.and.force_mag_pair(1,l).gt.zero))then	  
          l_pair_num = l_pair_num + 1
          endif	  
	enddo
	
        allocate( l_fc_p1(wdimn,l_pair_num))
        allocate( l_fc_p2(wdimn,l_pair_num))
        allocate( l_force_n(l_pair_num))
        allocate( l_force_t(l_pair_num))
        allocate( l_force_c(l_pair_num))
        
        l_fc_p1 = 0.d0
        l_fc_p2 = 0.d0
        l_force_n = 0.d0
        l_force_t = 0.d0         
	
        k=0
	do l=1, pair_num
	  if(force_mag_pair(1,l).gt.zero)then
	  k = k + 1
	  l_fc_p1(:,k) = des_pos_new(:,pairs(1,l))
	  l_fc_p2(:,k) = des_pos_new(:,pairs(2,l))
          l_force_n(k) = force_mag_pair(1,l)
          l_force_t(k) = force_mag_pair(2,l)
          l_force_c(k) = force_mag_pair(3,l)          
          endif	  
	enddo    
	
	
	lglocnt = 10
	llocalcnt = l_pair_num
	call global_sum(llocalcnt,lglocnt)
	allocate (dprocbuf(llocalcnt),drootbuf(lglocnt),iprocbuf(llocalcnt),irootbuf(lglocnt))
	allocate (ltemp_array(lglocnt,ltotvar))
	
	igath_sendcnt = llocalcnt
	lgathercnts = 0
	lgathercnts(mype) = llocalcnt
	call global_sum(lgathercnts,igathercnts)
	idispls(0) = 0
	do lproc = 1,numpes-1
	  idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
	end do

! gather information from all processor
         lcount = 1
         do k = 1,wDIMN     
            dprocbuf(1:llocalcnt) = l_fc_p1(k,1:llocalcnt)
            call desmpi_gatherv(ptype=2)
            ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         end do
         
         do k = 1,wDIMN
            dprocbuf(1:llocalcnt) = l_fc_p2(k,1:llocalcnt)
            call desmpi_gatherv(ptype=2)
            ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         end do
         
         dprocbuf(1:llocalcnt) = l_force_n(1:llocalcnt)
         call desmpi_gatherv(ptype=2)
         ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         
         dprocbuf(1:llocalcnt) = l_force_t(1:llocalcnt)
         call desmpi_gatherv(ptype=2)
         ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         
         dprocbuf(1:llocalcnt) = l_force_c(1:llocalcnt)
         call desmpi_gatherv(ptype=2)
         ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         
      
      if (mype.eq.pe_io)then	
          write(fname_vtp,'(A,"_FC_",I5.5,".vtp")') trim(run_name),vtp_findex-1 	  
	  des_unit =99
	  open(des_unit,file=fname_vtp,status='unknown') 
     
            write(des_unit,"(a)") '<?xml version="1.0"?>'
            write(des_unit,"(a,es24.16,a)") '<!-- time =',s_time,'s -->'
            write(des_unit,"(a)") '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
            write(des_unit,"(3x,a)") '<PolyData>'
            write(des_unit,"(6x,a,i10.10,a,a,i10.10,a)")&
               '<Piece NumberOfPoints="',2*lglocnt,'" NumberOfVerts="0" ',&
               'NumberOfLines="',lglocnt,'" NumberOfStrips="0" NumberOfPolys="0">'
            write(des_unit,"(9x,a)")'<PointData Scalars="f_mag">'
            write(des_unit,"(9x,a)")'<DataArray type="Float32" Name="f_mag_n" format="ascii">'
	    do l=1, lglocnt
		write(des_unit,"(15x,es12.5)")ltemp_array(l,2*wdimn+1)
		write(des_unit,"(15x,es12.5)")ltemp_array(l,2*wdimn+1)		
            enddo			
            write(des_unit,"(12x,a)") '</DataArray>'
			
            write(des_unit,"(9x,a)")'<DataArray type="Float32" Name="f_mag_t" format="ascii">'
	    do l=1, lglocnt
		write(des_unit,"(15x,es12.5)")ltemp_array(l,2*wdimn+2)
		write(des_unit,"(15x,es12.5)")ltemp_array(l,2*wdimn+2)		
            enddo				
            write(des_unit,"(12x,a)") '</DataArray>'	
            
            write(des_unit,"(9x,a)")'<DataArray type="Float32" Name="f_mag_c" format="ascii">'
	    do l=1, lglocnt
		write(des_unit,"(15x,es12.5)")ltemp_array(l,2*wdimn+3)
		write(des_unit,"(15x,es12.5)")ltemp_array(l,2*wdimn+3)		
            enddo				
            write(des_unit,"(12x,a)") '</DataArray>'            
			
            write(des_unit,"(9x,a)")'<DataArray type="Float32" Name="f_mag" format="ascii">'
	    do l=1, lglocnt
		write(des_unit,"(15x,es12.5)")sqrt(ltemp_array(l,2*wdimn+1)**2+ltemp_array(l,2*wdimn+2)**2)
		write(des_unit,"(15x,es12.5)")sqrt(ltemp_array(l,2*wdimn+1)**2+ltemp_array(l,2*wdimn+2)**2)		
            enddo			
            
            write(des_unit,"(12x,a)") '</DataArray>'
            write(des_unit,"(12x,a)") '</PointData>'
            write(des_unit,"(9x,a)") '<CellData></CellData>'
            write(des_unit,"(9x,a)") '<Points>'
            write(des_unit,"(12x,a)") '<DataArray type="Float32" Name="Position" NumberOfComponents="3" format="ascii">'
            
            if(wdimn.eq.2)then
	    do l=1, lglocnt	        
		write(des_unit,"(15x,es12.5)")ltemp_array(l,1),ltemp_array(l,2),0.0
		write(des_unit,"(15x,es12.5)")ltemp_array(l,3),ltemp_array(l,4),0.0		
            enddo
            else
            do l=1, lglocnt            	    
		write(des_unit,"(15x,es12.5)")ltemp_array(l,1),ltemp_array(l,2),ltemp_array(l,3)
		write(des_unit,"(15x,es12.5)")ltemp_array(l,4),ltemp_array(l,5),ltemp_array(l,6)
            enddo		
            endif

            write(des_unit,"(12x,a,/9x,a)")'</DataArray>','</Points>'
            
            write(des_unit,"(9x,a)") '<Lines>'
            write(des_unit,"(9x,a)") '<DataArray type="Int32" ', &
			   'Name="connectivity" format="ascii">'
	    do l=1, lglocnt
	        write(des_unit,"(5x,i10)")(l-1)*2,(l-1)*2+1
            enddo
            write(des_unit,"(12x,a,/9x,a)")'</DataArray>'
            
            write(des_unit,"(9x,a)") '<DataArray type="Int32" ', &
			   'Name="offsets" format="ascii">'
	    do l=1, lglocnt
	        write(des_unit,"(5x,i10)")l*2
            enddo
            write(des_unit,"(12x,a,/9x,a)")'</DataArray>','</Lines>'
            write(des_unit,"(9x,a,/9x,a,/9x,a,/9x,a)")'<Verts></Verts>',&
                '<Strips></Strips>','<Polys></Polys>'
            write(des_unit,"(6x,a,/3x,a,/a)")&
               '</Piece>','</PolyData>','</VTKFile>'
		
        close(des_unit)	
        endif
        deallocate (l_fc_p1,l_fc_p2,l_force_n,l_force_t,l_force_c,dprocbuf,drootbuf,iprocbuf,irootbuf,ltemp_array)        
	
		RETURN
      END SUBROUTINE WRITE_force_chain	

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: WRITE_DES_TECPLOT
!  Purpose: Writing DES output in TECPLOT format
!
!  Revision: For parallel runs added distributed and single IO
!  Comment: In earlier version the time instances are keep appended to
!           the tecplot file. This will make tecplot file so large for
!           large simulations. Hence seperate files are written for
!           each instances
!  Author : Pradeep G.
!
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE WRITE_DES_TECPLOT

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE discretelement
      USE run
      USE geometry
      USE physprop
      USE sendrecv
      USE des_bc
      USE compar
      USE cdist
      USE desmpi
      USE mpi_comm_des
      USE mpi_utility
      USE functions
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! file units for
! Pradeep remove parameter from following variables
      INTEGER:: DES_DATA,DES_EX,DES_EPS

! output file for basic DES variables including: position, velocity,
! radius, density, mark (flag)
      CHARACTER(LEN=50)     :: FNAME_DATA

! output file for extra DES variables including:
! solids time (S_TIME), maximum neighbor count, maximum overlap
! granular energy and granular temperature
      CHARACTER(LEN=50)     :: FNAME_EXTRA

! output file for axial solids volume fraction and granular temp
      CHARACTER(LEN=50)     :: FNAME_EPS

! dummy indices
      INTEGER L, K

! index to track accounted for particles
      INTEGER PC

! Variables related to gathering info at PE_IO
      integer llocalcnt,lglocnt,lgathercnts(0:numpes-1),lproc,ltotvar,lcount
      real,dimension(:,:), allocatable :: ltemp_array

      INTEGER :: wDIMN

! Set output dimnensions
      wDIMN = merge(2,3,NO_K)
! set the total variable based on dimension
      ltotvar = merge(8,9,NO_K)

! set the file name and unit number and open file
      des_data = 2000
      des_ex = 2100
      des_eps = 2200
      if (bdist_io) then
         write(fname_data,'(A,"_DES_DATA",I4.4,"_",I4.4,".dat")') trim(run_name),tecplot_findex,mype
         write(fname_extra,'(A,"_DES_EXTRA",I4.4,"_",I4.4,".dat")') trim(run_name),tecplot_findex,mype
         write(fname_eps,'(A,"_DES_EPS",I4.4,"_",I4.4,".dat")') trim(run_name),tecplot_findex,mype
         open(unit=des_data,file=fname_data,status='new',err=999)
         open(unit=des_ex,file=fname_extra,status='new',err=999)
         open(unit=des_eps,file=fname_eps,status='new',err=999)
      else
         if(mype.eq.pe_io) then
            write(fname_data,'(A,"_DES_DATA_",I4.4,".dat")') trim(run_name),tecplot_findex
            write(fname_extra,'(A,"_DES_EXTRA_",I4.4,".dat")') trim(run_name),tecplot_findex
            write(fname_eps,'(A,"_DES_EPS_",I4.4,".dat")') trim(run_name),tecplot_findex
            open(unit=des_data,file=fname_data,status='new',err=999)
            open(unit=des_ex,file=fname_extra,status='new',err=999)
            open(unit=des_eps,file=fname_eps,status='new',err=999)
         end if
      end if
      tecplot_findex = tecplot_findex + 1

! write header
      if (bdist_io .or. mype .eq. pe_io) then
         if(DO_K) then
            write (des_data,'(A)') &
            'variables = "x" "y" "z" "vx" "vy" "vz" "rad" "den" "mark"'
         else
            write (des_data, '(A)') &
            'variables = "x" "y" "vx" "vy" "omega" "rad" "den" "mark"'
         endif
         write (des_data, "(A,F15.7,A)")'zone T ="', s_time, '"'
      end if

! write data in ditributed mode or as a single file
      if(bdist_io) then
         pc = 1
         do l = 1,max_pip
            if(pc.gt.pip) exit
            if(.not.pea(l,1)) cycle
            pc = pc+1
            if(pea(l,4)) cycle
            if(DO_K) then
               write (des_data, '(8(2x,es12.5))')&
                  (des_pos_new(k,l),k=1,wDIMN),(des_vel_new(k,l),k=1,wDIMN), &
                   des_radius(l),ro_sol(l)
            else
               write (des_data, '(7(2x,es12.5))')&
                  (des_pos_new(k,l),k=1,wDIMN), (des_vel_new(k,l),k=1,wDIMN), &
                  omega_new(1,l), des_radius(l), ro_sol(l)
            endif
        end do

      else ! if bdist_io
! set parameters required for gathering info at PEIO and write as single file
         lglocnt = 10
         llocalcnt = pip - ighost_cnt
         call global_sum(llocalcnt,lglocnt)
         allocate (dprocbuf(llocalcnt),drootbuf(lglocnt),iprocbuf(llocalcnt),irootbuf(lglocnt))
         allocate (ltemp_array(lglocnt,ltotvar))
         igath_sendcnt = llocalcnt
         lgathercnts = 0
         lgathercnts(mype) = llocalcnt
         call global_sum(lgathercnts,igathercnts)
         idispls(0) = 0
         do lproc = 1,numpes-1
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
         end do

! gather information from all processor
         lcount = 1
         do k = 1,wDIMN
            call des_gather(des_pos_new(k,:))
            ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         end do
         do k = 1,wDIMN
            call des_gather(des_vel_new(k,:))
            ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         end do
         if(NO_K) then
            call des_gather(omega_new(1,:))
            ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         end if
         call des_gather(des_radius)
         ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1
         call des_gather(ro_sol)
         ltemp_array(:,lcount) = drootbuf(:); lcount=lcount+1

! write the data into file
         if (mype.eq.pe_io) then
            if(DO_K) then
               do l =1,lglocnt
                  write (des_data,'(8(2x,es12.5),I5)') (ltemp_array(l,k),k=1,8),int(ltemp_array(l,9))
               end do
            else
               do l =1,lglocnt
                  write (des_data,'(7(2x,es12.5),I5)') (ltemp_array(l,k),k=1,7),int(ltemp_array(l,8))
               end do
            end if
         end if
         deallocate (dprocbuf,drootbuf,iprocbuf,irootbuf,ltemp_array)
      end if
! close the files
      if (mype.eq.pe_io .or. bdist_io) then
         close(des_data)
         close(des_ex)
         close(des_eps)
      end if
      return

  999 write(*,"(/1x,70('*'),//,a,/,a,/1x,70('*'))")&
         ' From: write_des_tecplot ',&
         ' message: error opening des tecplot file. terminating run.'


      END SUBROUTINE WRITE_DES_TECPLOT
!-----------------------------------------------

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: WRITE_DES_BEDHEIGHT
!  Purpose: Writing DES output on bed height.

!  WARNING: This code is out-of-date and should be modified for consistency
!  with current DEM version.  Also this routine will be fairly specific
!  to a user needs and should probably be tailored as such

!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE WRITE_DES_BEDHEIGHT

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE discretelement
      USE run
      USE geometry
      USE physprop
      USE sendrecv
      USE des_bc
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! logical that identifies that the data file has been created
! and is already opened (initial checks/writes)
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
! logical used for testing is the data files already exists
      LOGICAL :: F_EXISTS
! output file for the bed height data
      CHARACTER(LEN=50)     :: FNAME_BH
! file unit for the bed height data
      INTEGER, PARAMETER :: BH_UNIT = 2010
! dummy index values
      INTEGER I, M
! variables for bed height calculation
      INTEGER, SAVE :: tcount = 1
      DOUBLE PRECISION :: height_avg, height_rms
      DOUBLE PRECISION, PARAMETER :: tmin = 5.d0
      DOUBLE PRECISION, DIMENSION(5000), SAVE :: bed_height_time

!-----------------------------------------------
! Functions
!-----------------------------------------------
!-----------------------------------------------


! after tmin start storing bed height. after enough measurements
! have been taken (i.e. tcount > 20) start to calculate a running
! average bed height and running rms bed height for solids phase 1 only
      height_avg = zero
      height_rms = zero

      if(time.gt.tmin) then
         if(tcount.le.5000)  then
            bed_height_time(tcount) = bed_height(1)
            !dt_time(tcount) = DT
            tcount = tcount + 1

            if(tcount.gt.20)  then
               do i = 1, tcount-1,1
                  height_avg = height_avg + bed_height_time(i)!*dt_time(i)
               enddo
               height_avg = height_avg/(tcount-1)
               do i = 1, tcount-1,1
                  height_rms = height_rms + ((bed_height_time(i)&
                       &-height_avg)**2)!*dt_time(i)
               enddo

               height_rms = sqrt(height_rms/(tcount-1))
            endif
         endif
      endif

      FNAME_BH = TRIM(RUN_NAME)//'_DES_BEDHEIGHT.dat'
      IF(FIRST_PASS) THEN
         F_EXISTS = .FALSE.
         INQUIRE(FILE=FNAME_BH,EXIST=F_EXISTS)
! If the file does not exist, then create it with the necessary
! header information.
         IF (.NOT.F_EXISTS) THEN
            OPEN(UNIT=BH_UNIT,FILE=FNAME_BH,&
               FORM="formatted",STATUS="new")
         ELSE
! To prevent overwriting existing files accidently, exit if the file
! exists and this is a NEW run.
            IF(RUN_TYPE .EQ. 'NEW') THEN
               WRITE(*,1000)
               WRITE(UNIT_LOG, 1000)
               CALL MFIX_EXIT(myPE)
            ELSE
! Open the file for appending of new data (RESTART_1 Case)
               OPEN(UNIT=BH_UNIT,FILE=FNAME_BH,POSITION="append")
            ENDIF
         ENDIF
         FIRST_PASS = .FALSE.
      ELSE
! Open the file and mark for appending
         OPEN(UNIT=BH_UNIT,FILE=FNAME_BH,POSITION="append")
      ENDIF

      WRITE(BH_UNIT, '(10(2X,E20.12))') s_time, &
         (bed_height(M), M=1,DES_MMAX), height_avg, height_rms
! Close the file and keep
      CLOSE(BH_UNIT, STATUS="KEEP")

      RETURN

 1000 FORMAT(/1X,70('*')//, ' From: WRITE_DES_BEDHEIGHT',/,&
         ' Message: bed_height.dat already exists in the run',&
         ' directory.',/10X, 'Run terminated to prevent',&
         ' accidental overwriting of files.',/1X,70('*')/)

      END SUBROUTINE WRITE_DES_BEDHEIGHT
!-----------------------------------------------



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: WRITE_DES_THETA
!  Purpose: The following code writes out des_theta to a file for each
!  ijk cell in the system each time des_granular_temperature is called.
!
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE WRITE_DES_THETA

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE discretelement
      USE run
      USE geometry
      USE physprop
      USE sendrecv
      USE des_bc
      USE functions
      IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! indices
      INTEGER I, J, K, IJK
!
      INTEGER M, NP
! logical that identifies that the data file has been created
! and is already opened (initial checks/writes)
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS
! file unit for the granular temperature data
      INTEGER, PARAMETER :: GT_UNIT = 2020
! output file for the granular temperature data
      CHARACTER(LEN=50)  :: FNAME_GT
!-----------------------------------------------

      FNAME_GT = TRIM(RUN_NAME)//'_DES_THETA.dat'
      IF (FIRST_PASS) THEN
         F_EXISTS = .FALSE.
         INQUIRE(FILE=FNAME_GT,EXIST=F_EXISTS)

         IF (.NOT.F_EXISTS) THEN
! If the file does not exist, then create it with the necessary
! header information.
            OPEN(UNIT=GT_UNIT,FILE=FNAME_GT,STATUS='NEW')
         ELSE
            IF(RUN_TYPE .EQ. 'NEW') THEN
! If the run is new and the GT file already exists replace it with a
! new file.
!               OPEN(UNIT=GT_UNIT,FILE=FNAME_GT,STATUS='REPLACE')
! Prevent overwriting an existing file by exiting if the file exists
! and this is a NEW run.
               WRITE(*,1001) FNAME_GT
               WRITE(UNIT_LOG,1001) FNAME_GT
               CALL MFIX_EXIT(myPE)
            ELSE
! Open the file for appending of new data (RESTART_1 Case)
               OPEN(UNIT=GT_UNIT, FILE=FNAME_GT, POSITION='APPEND')
            ENDIF
         ENDIF
         FIRST_PASS =  .FALSE.
      ELSE
! Open file and mark for appending
         OPEN(UNIT=GT_UNIT,FILE=FNAME_GT,POSITION='APPEND')
      ENDIF   ! endif (first_pass)

      WRITE(GT_UNIT,*) ''
      WRITE(GT_UNIT,'(A6,ES24.16)') 'Time=', S_TIME
      WRITE(GT_UNIT,'(A6,2X,3(A6,2X),A8)',ADVANCE="NO") 'IJK', &
         'I', 'J', 'K', 'NP'
      DO M = 1,DES_MMAX
         WRITE(GT_UNIT,'(7X,A6,I1)',ADVANCE="NO") 'THETA_',M
      ENDDO
      WRITE(GT_UNIT,*) ''
      DO IJK = IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            NP = PINC(IJK)
            WRITE(GT_UNIT,'(I6,2X,3(I6,2X),I8,(2X,ES15.5))') &
               IJK, I, J, K, NP, (DES_THETA(IJK,M), M = 1,DES_MMAX)
         ENDIF
      ENDDO

! Close the file and keep
      CLOSE(GT_UNIT, STATUS='KEEP')

      RETURN

 1001 FORMAT(/1X,70('*')//, ' From: WRITE_DES_THETA',/,&
         ' Message: ', A, ' already exists in the run',/10X,&
         'directory. Run terminated to prevent accidental overwriting',&
         /10X,'of files.',/1X,70('*')/)

      END SUBROUTINE WRITE_DES_THETA
!-----------------------------------------------
