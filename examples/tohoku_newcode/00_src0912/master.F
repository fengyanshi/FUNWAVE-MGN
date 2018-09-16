
PROGRAM MASTER

     USE GLOBAL
     USE LOCAL
     IMPLICIT NONE
     REAL(SP) :: DT_grid1
     INTEGER :: Kgrid,Isubgrid,Nsubgrid,VTYPE
     INTEGER :: Ksubgrid_internal,NumInternal,start_count

# if defined (PARALLEL)
     call MPI_INIT ( ier )
     call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ier)
     CALL MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs, ier)
# endif

     CALL READ_INPUT
	IF(INI_UVZ)THEN
        TIME=HotStartTime
      ENDIF

     CALL ALLOCATE_VAR_TILE      ! for all grids

# if defined (PARALLEL)
     CALL PARALLEL_CARTESIAN
# endif

     CALL ALLOCATE_VARIABLES ! using maximum dimension

     CALL INDEX_LOCAL(1)   ! grid(1): mloc,nloc,ibeg, iend, etc.

     CALL INITIALIZATION   ! for grid(1) only

     CALL INITIAL_GRID(1)   ! assign initial to grid(1)

# if defined (COUPLING)
     TIME=TIME_COUPLING_1
# endif

! # if defined (VESSEL)
!     CALL VESSEL_INITIAL   !consider later
! # endif

! time integration

      !record wall time
     call WALL_TIME_SECS(tbegin)

     START_GRID(1)=.TRUE.
     START_COUNT=1

     DO Kgrid=2,NumGrid
	  CALL INDEX_LOCAL(Kgrid)
	  CALL INITIAL_GRID(Kgrid)
     ENDDO
     
     DO WHILE (TIME<TOTAL_TIME)

! every grid calculation

        DO Kgrid=1,NumGrid

           NumInternal=TOTALRATIO_SPACING(Kgrid)
           IF(Kgrid.GT.1) DT=DT_grid1/TOTALRATIO_SPACING(Kgrid)
	       
           DO Ksubgrid_internal=1,NumInternal

              ! subgrid (initial)
			IF(Ksubgrid_internal == 1)THEN   !addition
                IF(Kgrid.LT.NumGrid)THEN
	            
                    CALL INDEX_LOCAL(Kgrid+1)
                    CALL GET_NESTING_BC(Kgrid+1,1)

                ENDIF
              ENDIF
              
			CALL INDEX_LOCAL(Kgrid)
			CALL LOAD_DATA(Kgrid)  !for print subgrid at initial

                CALL CALC_GRID(Kgrid,Ksubgrid_internal,NumInternal)
                IF(Kgrid.EQ.1)DT_grid1=DT

                IF(Ksubgrid_internal == NumInternal)THEN
                  IF(Kgrid .GT. 1) THEN
	              IF(Kgrid == NumGrid) THEN

	                DO Nsubgrid=NumGrid,2,-1

	                   CALL INDEX_LOCAL( Nsubgrid )
	                   CALL LOAD_DATA( Nsubgrid )
	                   CALL TWOWAY_NESTING( Nsubgrid )					   
					    
	                ENDDO

	              ENDIF
	            ENDIF
                ENDIF

              ! subgrid (boundary)
			IF(Ksubgrid_internal == NumInternal)THEN   !addition
                IF(Kgrid.LT.NumGrid)THEN

                    CALL INDEX_LOCAL(Kgrid+1)
                    CALL GET_NESTING_BC(Kgrid+1,2)

			  ENDIF
              ENDIF

           ENDDO ! end Ksubgrid_internal

      !===========print out=================
           IF(Kgrid.EQ.1)THEN
             PLOT_COUNT=PLOT_COUNT+DT

             IF(PLOT_COUNT>=PLOT_INTV)THEN
               PLOT_COUNT=PLOT_COUNT-PLOT_INTV
               DO Isubgrid=1,NumGrid
                  IsPrint(Isubgrid)=.TRUE.
               ENDDO
             ENDIF
           ENDIF

           IF(IsPrint(Kgrid))THEN
             CALL INDEX_LOCAL(Kgrid)  ! index back
	       CALL LOAD_DATA(Kgrid)    ! ykchoi 0720 
		   CALL PREVIEW_AMR(Kgrid)
             IsPrint(Kgrid)=.FALSE.
           ENDIF
	!=====================================

	!=====================================
		 CALL INITIALIZE_VARIABLES
	!=====================================

        ENDDO  ! end kgrid

        IF(START_COUNT.LT.NumGrid)THEN
          START_COUNT=START_COUNT+1
          START_GRID(START_COUNT)=.TRUE.
        ENDIF

! --------------------------------------------- 
     END DO  ! end time

# if defined (COUPLING)
     CLOSE(11)
# endif

! # if defined (VESSEL)
! DO Kves = 201,200+NumVessel
!     CLOSE(Kves)      !consider later
! ENDDO
! # endif
      
	!record wall time at the end
     call WALL_TIME_SECS(tend)

# if defined (PARALLEL)
     if (myid.eq.0) write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     if (myid.eq.0) write(3,*) 'Simulation takes',tend-tbegin,'seconds'
     if (myid.eq.0) WRITE(*,*)'Normal Termination!'
     if (myid.eq.0) WRITE(3,*)'Normal Termination!'
# else
     write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     write(3,*) 'Simulation takes',tend-tbegin,'seconds'
     WRITE(*,*)'Normal Termination!'
     WRITE(3,*)'Normal Termination!'
# endif

# if defined (PARALLEL)
     call MPI_FINALIZE ( ier )
# endif

END PROGRAM MASTER