
!===========================================================================
!===========================================================================

! --- For AMR modeling !ykchoi
SUBROUTINE INITIAL_GRID(ng)
     USE GLOBAL  
     USE LOCAL   
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng
     INTEGER :: mbox1,nbox1,isk,mbox0,nbox0,m_move,n_move

     REAL(SP) :: Dphi_r,Dtheta_r  !ykchoi
     INTEGER :: ParMloc, ParNloc, ParMloc1, ParNloc1  !ykchoi

     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: rMASKParent, rMASKChild  !ykchoi  07/24/2018

# if defined (PARALLEL)
     INTEGER :: parent_area_i, parent_area_j
     INTEGER :: area_i, area_j
     INTEGER :: II,JJ
# endif
     
! for ng=1
     IF(ng==1)THEN
# if defined(CARTESIAN)
	 DX_Tile(ng) = DX
	 DY_Tile(ng) = DY
# else
       DX_Tile(1:Mloc,1:Nloc,ng) = DX(1:Mloc,1:Nloc)
       DY_Tile(1:Mloc,1:Nloc,ng) = DY(1:Mloc,1:Nloc)
# endif
       Depth_Tile(1:Mloc,1:Nloc,ng) = Depth(1:Mloc,1:Nloc)
       DepthX_Tile(1:Mloc1,1:Nloc,ng) = DepthX(1:Mloc1,1:Nloc)
       DepthY_Tile(1:Mloc,1:Nloc1,ng) = DepthY(1:Mloc,1:Nloc1)

       MASK_STRUC_Tile(1:Mloc,1:Nloc,ng) = MASK_STRUC(1:Mloc,1:Nloc)  !ykchoi 07/24/2018
	 MASK_Tile(1:Mloc,1:Nloc,ng) = MASK(1:Mloc,1:Nloc)
       MASK9_Tile(1:Mloc,1:Nloc,ng) = MASK9(1:Mloc,1:Nloc)
       CD_Tile(1:Mloc,1:Nloc,ng) = CD(1:Mloc,1:Nloc)                  !ykchoi 08/10/2018

       U_Tile(1:Mloc,1:Nloc,ng) = U(1:Mloc,1:Nloc)
       V_Tile(1:Mloc,1:Nloc,ng) = V(1:Mloc,1:Nloc)
       Ubar_Tile(1:Mloc,1:Nloc,ng) = Ubar(1:Mloc,1:Nloc)
       Vbar_Tile(1:Mloc,1:Nloc,ng) = Vbar(1:Mloc,1:Nloc)
       Eta_Tile(1:Mloc,1:Nloc,ng) = Eta(1:Mloc,1:Nloc)
# if defined(CARTESIAN)
!      We don't need latitude
# else
       Lat_Theta_Tile(1:Mloc,1:Nloc,ng) = Lat_Theta(1:Mloc,1:Nloc)
# endif

     ELSE

! for others
! mbox is box in a coarser grid, include ghostcells, m represents m+nghost

       !========================
	 ! case : ng >= 2
	 !========================

       isk=RATIO_SPACING(ng)
	 
       mbox1=MboxRef(ng)
       nbox1=NboxRef(ng)

# if defined (PARALLEL)
! ykchoi delete 0522
	 EastParentID(:,:,ng)=-99;   WestParentID(:,:,ng)=-99;
	 SouthParentID(:,:,ng)=-99;  NorthParentID(:,:,ng)=-99;

	 !ParentDomainID(:,:,ng)=-99;

	 parent_area_i = GridDimX( ng-1 )/px;
	 parent_area_j = GridDimY( ng-1 )/py;

	 if( n_west .eq. MPI_PROC_NULL ) THEN
	   DO J=1,Nloc
	      DO I=1,Nghost
               II = mbox1 + ( I-1+(Mloc-2*Nghost)*npx )/isk !actually npx=0 at west
               JJ = nbox1 + ( J-1+(Nloc-2*Nghost)*npy )/isk
			 !--- ykchoi(22/June/2017)
			 !These routines (II,JJ) are also used in subourtine INTERP_BC (bc_nesting.F).
		     !Thus, in ProcessorID, i,j are computed by using these II, JJ.
	         !Please see subourtine INTERP_BC before modifying these routines.

			 CALL FIND_PROCESSOR( parent_area_i, PX, II, &
	                              parent_area_j, PY, JJ, area_i, area_j )
	         
			 WestParentID(I,J,ng) = ProcessorID(area_i,area_j)
	      ENDDO
	   ENDDO
	 endif

	 if( n_east .eq. MPI_PROC_NULL ) THEN
	   DO J=1,Nloc
   	      DO I=Mloc-Nghost+1,Mloc
		     II = mbox1 + ( I-1+(Mloc-2*Nghost)*npx )/isk
               JJ = nbox1 + ( J-1+(Nloc-2*Nghost)*npy )/isk

			 CALL FIND_PROCESSOR( parent_area_i, PX, II, &
	                              parent_area_j, PY, JJ, area_i, area_j )

			 EastParentID(I+Nghost-Mloc,J,ng) = ProcessorID(area_i,area_j)  !be careful
	      ENDDO
	   ENDDO
	 endif

	 if( n_suth .eq. MPI_PROC_NULL ) then
	   DO I=1,Mloc
	      DO J=1,Nghost
               II = mbox1 + ( I-1+(Mloc-2*Nghost)*npx )/isk
               JJ = nbox1 + ( J-1+(Nloc-2*Nghost)*npy )/isk !actually npy=0 at south

			 CALL FIND_PROCESSOR( parent_area_i, PX, II, &
	                              parent_area_j, PY, JJ, area_i, area_j )
	         
			 SouthParentID(I,J,ng) = ProcessorID(area_i,area_j)
	      ENDDO
	   ENDDO
	 endif

	 if( n_nrth .eq. MPI_PROC_NULL ) then
	   DO I=1,Mloc
	      DO J=Nloc-Nghost+1,Nloc
		     II = mbox1 + ( I-1+(Mloc-2*Nghost)*npx )/isk
               JJ = nbox1 + ( J-1+(Nloc-2*Nghost)*npy )/isk

			 CALL FIND_PROCESSOR( parent_area_i, PX, II, &
	                              parent_area_j, PY, JJ, area_i, area_j )

	         NorthParentID(I,J+Nghost-Nloc,ng) = ProcessorID(area_i,area_j)  !be careful
	      ENDDO
	   ENDDO
	 endif

!	 DO J=1,Nloc
!	    DO I=1,Mloc
!	       II = mbox1 + ( I-1+(Mloc-2*Nghost)*npx )/isk !actually npx=0 at west
!		   JJ = nbox1 + ( J-1+(Nloc-2*Nghost)*npy )/isk !actually npy=0 at south

!		   CALL FIND_PROCESSOR( parent_area_i, PX, II, &
!	                            parent_area_j, PY, JJ, area_i, area_j )

!	       ParentDomainID(I,J,ng) = ProcessorID(area_i,area_j)
!	    ENDDO
!	 ENDDO
# endif

! static
       !Careful when "Mloc" on grid1 < "Mloc" on grid2
	 ! --> Dx_Tile in some index = 0
       !DX_Tile(1:Mloc,1:Nloc,ng)=DX_Tile(1:Mloc,1:Nloc,ng-1)/isk
       !DY_Tile(1:Mloc,1:Nloc,ng)=DY_Tile(1:Mloc,1:Nloc,ng-1)/isk

# if defined (CARTESIAN)
       DX_Tile(ng) = DX_Tile(ng-1)/RATIO_SPACING(ng)
	 DY_Tile(ng) = DY_Tile(ng-1)/RATIO_SPACING(ng)
# else
       Dphi_r=Dphi*pi/180.0_SP
       Dtheta_r=Dtheta*pi/180.0_SP

       Do J=1,Nloc
       Do I=1,Mloc
# if defined (PARALLEL)
          Lat_theta(I,J)=Lat_South*pi/180.0_SP-Nghost*Dtheta_r &
                         +(npy*Nglob/py+J-1)*Dtheta_r
# else
          Lat_theta(I,J)=Lat_South*pi/180.0_SP-Nghost*Dtheta_r &
                        +(J-1)*Dtheta_r
# endif
          !DX_Tile(I,J,ng) = R_earth*Dphi_r*COS(Lat_theta(I,J))/isk
          !DY_Tile(I,J,ng) = R_earth*Dtheta_r/isk
	    DX_Tile(I,J,ng) = R_earth*Dphi_r*COS(Lat_theta(I,J))/TOTALRATIO_SPACING(ng)
	    DY_Tile(I,J,ng) = R_earth*Dtheta_r/TOTALRATIO_SPACING(ng)
       ENDDO
       ENDDO
# endif

! depth
# if defined (PARALLEL)
	 ParMloc = GridDimX(ng-1)/px + 2*Nghost
	 ParNloc = GridDimY(ng-1)/py + 2*Nghost
	 ParMloc1 = ParMloc+1
	 ParNloc1 = ParNloc+1

# else
! not parallel
       !for serial !ykchoi
	 ParMloc = GridDimX(ng-1) + 2*Nghost
	 ParNloc = GridDimY(ng-1) + 2*Nghost
	 ParMloc1 = ParMloc+1
	 ParNloc1 = ParNloc+1

# endif
! Eta
       CALL LINEAR_INTERP(ParMloc,ParNloc,Mloc,Nloc,mbox1,nbox1,isk, &
                  Eta_Tile(1:ParMloc,1:ParNloc,ng-1),Eta_Tile(1:Mloc,1:Nloc,ng),ng)

	 !call mpi_finalize(ier) !ykchoi test

	 IF( ADJUSTL(trim(DEPTH_TYPE_SUB(1:3))) == 'DAT' )THEN

# if defined (PARALLEL)
         call GetFileChild(DEPTH_FILE_SUB(ng),Depth_Tile(1:Mloc,1:Nloc,ng) )

         ! re-construct Depth
         DO J=1,Nloc
         DO I=2,Mloc
            DepthX_Tile(I,J,ng)=0.5_SP*( Depth_Tile(I-1,J,ng) + Depth_Tile(I,J,ng) )
         ENDDO
         ENDDO
         DO J=1,Nloc
            DepthX_Tile(1,J,ng)=0.5_SP*(3.0_SP*Depth_Tile(1,J,ng)-Depth_Tile(2,J,ng))
            DepthX_Tile(Mloc1,J,ng)=0.5_SP*(3.0_SP*Depth_Tile(Mloc,J,ng)-Depth_Tile(Mloc-1,J,ng))
         ENDDO

         DO J=2,Nloc
         DO I=1,Mloc
            DepthY_Tile(I,J,ng)=0.5_SP*(Depth_Tile(I,J-1,ng)+Depth_Tile(I,J,ng))
         ENDDO
         ENDDO
         DO I=1,Mloc
            DepthY_Tile(I,1,ng)=0.5_SP*(3.0_SP*Depth_Tile(I,1,ng)-Depth_Tile(I,2,ng))
            DepthY_Tile(I,Nloc1,ng)=0.5_SP*(3.0_SP*Depth_Tile(I,Nloc,ng)-Depth_Tile(I,Nloc-1,ng))
         ENDDO

	   call phi_exch( DepthX_Tile(1:Mloc,1:Nloc,ng) )
	   call phi_exch( DepthY_Tile(1:Mloc,1:Nloc,ng) )

         ! For consistency with init.F
	   DO J=1,Nloc
         DO I=1,Mloc

           Depth_Tile(I,J,ng) = 0.25_SP*( DepthX_Tile(I,J,ng) + DepthX_Tile(I+1,J,ng) + &
	                          DepthY_Tile(I,J,ng) + DepthY_Tile(I,J+1,ng) )

         ENDDO
         ENDDO
	   !  consider "waterlevel" later

# else
         OPEN( 1, FILE=ADJUSTL(TRIM(DEPTH_FILE_SUB(ng))) )
	   DO J=1,Nloc
	      READ(1,*)( Depth_Tile(I,J,ng), I=1,Mloc )
	   ENDDO
         CLOSE(1)

         ! re-construct Depth
         DO J=1,Nloc
         DO I=2,Mloc
            DepthX_Tile(I,J,ng)=0.5_SP*( Depth_Tile(I-1,J,ng) + Depth_Tile(I,J,ng) )
         ENDDO
         ENDDO
         DO J=1,Nloc
            DepthX_Tile(1,J,ng)=0.5_SP*(3.0_SP*Depth_Tile(1,J,ng)-Depth_Tile(2,J,ng))
            DepthX_Tile(Mloc1,J,ng)=0.5_SP*(3.0_SP*Depth_Tile(Mloc,J,ng)-Depth_Tile(Mloc-1,J,ng))
         ENDDO

         DO J=2,Nloc
         DO I=1,Mloc
            DepthY_Tile(I,J,ng)=0.5_SP*(Depth_Tile(I,J-1,ng)+Depth_Tile(I,J,ng))
         ENDDO
         ENDDO
         DO I=1,Mloc
            DepthY_Tile(I,1,ng)=0.5_SP*(3.0_SP*Depth_Tile(I,1,ng)-Depth_Tile(I,2,ng))
            DepthY_Tile(I,Nloc1,ng)=0.5_SP*(3.0_SP*Depth_Tile(I,Nloc,ng)-Depth_Tile(I,Nloc-1,ng))
         ENDDO

         ! For consistency with init.F
	   DO J=1,Nloc
         DO I=1,Mloc

           Depth_Tile(I,J,ng) = 0.25_SP*( DepthX_Tile(I,J,ng) + DepthX_Tile(I+1,J,ng) + &
	                          DepthY_Tile(I,J,ng) + DepthY_Tile(I,J+1,ng) )

         ENDDO
         ENDDO
	   !  consider "waterlevel" later

# endif
	 ELSE

         CALL LINEAR_INTERP(ParMloc,ParNloc,Mloc,Nloc,mbox1,nbox1,isk, &
                  Depth_Tile(1:ParMloc,1:ParNloc,ng-1),Depth_Tile(1:Mloc,1:Nloc,ng),ng)

         CALL LINEAR_INTERP(ParMloc1,ParNloc,Mloc1,Nloc,mbox1,nbox1,isk, &
                  DepthX_Tile(1:ParMloc1,1:ParNloc,ng-1),DepthX_Tile(1:Mloc1,1:Nloc,ng),ng)

         CALL LINEAR_INTERP(ParMloc,ParNloc1,Mloc,Nloc1,mbox1,nbox1,isk, &
                  DepthY_Tile(1:ParMloc,1:ParNloc1,ng-1),DepthY_Tile(1:Mloc,1:Nloc1,ng),ng)
	 ENDIF
!======ykchoi 08/10/2018, 08/15/2018
       IF(IN_Cd) THEN
       CALL LINEAR_INTERP(ParMloc,ParNloc,Mloc,Nloc,mbox1,nbox1,isk, &
                  CD_Tile(1:ParMloc,1:ParNloc,ng-1),CD_Tile(1:Mloc,1:Nloc,ng),ng)
       ELSE
	 CD_Tile(:,:,ng)=Cd_fixed
       ENDIF	  

!------------------------------------------------------------------
!------------------------------------------------ ykchoi 07/24/2018
! Young-Kwang Choi
! In this version, the easiest way is applied to see propagation of water waves in water-land cells.
!==============MASK_STRUC
	 MASK_STRUC_Tile(:,:,ng) = 1  !1 - wet;  0 - dry
       ALLOCATE( rMASKParent(1:ParMloc,1:ParNloc), rMASKChild(1:Mloc,1:Nloc) )

	 rMASKParent(1:ParMloc,1:ParNloc) = real( MASK_STRUC_Tile(1:ParMloc,1:ParNloc,ng-1) )
	 rMASKChild(1:Mloc,1:Nloc) = real( MASK_STRUC_Tile(1:Mloc,1:Nloc,ng) )

       CALL LINEAR_INTERP(ParMloc,ParNloc,Mloc,Nloc,mbox1,nbox1,isk, &
                  rMASKParent(1:ParMloc,1:ParNloc),rMASKChild(1:Mloc,1:Nloc),ng)
	 
	 MASK_STRUC_Tile(1:Mloc,1:Nloc,ng) = INT( rMASKChild(1:Mloc,1:Nloc) )
	 DEALLOCATE( rMASKParent, rMASKChild )

!==============MASK
	 MASK_Tile(:,:,ng) = 1  !1 - wet;  0 - dry
       ALLOCATE( rMASKParent(1:ParMloc,1:ParNloc), rMASKChild(1:Mloc,1:Nloc) )

	 rMASKParent(1:ParMloc,1:ParNloc) = real( MASK_Tile(1:ParMloc,1:ParNloc,ng-1) )
	 rMASKChild(1:Mloc,1:Nloc) = real( MASK_Tile(1:Mloc,1:Nloc,ng) )

       CALL LINEAR_INTERP(ParMloc,ParNloc,Mloc,Nloc,mbox1,nbox1,isk, &
                  rMASKParent(1:ParMloc,1:ParNloc),rMASKChild(1:Mloc,1:Nloc),ng)
       
	 MASK_Tile(1:Mloc,1:Nloc,ng) = INT( rMASKChild(1:Mloc,1:Nloc) )
	 DEALLOCATE( rMASKParent, rMASKChild )
	 
       !DO J=1,Nloc
       !DO I=1,Mloc
       ! IF( Eta_Tile(I,J,ng) < -DEPTH_Tile(I,J,ng) ) THEN
       !  MASK_Tile(I,J,ng)=0
       !  Eta_Tile(I,J,ng)=MinDepth-Depth_Tile(I,J,ng)
       ! ELSE
       !  MASK_Tile(I,J,ng)=1
       ! ENDIF
       !ENDDO
       !ENDDO

!==============MASK9
	MASK9_Tile(:,:,ng) = 1  !1 - Boussinesq;  0 - SWE
	DO J=2,Nloc-1  !ykchoi for nesting 0526
	DO I=2,Mloc-1
      
	MASK9_Tile(I,J,ng)=MASK_Tile(I,J,ng)*MASK_Tile(I-1,J,ng)*MASK_Tile(I+1,J,ng)  &
                *MASK_Tile(I+1,J+1,ng)*MASK_Tile(I,J+1,ng)*MASK_Tile(I-1,J+1,ng) &
                *MASK_Tile(I+1,J-1,ng)*MASK_Tile(I,J-1,ng)*MASK_Tile(I-1,J-1,ng) 
      IF(ABS(Eta_Tile(I,J,ng))/MAX(DEPTH_Tile(I,J,ng),MinDepthFrc)>SWE_ETA_DEP)THEN
       MASK9_Tile(I,J,ng)=ZERO
      ENDIF

      ENDDO
      ENDDO

# if defined(PARALLEL)
      CALL PHI_INT_EXCH( MASK_STRUC_Tile(1:Mloc,1:Nloc,ng) )
	CALL PHI_INT_EXCH( MASK_Tile(1:Mloc,1:Nloc,ng) )
      CALL PHI_INT_EXCH( MASK9_Tile(1:Mloc,1:Nloc,ng) )
# endif
!------------------------------------------------ ykchoi 07/24/2018
!------------------------------------------------------------------

! U
       CALL LINEAR_INTERP(ParMloc,ParNloc,Mloc,Nloc,mbox1,nbox1,isk, &
                  U_Tile(1:ParMloc,1:ParNloc,ng-1),U_Tile(1:Mloc,1:Nloc,ng),ng)

! V
       CALL LINEAR_INTERP(ParMloc,ParNloc,Mloc,Nloc,mbox1,nbox1,isk, &
                  V_Tile(1:ParMloc,1:ParNloc,ng-1),V_Tile(1:Mloc,1:Nloc,ng),ng)

! Ubar
       CALL LINEAR_INTERP(ParMloc,ParNloc,Mloc,Nloc,mbox1,nbox1,isk, &
                  Ubar_Tile(1:ParMloc,1:ParNloc,ng-1),Ubar_Tile(1:Mloc,1:Nloc,ng),ng)

! Vbar
       CALL LINEAR_INTERP(ParMloc,ParNloc,Mloc,Nloc,mbox1,nbox1,isk, &
                  Vbar_Tile(1:ParMloc,1:ParNloc,ng-1),Vbar_Tile(1:Mloc,1:Nloc,ng),ng)

! Lat_theta
# if defined (CARTESIAN)
       
# else
	 CALL LINEAR_INTERP(ParMloc,ParNloc,Mloc,Nloc,mbox1,nbox1,isk, &
                  Lat_Theta_Tile(1:ParMloc,1:ParNloc,ng-1),Lat_Theta_Tile(1:Mloc,1:Nloc,ng),ng)
# endif
! exchange for parallel !!!

     ENDIF  ! end ng

END SUBROUTINE INITIAL_GRID

!===========================================================================
!===========================================================================

# if defined (PARALLEL)
SUBROUTINE FIND_PROCESSOR( parent_area_i, PX, II, &
	                      parent_area_j, PY, JJ, area_i, area_j )
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: parent_area_i, PX, II, parent_area_j, PY, JJ
      INTEGER,INTENT(OUT) :: area_i, area_j     
      INTEGER :: K

	DO K=1,PX
	   IF( parent_area_i*(K-1) < II .and. II <= parent_area_i*K ) THEN
	       area_i = K
	   ENDIF
	ENDDO
	DO K=1,PY
	   IF( parent_area_j*(K-1) < JJ .and. JJ <= parent_area_j*K ) THEN
	       area_j = K
	   ENDIF
	ENDDO	         
ENDSUBROUTINE FIND_PROCESSOR
# endif

!===========================================================================
!===========================================================================

! --- For AMR modeling !ykchoi
SUBROUTINE LINEAR_INTERP(MaxM,MaxN,M,N,mb,nb,isk,Fin,Fout,ng)
      USE PARAM
!ykchoi
# if defined (PARALLEL)
	USE GLOBAL, ONLY : GridDimX,GridDimY,Nghost,px,py,npx,npy,myid
# endif

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: M,N,mb,nb,isk,MaxM,MaxN,ng   !ykchoi
      INTEGER :: II,JJ
      REAL(SP) :: rII,rJJ      
      REAL(SP),DIMENSION(MaxM,MaxN),INTENT(IN) :: Fin
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Fout

!ykchoi
# if defined (PARALLEL)
      REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VarGrid1  ! global including ghost
      INTEGER :: mm1, nn1, Mloc_grid1, Nloc_grid1
# endif

!ykchoi
# if defined (PARALLEL)
      mm1 = GridDimX(ng-1) + 2*Nghost
      nn1 = GridDimY(ng-1) + 2*Nghost
      ALLOCATE( VarGrid1(mm1,nn1) )

      Mloc_grid1 = GridDimX(ng-1)/px + 2*Nghost
      Nloc_grid1 = GridDimY(ng-1)/py + 2*Nghost 
      CALL GATHER_GRID( VarGrid1, Fin(1:Mloc_grid1,1:Nloc_grid1), &
               Mloc_grid1, Nloc_grid1, mm1, nn1, Nghost)

	!if(myid.eq.0) then
	!  print *,mm1, nn1
	!  open(2005,file='tmp.out')

	!  do j=1,nn1
	!     write( 2005,'(5000E16.6)')( real(VarGrid1(i,j)), i=1,mm1 )
	!  enddo

	!  close(2005)
	!endif
# endif


# if defined (PARALLEL)

! ykchoi delete 0522
       DO J=1,N
        DO I=1,M
         II = mb + (I-1+(M-2*Nghost)*npx)/isk
         JJ = nb + (J-1+(N-2*Nghost)*npy)/isk
         rII = REAL(mb) + REAL(I-1+(M-2*Nghost)*npx)/REAL(isk) - REAL(II)
         rJJ = REAL(nb) + REAL(J-1+(N-2*Nghost)*npy)/REAL(isk) - REAL(JJ)
         Fout(I,J)  &
	     = ( (1.0_SP-rII)*VarGrid1(II,JJ) + rII*VarGrid1(II+1,JJ) )*( 1.0_SP-rJJ ) + &
             ( (1.0_SP-rII)*VarGrid1(II,JJ+1) + rII*VarGrid1(II+1,JJ+1) )*rJJ
        ENDDO
       ENDDO
	 DEALLOCATE( VarGrid1 )

# else
! not parallel

       DO J=1,N
       DO I=1,M
        II=mb+(I-1)/isk
        JJ=nb+(J-1)/isk
        rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
        rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
        Fout(I,J)=((1.0_SP-rII)*Fin(II,JJ)+rII*Fin(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1)+rII*Fin(II+1,JJ+1))*rJJ
       ENDDO
       ENDDO
# endif

END SUBROUTINE LINEAR_INTERP

!===========================================================================
!===========================================================================

! --- For AMR modeling !ykchoi

SUBROUTINE LOAD_DATA(ng)

     USE GLOBAL  
     USE LOCAL

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng

# if defined(CARTESIAN)
      ! No variables "Lat_Theta"
# else
     Lat_Theta(1:Mloc,1:Nloc) = Lat_Theta_Tile(1:Mloc,1:Nloc,ng)
# endif

     Eta(1:Mloc,1:Nloc) = Eta_Tile(1:Mloc,1:Nloc,ng)
     Ubar(1:Mloc,1:Nloc) = Ubar_Tile(1:Mloc,1:Nloc,ng)
     Vbar(1:Mloc,1:Nloc) = Vbar_Tile(1:Mloc,1:Nloc,ng)
     
     CD(1:Mloc,1:Nloc) = CD_Tile(1:Mloc,1:Nloc,ng)    !ykchoi 08/10/2018

     U(1:Mloc,1:Nloc) = U_Tile(1:Mloc,1:Nloc,ng)
     V(1:Mloc,1:Nloc) = V_Tile(1:Mloc,1:Nloc,ng)

     MASK_STRUC(1:Mloc,1:Nloc) = MASK_STRUC_Tile(1:Mloc,1:Nloc,ng)  !ykchoi 07/24/2018
     
     MASK(1:Mloc,1:Nloc) = MASK_Tile(1:Mloc,1:Nloc,ng)
     MASK9(1:Mloc,1:Nloc) = MASK9_Tile(1:Mloc,1:Nloc,ng)
     Depth(1:Mloc,1:Nloc) = Depth_Tile(1:Mloc,1:Nloc,ng)
     DepthX(1:Mloc1,1:Nloc) = DepthX_Tile(1:Mloc1,1:Nloc,ng)
     DepthY(1:Mloc,1:Nloc1) = DepthY_Tile(1:Mloc,1:Nloc1,ng)

# if defined(CARTESIAN)
     DX = DX_Tile(ng)
     DY = DY_Tile(ng)
# else
     DX(1:Mloc,1:Nloc) = DX_Tile(1:Mloc,1:Nloc,ng)
     DY(1:Mloc,1:Nloc) = DY_Tile(1:Mloc,1:Nloc,ng)
# endif

END SUBROUTINE LOAD_DATA

!===========================================================================
!===========================================================================

! --- For AMR modeling !ykchoi
SUBROUTINE CALC_GRID(ng,NestStep,NestTotal)

     USE GLOBAL  
     USE LOCAL   
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng,NestStep,NestTotal

     INTEGER::ISTAGE

! INTEGER :: mm1, nn1

! Consider later (ykchoi)
! # if defined (SAMPLES)
!     IF(WaveMaker(1:7)=='LEF_SOL')THEN
!        CALL SOLITARY_WAVE_LEFT_BOUNDARY
!     ENDIF  
! # endif 

   ! update three variables
     Eta0=Eta
     Ubar0=Ubar
     Vbar0=Vbar

     CALL UPDATE_MASK

     CALL EXCHANGE(ng)    !AMR routine
     IF(ng.GT.1)THEN
       CALL USE_NESTING_BC(NestStep,NestTotal,ng)
     ENDIF

# if defined (COUPLING)
     CALL OneWayCoupling
# endif

! calculate other variables for fluxes
     H=Eta*Gamma3+Depth
     HU=H*U*MASK
     HV=H*V*MASK

     IF(ng.EQ.1)THEN
# if defined(CARTESIAN)
       CALL ESTIMATE_DT(Mloc,Nloc,DX,DY,  &
	      U(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),H(1:Mloc,1:Nloc),MinDepthFrc,DT,CFL,TIME)
# else
       CALL ESTIMATE_DT(Mloc,Nloc,DX(1:Mloc,1:Nloc),DY(1:Mloc,1:Nloc),  &
	      U(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),H(1:Mloc,1:Nloc),MinDepthFrc,DT,CFL,TIME)
# endif
     ENDIF

     U0=U
     V0=V

! # if defined (VESSEL)
!       CALL VESSEL_FORCING    !consider later
! # endif

! 3-ORDER RUNGE-KUTTA TIME STEPPING
     DO ISTAGE=1,3

        IF(DISPERSION)THEN
          CALL Cal_Dispersion(ng)    !AMR routine
        ENDIF
	  
	  CALL FLUXES(ng)
	  
! # if defined (WIND)
!        IF(WindForce)THEN
!          CALL WIND_EFFECT    !consider later
!        ENDIF
! # endif

	  CALL SourceTerms(ng)   ! put sourceterms after fluxes in order to get eta_t 

	  CALL ESTIMATE_HUV(ISTAGE, ng)    !AMR routine

        CALL WAVE_BREAKING

! # if defined (FILTERING)
! # if defined (CARTESIAN)
!======================================================
!====== Filtering process
! To conduct Vincent and Briggs experiment
! spherical mode?
!        CALL SHAPIRO_3(Mloc,Nloc,Nghost,MASK(1:Mloc,1:Nloc),Eta(1:Mloc,1:Nloc))
        !CALL SHAPIRO_3(Mloc,Nloc,Nghost,MASK(1:Mloc,1:Nloc),U(1:Mloc,1:Nloc))
        !CALL SHAPIRO_3(Mloc,Nloc,Nghost,MASK(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc))

!        H=Eta*Gamma3+Depth

!        DO J=Jbeg,Jend
!        DO I=Ibeg,Iend   
!	     IF(MASK(I,J)==1)THEN
!		   HU(I,J)=Max(H(I,J),MinDepthFrc)*U(I,J)
!		   HV(I,J)=Max(H(I,J),MinDepthFrc)*V(I,J)
!	     ENDIF
!        ENDDO
!        ENDDO     
!====== Filtering process
!======================================================
! # endif
! # endif

        CALL EXCHANGE(ng)    !AMR routine

# if defined (COUPLING)
        CALL OneWayCoupling
# endif

        IF(ng.GT.1)THEN
          CALL USE_NESTING_BC(NestStep,NestTotal,ng)
        ENDIF

        IF(ng.Eq.1)THEN
          IF(WaveMaker(1:3)=='ABS') THEN
            CALL ABSORBING_GENERATING_BC
          ENDIF

          IF(DIRECT_SPONGE)THEN
            CALL SPONGE_DAMPING
          ENDIF
        ENDIF

     ENDDO

! # if defined (FILTERING)
! # if defined (CARTESIAN)
!======================================================
!====== Filtering process
! spherical mode?
!        CALL SHAPIRO_3(Mloc,Nloc,Nghost,MASK(1:Mloc,1:Nloc),Eta(1:Mloc,1:Nloc))
        !CALL SHAPIRO_3(Mloc,Nloc,Nghost,MASK(1:Mloc,1:Nloc),U(1:Mloc,1:Nloc))
        !CALL SHAPIRO_3(Mloc,Nloc,Nghost,MASK(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc))

!        H=Eta*Gamma3+Depth

!        DO J=Jbeg,Jend
!        DO I=Ibeg,Iend   
!	     IF(MASK(I,J)==1)THEN
!		   HU(I,J)=Max(H(I,J),MinDepthFrc)*U(I,J)
!		   HV(I,J)=Max(H(I,J),MinDepthFrc)*V(I,J)
!	     ENDIF
!        ENDDO
!        ENDDO
!	  CALL EXCHANGE(ng)    !AMR routine
!====== Filtering process
!======================================================
! # endif
! # endif     

     CALL UPDATE_MASK  !ykchoi 07/23/2018

!==========================
! save for next time step
!==========================
     MASK_TILE(1:Mloc,1:Nloc,ng) = MASK(1:Mloc,1:Nloc)
     MASK9_TILE(1:Mloc,1:Nloc,ng) = MASK9(1:Mloc,1:Nloc)
     
     MASK_STRUC_Tile(1:Mloc,1:Nloc,ng) = MASK_STRUC(1:Mloc,1:Nloc)   !ykchoi 07/24/2018

     U_TILE(1:Mloc,1:Nloc,ng) = U(1:Mloc,1:Nloc)
     V_TILE(1:Mloc,1:Nloc,ng) = V(1:Mloc,1:Nloc)
     Ubar_TILE(1:Mloc,1:Nloc,ng) = Ubar(1:Mloc,1:Nloc)
     Vbar_TILE(1:Mloc,1:Nloc,ng) = Vbar(1:Mloc,1:Nloc)
       
     Eta_TILE(1:Mloc,1:Nloc,ng) = Eta(1:Mloc,1:Nloc)

!     CALL MIXING_STUFF    !consider later

!  find maximum eta velocity
	IF(ng.eq.1)THEN 
       
	 IF (OUT_Hmax.OR.OUT_Hmin.OR.OUT_Umax.OR.OUT_MFmax.OR.OUT_VORmax)THEN
         CALL MAX_MIN_PROPERTY
       ENDIF

	 CALL OUTPUT   !STATISTICS, STATIONS
       
	 CALL CHECK_BLOWUP

	ENDIF

END SUBROUTINE CALC_GRID

!===========================================================================
!===========================================================================

SUBROUTINE INITIALIZE_VARIABLES

     USE GLOBAL
     IMPLICIT NONE

     U_NESTING_EAST_PAR = U_NESTING_EAST
     V_NESTING_EAST_PAR = V_NESTING_EAST
     Z_NESTING_EAST_PAR = Z_NESTING_EAST

     U_NESTING_WEST_PAR = U_NESTING_WEST
     V_NESTING_WEST_PAR = V_NESTING_WEST
     Z_NESTING_WEST_PAR = Z_NESTING_WEST

     U_NESTING_SOUTH_PAR = U_NESTING_SOUTH
     V_NESTING_SOUTH_PAR = V_NESTING_SOUTH
     Z_NESTING_SOUTH_PAR = Z_NESTING_SOUTH

     U_NESTING_NORTH_PAR = U_NESTING_NORTH
     V_NESTING_NORTH_PAR = V_NESTING_NORTH
     Z_NESTING_NORTH_PAR = Z_NESTING_NORTH

      !Initialize
     Eta = ZERO
     Ubar = ZERO
     Vbar = ZERO

     CD = ZERO       !ykchoi 08/10/2018

     U = ZERO
     V = ZERO
     
     MASK = 1        !ykchoi 07/24/2018
     MASK9 = 1       !ykchoi 07/24/2018
     MASK_STRUC = 1  !ykchoi 07/24/2018
	     
     Depth = ZERO
     DepthX = ZERO
     DepthY = ZERO
     DX = ZERO
     DY = ZERO

     U_NESTING_EAST = ZERO
     V_NESTING_EAST = ZERO
     Z_NESTING_EAST = ZERO

     U_NESTING_WEST = ZERO
     V_NESTING_WEST = ZERO
     Z_NESTING_WEST = ZERO

     U_NESTING_SOUTH = ZERO
     V_NESTING_SOUTH = ZERO
     Z_NESTING_SOUTH = ZERO

     U_NESTING_NORTH = ZERO
     V_NESTING_NORTH = ZERO
     Z_NESTING_NORTH = ZERO

END SUBROUTINE INITIALIZE_VARIABLES

!===========================================================================
!===========================================================================

SUBROUTINE TWOWAY_NESTING(ng)

     USE GLOBAL  
     USE LOCAL   
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng

     INTEGER :: mbox1,nbox1,isk
     INTEGER :: mloc_grid1, nloc_grid1

     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: rMASKParent, rMASKChild  !ykchoi  07/24/2018

	isk=RATIO_SPACING(ng)
	mbox1=MboxRef(ng)
	nbox1=NboxRef(ng)

# if defined(PARALLEL)
      mloc_grid1 = GridDimX(ng-1)/px + 2*Nghost
      nloc_grid1 = GridDimY(ng-1)/py + 2*Nghost
# else
      mloc_grid1 = GridDimX(ng-1) + 2*Nghost
      nloc_grid1 = GridDimY(ng-1) + 2*Nghost
# endif

	CALL LINEAR_INTERP_INV( mloc_grid1,nloc_grid1,Mloc,Nloc,mbox1,nbox1, &
		        isk,Eta_Tile(1:Mloc,1:Nloc,ng),Eta_Tile(1:mloc_grid1,1:nloc_grid1,ng-1),ng )

	CALL LINEAR_INTERP_INV( mloc_grid1,nloc_grid1,Mloc,Nloc,mbox1,nbox1, &
		        isk,U_Tile(1:Mloc,1:Nloc,ng),U_Tile(1:mloc_grid1,1:nloc_grid1,ng-1),ng )

	CALL LINEAR_INTERP_INV( mloc_grid1,nloc_grid1,Mloc,Nloc,mbox1,nbox1, &
		        isk,V_Tile(1:Mloc,1:Nloc,ng),V_Tile(1:mloc_grid1,1:nloc_grid1,ng-1),ng )

# if defined (PARALLEL)
	call phi_exch( Ubar_Tile(1:Mloc,1:Nloc,ng) )  !AMR 1:Mloc, 1:Nloc
# endif
	CALL LINEAR_INTERP_INV( mloc_grid1,nloc_grid1,Mloc,Nloc,mbox1,nbox1, &
		        isk,Ubar_Tile(1:Mloc,1:Nloc,ng),Ubar_Tile(1:mloc_grid1,1:nloc_grid1,ng-1),ng )

# if defined (PARALLEL)	
	call phi_exch( Vbar_Tile(1:Mloc,1:Nloc,ng) )  !AMR 1:Mloc, 1:Nloc
# endif
	CALL LINEAR_INTERP_INV( mloc_grid1,nloc_grid1,Mloc,Nloc,mbox1,nbox1, &
		        isk,Vbar_Tile(1:Mloc,1:Nloc,ng),Vbar_Tile(1:mloc_grid1,1:nloc_grid1,ng-1),ng )

!------- ykchoi 07/24/2018
!==============MASK
	ALLOCATE( rMASKChild(1:Mloc,1:Nloc), rMASKParent(1:mloc_grid1,1:nloc_grid1) )
	
	rMASKChild(1:Mloc,1:Nloc) = real( MASK_Tile(1:Mloc,1:Nloc,ng) )
	rMASKParent(1:mloc_grid1,1:nloc_grid1) = real( MASK_Tile(1:mloc_grid1,1:nloc_grid1,ng-1) )
	
	CALL LINEAR_INTERP_INV( mloc_grid1,nloc_grid1,Mloc,Nloc,mbox1,nbox1, &
		        isk,rMASKChild(1:Mloc,1:Nloc),rMASKParent(1:mloc_grid1,1:nloc_grid1),ng )

	MASK_Tile(1:mloc_grid1,1:nloc_grid1,ng-1) = INT( rMASKParent(1:mloc_grid1,1:nloc_grid1) )
	DEALLOCATE( rMASKChild, rMASKParent )

!==============MASK9
	ALLOCATE( rMASKChild(1:Mloc,1:Nloc), rMASKParent(1:mloc_grid1,1:nloc_grid1) )

	rMASKChild(1:Mloc,1:Nloc) = real( MASK9_Tile(1:Mloc,1:Nloc,ng) )
	rMASKParent(1:mloc_grid1,1:nloc_grid1) = real( MASK9_Tile(1:mloc_grid1,1:nloc_grid1,ng-1) )

	CALL LINEAR_INTERP_INV( mloc_grid1,nloc_grid1,Mloc,Nloc,mbox1,nbox1, &
		        isk,rMASKChild(1:Mloc,1:Nloc),rMASKParent(1:mloc_grid1,1:nloc_grid1),ng )

	MASK9_Tile(1:mloc_grid1,1:nloc_grid1,ng-1) = INT( rMASKParent(1:mloc_grid1,1:nloc_grid1) )
	DEALLOCATE( rMASKChild, rMASKParent )
!------- ykchoi 07/24/2018


# if defined (PARALLEL)
	!call phi_exch( Ubar_Tile(1:mloc_grid1,1:nloc_grid1,ng-1) )
	!call phi_exch( Vbar_Tile(1:mloc_grid1,1:nloc_grid1,ng-1) )
# endif

!	CALL LINEAR_INTERP_INV(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1, &
!		        isk,Eta_Tile(1:Mloc,1:Nloc,ng),Eta_Tile(:,:,ng-1),ng)

!	CALL LINEAR_INTERP_INV(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1, &
!		        isk,U_Tile(1:Mloc,1:Nloc,ng),U_Tile(:,:,ng-1),ng)
                    
!	CALL LINEAR_INTERP_INV(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1, &
!		        isk,V_Tile(1:Mloc,1:Nloc,ng),V_Tile(:,:,ng-1),ng)

!	CALL LINEAR_INTERP_INV(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1, &
!		        isk,Ubar_Tile(1:Mloc,1:Nloc,ng),Ubar_Tile(:,:,ng-1),ng)

!	CALL LINEAR_INTERP_INV(MaxDimX,MaxDimY,Mloc,Nloc,mbox1,nbox1, &
!		        isk,Vbar_Tile(1:Mloc,1:Nloc,ng),Vbar_Tile(:,:,ng-1),ng)

! calculate MASK and MASK9 separately
! ykchoi - temporary
!      DO J=1,nloc_grid1
!      DO I=1,mloc_grid1
!        IF( Eta_Tile(I,J,ng-1) < -DEPTH_Tile(I,J,ng-1) ) THEN
!         MASK_Tile(I,J,ng-1)=0
!         Eta_Tile(I,J,ng-1)=MinDepth-Depth_Tile(I,J,ng-1)
!        ELSE
!         MASK_Tile(I,J,ng-1)=1
!        ENDIF
!      ENDDO
!      ENDDO

!      DO J=2,nloc_grid1-1
!      DO I=2,mloc_grid1-1
!        MASK9_Tile(I,J,ng-1)=MASK_Tile(I,J,ng-1)*MASK_Tile(I-1,J,ng-1)*MASK_Tile(I+1,J,ng-1)  &
!                *MASK_Tile(I+1,J+1,ng-1)*MASK_Tile(I,J+1,ng-1)*MASK_Tile(I-1,J+1,ng-1) &
!                *MASK_Tile(I+1,J-1,ng-1)*MASK_Tile(I,J-1,ng-1)*MASK_Tile(I-1,J-1,ng-1) 
!        IF(ABS(Eta_Tile(I,J,ng-1))/MAX(DEPTH_Tile(I,J,ng-1),MinDepthFrc)>SWE_ETA_DEP)THEN
!          MASK9_Tile(I,J,ng-1)=ZERO
!        ENDIF

!      ENDDO
!      ENDDO

ENDSUBROUTINE TWOWAY_NESTING

!===========================================================================
!===========================================================================

SUBROUTINE LINEAR_INTERP_INV(MaxM,MaxN,M,N,mb,nb,isk,Fin,Fout,ng)

      USE PARAM

	USE GLOBAL, ONLY : Nghost

# if defined (PARALLEL)
!	USE GLOBAL, ONLY : GridDimX, GridDimY, px, py, npx, npy,   &
!	                   ParentDomainID, myid, ProcessorID, ier
	USE GLOBAL, ONLY : GridDimX, GridDimY, px, py, npx, npy,   &
	                   myid, ProcessorID, ier
# endif

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: MaxM,MaxN,M,N,mb,nb,isk,ng
      INTEGER :: II,JJ
      !REAL(SP) :: rII,rJJ      
      REAL(SP),DIMENSION(M,N),INTENT(IN) :: Fin
      REAL(SP),DIMENSION(MaxM,MaxN),INTENT(OUT) :: Fout

      !----ykchoi 0705
	INTEGER :: start_pts
      !----ykchoi 0705

# if defined (PARALLEL)
	INTEGER :: parent_area_i, parent_area_j

      !----ykchoi 0624
      !REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VarGridChild  ! global including ghost
      INTEGER :: MglobChild, NglobChild
	!INTEGER :: Mloc_child, Nloc_child
	INTEGER :: globi, globj, loci, locj
      !----ykchoi 0624
	
	!----ykchoi 0707
	INTEGER :: MxCompDomChild, NyCompDomChild, MglobChildIsk, NglobChildIsk
	INTEGER :: start_pts_outer
	REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VarGridChildIsk
	!----ykchoi 0707
# endif

# if defined (PARALLEL)
	parent_area_i = GridDimX( ng-1 )/px;
	parent_area_j = GridDimY( ng-1 )/py;

      !----ykchoi 0624
      MglobChild = GridDimX(ng) + 2*Nghost
      NglobChild = GridDimY(ng) + 2*Nghost
      !ALLOCATE( VarGridChild(MglobChild, NglobChild) )

	start_pts=2*isk+1	
	if( isk .ge. Nghost ) start_pts=isk+1
	if( isk .eq. 1 ) start_pts=Nghost+1	

	MxCompDomChild = MglobChild - ( start_pts - 1 ) - Nghost
	NyCompDomChild = NglobChild - ( start_pts - 1 ) - Nghost

	MglobChildIsk = ( MxCompDomChild - 1 )/isk + 1
	NglobChildIsk = ( NyCompDomChild - 1 )/isk + 1

	ALLOCATE( VarGridChildIsk(MglobChildIsk, NglobChildIsk) )

      !Mloc_child = GridDimX(ng)/px + 2*Nghost
      !Nloc_child = GridDimY(ng)/py + 2*Nghost 
      !CALL GATHER_GRID( VarGridChild, Fin(1:Mloc_child,1:Nloc_child), &
      !         Mloc_child, Nloc_child, Mglob_child, Nglob_child, Nghost)
      
	!CALL GATHER_GRID( VarGridChild, Fin, &
      !         M, N, MglobChild, NglobChild, Nghost)
	CALL GATHER_CHILDGRID( VarGridChildIsk, Fin,          &
               M, N, MglobChild, NglobChild, Nghost, isk,   & 
			 MglobChildIsk, NglobChildIsk, start_pts )

!	if(myid.eq.1) then
!	  print *,MglobChildIsk, NglobChildIsk
!	  open(2005,file='tmp0.out',status='unknown')
!	  do j=1,NglobChildIsk
!	     write( 2005,'(5000E16.6)')( real(VarGridChildIsk(i,j)), i=1,MglobChildIsk )
!	  enddo
!	  close(2005)
!	endif
!	call mpi_finalize(ier) !ykchoi test
      !----ykchoi 0624

      !call PutFile('tmp',VarGridChild(1:mm1,1:nn1))
# endif

# if defined (PARALLEL)
	!DO J=Nghost+1, Nghost+GridDimY(ng), isk
	! DO I=Nghost+1, Nghost+GridDimX(ng), isk    !Child M, N
!	DO J=Nghost+1, NglobChild-Nghost, isk
!	 DO I=Nghost+1, MglobChild-Nghost, isk    !Child M, N

	start_pts_outer=2
	if( isk .ge. Nghost ) start_pts_outer=1
	if( isk .eq. 1 ) start_pts_outer=Nghost
	
	DO J=1, NglobChildIsk
	 DO I=1, MglobChildIsk    !Child M, N

! if px=1, py=1
!	    globi = mb + (I-1)/isk
!	    globj = nb + (J-1)/isk

!	    Fout(globi,globj) = VarGridChild(I,J)

	    globi = mb + start_pts_outer + (I-1)
	    globj = nb + start_pts_outer + (J-1)
	    !loci = globi
	    !locj = globj
	    loci = globi - parent_area_i*npx
	    locj = globj - parent_area_j*npy

         !II=mb+(I-1)/isk
	   !JJ=nb+(J-1)/isk
	   !rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
	   !rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)

	    !IF( (loci .ge. Nghost+1)  .and.  (loci .le. Nghost+parent_area_i) ) THEN
	      !IF( (locj .ge. Nghost+1)  .and.  (locj .le. Nghost+parent_area_j) ) THEN
	    IF( (loci .ge. 1)  .and.  (loci .le. 2*Nghost+parent_area_i) ) THEN
	      IF( (locj .ge. 1)  .and.  (locj .le. 2*Nghost+parent_area_j) ) THEN
	    
		!if(myid.eq.0) then
	   !Fout(II,JJ) = ( (1.0_SP-rII)*VarGridChild(I,J) + rII*VarGridChild(I+1,J) )*(1.0_SP-rJJ)+ &
         !              ( (1.0_SP-rII)*VarGridChild(I,J+1) + rII*VarGridChild(I+1,J+1) )*rJJ
	        !Fout(loci,locj) = VarGridChild(I,J)
	        Fout(loci,locj) = VarGridChildIsk(I,J)
		!endif
	      ENDIF
	    ENDIF

	 ENDDO
	ENDDO
	DEALLOCATE( VarGridChildIsk )

# else
! not parallel

	start_pts=2*isk+1	
	if( isk .ge. Nghost ) start_pts=isk+1
	if( isk .eq. 1 ) start_pts=Nghost+1	

	! DO J=Nghost+1,N-Nghost,isk
	! DO I=Nghost+1,M-Nghost,isk
	DO J=start_pts, N-Nghost, isk
	DO I=start_pts, M-Nghost, isk

         II=mb+(I-1)/isk
	   JJ=nb+(J-1)/isk
	   Fout(II,JJ) = Fin(I,J)
	   !rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
	   !rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)     
	   !Fout(II,JJ) = ( (1.0_SP-rII)*Fin(I,J) + rII*Fin(I+1,J) )*(1.0_SP-rJJ)+ &
         !              ( (1.0_SP-rII)*Fin(I,J+1) + rII*Fin(I+1,J+1) )*rJJ	     
      ENDDO
      ENDDO
# endif

ENDSUBROUTINE LINEAR_INTERP_INV

! --- For AMR modeling !ykchoi
# if defined (PARALLEL)
! --------------------------------------------
!  Gather 2D variables from all processors
!  this gathering includes ghost cells
!  08/20/2013, fyshi
!  07/07/2018, Young-Kwang Choi : work for two-way nesting
! --------------------------------------------
SUBROUTINE GATHER_CHILDGRID( phi_out2, phi_in, &
           Mloc, Nloc, mm, nn, Nghost, isk, mm_isk, nn_isk, start_pts )
! mm and nn are global but include ghost cells

    USE PARAM
    USE GLOBAL, ONLY : nprocs, npx, npy, myid, ier
    IMPLICIT NONE

    integer,intent(in) :: Mloc, Nloc, mm, nn, Nghost, isk
    integer,intent(in) :: mm_isk, nn_isk, start_pts
    real(SP),dimension(Mloc,Nloc),intent(in) :: phi_in
    real(SP),dimension(mm_isk,nn_isk),intent(out) :: phi_out2
    integer,dimension(nprocs) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(Mloc,Nloc) :: xx
    real(SP),dimension(nprocs) :: xxx
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: iglob,jglob,len,n,l

    real(SP),dimension(:,:),allocatable :: phi_out

    integer :: ii,jj

    allocate( phi_out(mm,nn) )
    phi_out=0.0_SP

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    
    !put the data in master processor into the global var
    if(myid==0) then
      phi_out(1:Mloc,1:Nloc) = Phi_in(1:Mloc,1:Nloc) 
    endif
    
    !collect data from other processors into the master processor
    len = Mloc*Nloc
    do n = 1,nprocs-1
       if( myid == 0 ) then
	   call mpi_irecv(xx, len, MPI_SP, n, 0, MPI_COMM_WORLD, req(1), ier)
	   call mpi_waitall(1, req, status, ier)
	   do j=1, nloc
	   do i=1, mloc
	      iglob = npxs(n+1)*(Mloc-2*Nghost)+i
	      jglob = npys(n+1)*(Nloc-2*Nghost)+j
	      phi_out( iglob, jglob ) = xx(i,j)
	   enddo
	   enddo
	 endif
       
	 if( myid == n ) then
	   call mpi_send(phi_in, len, MPI_SP, 0, 0, MPI_COMM_WORLD, ier)	    
	 endif
    enddo

!---------ykchoi
!    if( myid == 0 ) then
!      ii=0;
!      jj=0;
!      do j = Nghost+1, nn-Nghost, isk
!	  jj=jj+1 
!        do i = Nghost+1, mm-Nghost, isk
!	     ii=ii+1
!	     phi_out2(ii,jj)=phi_out(i,j)
!        enddo
!	  ii=0
!      enddo
!    endif
!    call MPI_Barrier(MPI_COMM_WORLD,ier)	
!---------ykchoi

    !scattering to every processors
      ii=0;
      jj=0;
      !do j = Nghost+1, nn-Nghost, isk
	do j = start_pts, nn-Nghost, isk
         jj=jj+1
         !do i = Nghost+1, mm-Nghost, isk
         do i = start_pts, mm-Nghost, isk
            ii=ii+1
            if (myid.eq.0) then
	         do l=1,nprocs
	            xxx(l) = phi_out(i,j)
		     enddo
            endif
            call MPI_Scatter(xxx,1,MPI_SP,&
                     phi_out2(ii,jj),1,MPI_SP,0,MPI_COMM_WORLD,ier)
         enddo
         ii=0
      enddo
    deallocate( phi_out )
END SUBROUTINE GATHER_CHILDGRID
# endif