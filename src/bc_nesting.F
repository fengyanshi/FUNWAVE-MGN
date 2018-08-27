! --- For AMR modeling !ykchoi
SUBROUTINE GET_NESTING_BC(ng,istep)
    USE GLOBAL
    USE LOCAL
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng,istep
    INTEGER :: mbox1,nbox1,isk

    isk=RATIO_SPACING(ng)
    mbox1=MboxRef(ng)
    nbox1=NboxRef(ng)

       CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Nghost,mbox1,nbox1,isk,U_Tile(:,:,ng-1), &
               U_NESTING_WEST(1:Nloc,1:Nghost,istep),U_NESTING_EAST(1:Nloc,1:Nghost,istep),  &
               U_NESTING_SOUTH(1:Mloc,1:Nghost,istep),U_NESTING_NORTH(1:Mloc,1:Nghost,istep), ng)
       CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Nghost,mbox1,nbox1,isk,V_Tile(:,:,ng-1), &
               V_NESTING_WEST(1:Nloc,1:Nghost,istep),V_NESTING_EAST(1:Nloc,1:Nghost,istep),  &
               V_NESTING_SOUTH(1:Mloc,1:Nghost,istep),V_NESTING_NORTH(1:Mloc,1:Nghost,istep), ng)
       CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Nghost,mbox1,nbox1,isk,Eta_Tile(:,:,ng-1), &
               Z_NESTING_WEST(1:Nloc,1:Nghost,istep),Z_NESTING_EAST(1:Nloc,1:Nghost,istep),  &
               Z_NESTING_SOUTH(1:Mloc,1:Nghost,istep),Z_NESTING_NORTH(1:Mloc,1:Nghost,istep), ng)

END SUBROUTINE GET_NESTING_BC

!==================================================================

! --- For AMR modeling !ykchoi
SUBROUTINE INTERP_BC(MaxM,MaxN,M,N,Nghost,mb,nb,isk,Fin,  &
                    Fout_west,Fout_east,Fout_south,Fout_north, ng)
      USE PARAM

! ykchoi
# if defined (PARALLEL)
    USE GLOBAL, ONLY : n_west, n_east, n_suth, n_nrth, &
                       px, py, npx, npy, GridDimX, GridDimY, &
					 myid, ProcessorID, ier, &
					 WestParentID, EastParentID, &
					 SouthParentID, NorthParentID
# endif

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: M,N,mb,nb,isk,MaxM,MaxN,Nghost, ng
      INTEGER :: II,JJ,IN,JN
      REAL(SP) :: rII,rJJ
      REAL(SP),DIMENSION(MaxM,MaxN),INTENT(IN) :: Fin
      REAL(SP),DIMENSION(N,Nghost),INTENT(OUT) :: Fout_west,Fout_east
      REAL(SP),DIMENSION(M,Nghost),INTENT(OUT) :: Fout_south,Fout_north

! ykchoi
# if defined (PARALLEL)
!    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VarGrid1  ! global including ghost
!    INTEGER :: mm1, nn1, Mloc_grid1, Nloc_grid1

    INTEGER :: parent_area_i, parent_area_j
    INTEGER :: len,loci,locj,ireq
    INTEGER :: istatus(mpi_status_size)
    INTEGER,DIMENSION(:,:),ALLOCATABLE :: ParentLocalID
    REAL(SP),DIMENSION(:),ALLOCATABLE :: xx
# endif

# if defined (PARALLEL)
!      mm1 = GridDimX(ng-1) + 2*Nghost
!      nn1 = GridDimY(ng-1) + 2*Nghost
!      ALLOCATE( VarGrid1(mm1,nn1) )

!      Mloc_grid1 = GridDimX(ng-1)/px + 2*Nghost
!      Nloc_grid1 = GridDimY(ng-1)/py + 2*Nghost 
!      CALL GATHER_GRID( VarGrid1, Fin(1:Mloc_grid1,1:Nloc_grid1), &
!               Mloc_grid1, Nloc_grid1, mm1, nn1, Nghost)
# endif

! ykchoi - parallel
# if defined (PARALLEL)
	parent_area_i = GridDimX( ng-1 )/px;
	parent_area_j = GridDimY( ng-1 )/py;
# endif

! --- West
# if defined (PARALLEL)
	len = N*Nghost
	ALLOCATE( ParentLocalID(1:Nghost,1:N), xx(4) ); ParentLocalID=-99;  xx=0.0_SP

	DO K=1,py

	   IF( myid == ProcessorID(1,K) ) THEN
	     ParentLocalID(1:Nghost,1:N)=WestParentID(1:Nghost,1:N,ng)
	     CALL MPI_BCAST(ParentLocalID,len,MPI_INTEGER,ProcessorID(1,K),MPI_COMM_WORLD,ier)

         ELSE
		 CALL MPI_BCAST(ParentLocalID,len,MPI_INTEGER,ProcessorID(1,K),MPI_COMM_WORLD,ier)
	   ENDIF

	   DO J=1,N
	      DO I=1,Nghost
	         II = mb + ( I-1 )/isk
	         JJ = nb + ( J-1+(N-2*Nghost)*(K-1) )/isk
               rII = REAL(mb) + REAL(I-1)/REAL(isk) - REAL(II)
               rJJ = REAL(nb) + REAL(J-1+(N-2*Nghost)*(K-1))/REAL(isk) - REAL(JJ)

			 IF( myid == ParentLocalID(I,J) ) THEN
	           loci = II - parent_area_i*npx 
			   locj = JJ - parent_area_j*npy
	           
			   xx(1) = Fin(loci, locj);
			   xx(2) = Fin(loci+1, locj);
			   xx(3) = Fin(loci, locj+1);
			   xx(4) = Fin(loci+1, locj+1);

			   CALL MPI_ISEND(xx,4,MPI_SP,ProcessorID(1,K),10,MPI_COMM_WORLD,ireq,ier)
			   CALL MPI_WAIT( ireq, istatus, ier )
	
	         ENDIF

			 IF( myid == ProcessorID(1,K) ) THEN
	           CALL MPI_IRECV(xx,4,MPI_SP,ParentLocalID(I,J),10,MPI_COMM_WORLD,ireq,ier)
			   CALL MPI_WAIT( ireq, istatus, ier )
			   
			   Fout_west(J,I)  &
	           = ( (1.0_SP-rII)*xx(1) + rII*xx(2) )*( 1.0_SP-rJJ ) + &
                   ( (1.0_SP-rII)*xx(3) + rII*xx(4) )*rJJ

			 ENDIF
	      ENDDO
	   ENDDO

	ENDDO
	DEALLOCATE( ParentLocalID, xx )
			
!     if( n_west .eq. MPI_PROC_NULL ) then
!       DO J=1,N
!        DO I=1,Nghost
!         II = mb + (I-1+(M-2*Nghost)*npx)/isk  ! actually npx=0 at west
!         JJ = nb + (J-1+(N-2*Nghost)*npy)/isk
!         rII = REAL(mb) + REAL(I-1+(M-2*Nghost)*npx)/REAL(isk) - REAL(II)
!         rJJ = REAL(nb) + REAL(J-1+(N-2*Nghost)*npy)/REAL(isk) - REAL(JJ)
!         Fout_west(J,I)  &
!	     = ( (1.0_SP-rII)*VarGrid1(II,JJ) + rII*VarGrid1(II+1,JJ) )*( 1.0_SP-rJJ ) + &
!             ( (1.0_SP-rII)*VarGrid1(II,JJ+1) + rII*VarGrid1(II+1,JJ+1) )*rJJ
!        ENDDO
!       ENDDO
!     endif

# else
! not parallel
       DO J=1,N
        DO I=1,Nghost
         II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_west(J,I)=((1.0_SP-rII)*Fin(II,JJ)+rII*Fin(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1)+rII*Fin(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO
# endif


! --- East
# if defined (PARALLEL)
	len = N*Nghost
	ALLOCATE( ParentLocalID(1:Nghost,1:N), xx(4) ); ParentLocalID=-99;  xx=0.0_SP
	
	DO K=1,py

	   IF( myid == ProcessorID(px,K) ) THEN
	     ParentLocalID(1:Nghost,1:N)=EastParentID(1:Nghost,1:N,ng)
	     CALL MPI_BCAST(ParentLocalID,len,MPI_INTEGER,ProcessorID(px,K),MPI_COMM_WORLD,ier)

         ELSE
		 CALL MPI_BCAST(ParentLocalID,len,MPI_INTEGER,ProcessorID(px,K),MPI_COMM_WORLD,ier)
	   ENDIF
 
	   DO J=1,N
	      DO I=M-Nghost+1,M
	         II = mb + ( I-1+(M-2*Nghost)*(px-1) )/isk
	         JJ = nb + ( J-1+(N-2*Nghost)*(K-1) )/isk
               rII = REAL(mb) + REAL(I-1+(M-2*Nghost)*(px-1))/REAL(isk) - REAL(II)
               rJJ = REAL(nb) + REAL(J-1+(N-2*Nghost)*(K-1))/REAL(isk) - REAL(JJ)

			 IN=I+Nghost-M
			 IF( myid == ParentLocalID(IN,J) ) THEN  !be careful
	           loci = II - parent_area_i*npx 
			   locj = JJ - parent_area_j*npy
           
			   xx(1) = Fin(loci, locj);
			   xx(2) = Fin(loci+1, locj);
			   xx(3) = Fin(loci, locj+1);
			   xx(4) = Fin(loci+1, locj+1);

			   CALL MPI_ISEND(xx,4,MPI_SP,ProcessorID(px,K),11,MPI_COMM_WORLD,ireq,ier)
			   CALL MPI_WAIT( ireq, istatus, ier )

			 ENDIF

			 IF( myid == ProcessorID(px,K) ) THEN
	           CALL MPI_IRECV(xx,4,MPI_SP,ParentLocalID(IN,J),11,MPI_COMM_WORLD,ireq,ier)
			   CALL MPI_WAIT( ireq, istatus, ier )
 
			   Fout_east(J,I-M+Nghost)  &
	           = ( (1.0_SP-rII)*xx(1) + rII*xx(2) )*( 1.0_SP-rJJ ) + &
                   ( (1.0_SP-rII)*xx(3) + rII*xx(4) )*rJJ
	
			 ENDIF
	      ENDDO
	   ENDDO

	ENDDO
	DEALLOCATE( ParentLocalID, xx )
	
!     if( n_east .eq. MPI_PROC_NULL ) then
!       DO J=1,N
!        DO I=M-Nghost+1,M
!         II = mb + (I-1+(M-2*Nghost)*npx)/isk
!         JJ = nb + (J-1+(N-2*Nghost)*npy)/isk
!         rII = REAL(mb) + REAL(I-1+(M-2*Nghost)*npx)/REAL(isk) - REAL(II)
!         rJJ = REAL(nb) + REAL(J-1+(N-2*Nghost)*npy)/REAL(isk) - REAL(JJ)
!         Fout_east(J,I-M+Nghost)  &
!	     = ( (1.0_SP-rII)*VarGrid1(II,JJ) + rII*VarGrid1(II+1,JJ) )*( 1.0_SP-rJJ ) + &
!             ( (1.0_SP-rII)*VarGrid1(II,JJ+1) + rII*VarGrid1(II+1,JJ+1) )*rJJ
!        ENDDO
!       ENDDO 
!     endif

# else
! not parallel
       DO J=1,N
        DO I=M-Nghost+1,M
         II = mb + (I-1)/isk
	   JJ = nb + (J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_east(J,I-M+Nghost)=((1.0_SP-rII)*Fin(II,JJ)+rII*Fin(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1)+rII*Fin(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO
# endif

! --- South
# if defined (PARALLEL)
	len = M*Nghost
	ALLOCATE( ParentLocalID(1:M,1:Nghost), xx(4) ); ParentLocalID=-99;  xx=0.0_SP

	DO K=1,px

	   IF( myid == ProcessorID(K,1) ) THEN
	     ParentLocalID(1:M,1:Nghost)=SouthParentID(1:M,1:Nghost,ng)
	     CALL MPI_BCAST(ParentLocalID,len,MPI_INTEGER,ProcessorID(K,1),MPI_COMM_WORLD,ier)

         ELSE
		 CALL MPI_BCAST(ParentLocalID,len,MPI_INTEGER,ProcessorID(K,1),MPI_COMM_WORLD,ier)
	   ENDIF

	   DO I=1,M
	      DO J=1,Nghost
	         II = mb + ( I-1+(M-2*Nghost)*(K-1) )/isk
			 JJ = nb + ( J-1 )/isk
			 rII = REAL(mb) + REAL( I-1+(M-2*Nghost)*(K-1) )/REAL(isk) - REAL(II)
			 rJJ = REAL(nb) + REAL( J-1 )/REAL(isk) - REAL(JJ)

			 IF( myid == ParentLocalID(I,J) ) THEN
	           loci = II - parent_area_i*npx
			   locj = JJ - parent_area_j*npy

			   xx(1) = Fin(loci, locj);
			   xx(2) = Fin(loci+1, locj);
			   xx(3) = Fin(loci, locj+1);
			   xx(4) = Fin(loci+1, locj+1);
	
			   CALL MPI_ISEND(xx,4,MPI_SP,ProcessorID(K,1),12,MPI_COMM_WORLD,ireq,ier)
			   CALL MPI_WAIT( ireq, istatus, ier )
	         ENDIF

	         IF( myid == ProcessorID(K,1) ) THEN
	           CALL MPI_IRECV(xx,4,MPI_SP,ParentLocalID(I,J),12,MPI_COMM_WORLD,ireq,ier)
			   CALL MPI_WAIT( ireq, istatus, ier )

			   Fout_south(I,J)  &
	           = ( (1.0_SP-rII)*xx(1) + rII*xx(2) )*( 1.0_SP-rJJ ) + &
	             ( (1.0_SP-rII)*xx(3) + rII*xx(4) )*rJJ
	         ENDIF

	      ENDDO
	   ENDDO

	ENDDO
	DEALLOCATE( ParentLocalID, xx )

!     if( n_suth .eq. MPI_PROC_NULL ) then
!       DO I=1,M
!        DO J=1,Nghost
!         II = mb+(I-1+(M-2*Nghost)*npx)/isk
!         JJ = nb+(J-1+(N-2*Nghost)*npy)/isk  ! actually npy=0 at south
!         rII = REAL(mb) + REAL(I-1+(M-2*Nghost)*npx)/REAL(isk) - REAL(II)
!         rJJ = REAL(nb) + REAL(J-1+(N-2*Nghost)*npy)/REAL(isk) - REAL(JJ)
!         Fout_south(I,J)  &
!	     = ( (1.0_SP-rII)*VarGrid1(II,JJ) + rII*VarGrid1(II+1,JJ) )*( 1.0_SP-rJJ ) + &
!             ( (1.0_SP-rII)*VarGrid1(II,JJ+1) + rII*VarGrid1(II+1,JJ+1) )*rJJ
!        ENDDO
!       ENDDO
!     endif

# else
! not parallel
       DO I=1,M
        DO J=1,Nghost
         II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_south(I,J)=((1.0_SP-rII)*Fin(II,JJ)+rII*Fin(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1)+rII*Fin(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO 
# endif

! --- North
# if defined (PARALLEL)
	len = M*Nghost
	ALLOCATE( ParentLocalID(1:M,1:Nghost), xx(4) ); ParentLocalID=-99;  xx=0.0_SP

	DO K=1,px

	   IF( myid == ProcessorID(K,py) ) THEN
	     ParentLocalID(1:M,1:Nghost)=NorthParentID(1:M,1:Nghost,ng)
	     CALL MPI_BCAST(ParentLocalID,len,MPI_INTEGER,ProcessorID(K,py),MPI_COMM_WORLD,ier)

         ELSE
		 CALL MPI_BCAST(ParentLocalID,len,MPI_INTEGER,ProcessorID(K,py),MPI_COMM_WORLD,ier)
	   ENDIF

	   DO I=1,M
	      DO J=N-Nghost+1,N
	         II = mb + ( I-1+(M-2*Nghost)*(K-1) )/isk
			 JJ = nb + ( J-1+(N-2*Nghost)*(py-1) )/isk
			 rII = REAL(mb) + REAL( I-1+(M-2*Nghost)*(K-1) )/REAL(isk) - REAL(II)
			 rJJ = REAL(nb) + REAL( J-1+(N-2*Nghost)*(py-1) )/REAL(isk) - REAL(JJ)

			 JN=J+Nghost-N
			 IF( myid == ParentLocalID(I,JN) ) THEN  !be careful
	           loci = II - parent_area_i*npx
			   locj = JJ - parent_area_j*npy

			   xx(1) = Fin(loci, locj);
			   xx(2) = Fin(loci+1, locj);
			   xx(3) = Fin(loci, locj+1);
			   xx(4) = Fin(loci+1, locj+1);
	
			   CALL MPI_ISEND(xx,4,MPI_SP,ProcessorID(K,py),13,MPI_COMM_WORLD,ireq,ier)
			   CALL MPI_WAIT( ireq, istatus, ier )
	         ENDIF

	         IF( myid == ProcessorID(K,py) ) THEN
	           CALL MPI_IRECV(xx,4,MPI_SP,ParentLocalID(I,JN),13,MPI_COMM_WORLD,ireq,ier)
			   CALL MPI_WAIT( ireq, istatus, ier )

			   Fout_north(I,J-N+Nghost)  &
	           = ( (1.0_SP-rII)*xx(1) + rII*xx(2) )*( 1.0_SP-rJJ ) + &
	             ( (1.0_SP-rII)*xx(3) + rII*xx(4) )*rJJ
	         ENDIF

	      ENDDO
	   ENDDO
	ENDDO
	DEALLOCATE( ParentLocalID, xx )

!     if( n_nrth .eq. MPI_PROC_NULL ) then
!       DO I=1,M
!        DO J=N-Nghost+1,N
!         II=mb+(I-1+(M-2*Nghost)*npx)/isk
!         JJ=nb+(J-1+(N-2*Nghost)*npy)/isk
!         rII=REAL(mb)+REAL(I-1+(M-2*Nghost)*npx)/REAL(isk)-REAL(II)
!         rJJ=REAL(nb)+REAL(J-1+(N-2*Nghost)*npy)/REAL(isk)-REAL(JJ)
!         Fout_north(I,J-N+Nghost)  &
!	     = ( (1.0_SP-rII)*VarGrid1(II,JJ) + rII*VarGrid1(II+1,JJ) )*( 1.0_SP-rJJ )+ &
!             ( (1.0_SP-rII)*VarGrid1(II,JJ+1) + rII*VarGrid1(II+1,JJ+1) )*rJJ
!        ENDDO
!       ENDDO 
!     endif

# else
! not parallel
       DO I=1,M
        DO J=N-Nghost+1,N
	   II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_north(I,J-N+Nghost)=((1.0_SP-rII)*Fin(II,JJ)+rII*Fin(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1)+rII*Fin(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO 
# endif

END SUBROUTINE INTERP_BC

!=====================================================================

! --- For AMR modeling !ykchoi
# if defined (PARALLEL)
! --------------------------------------------
!  Gather 2D variables from all processors
!  this gathering includes ghost cells
!  08/20/2013, fyshi
! --------------------------------------------
SUBROUTINE GATHER_GRID( phi_out, phi_in, &
           Mloc, Nloc, mm, nn, Nghost )
! mm and nn are global but include ghost cells

    USE PARAM
    USE GLOBAL, ONLY : nprocs, npx, npy, myid, ier
    IMPLICIT NONE

    integer,intent(in) :: Mloc, Nloc, mm, nn, Nghost
    real(SP),dimension(Mloc,Nloc),intent(in) :: phi_in
    real(SP),dimension(mm,nn),intent(out) :: phi_out
    integer,dimension(nprocs) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(Mloc,Nloc) :: xx
    real(SP),dimension(nprocs) :: xxx
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: iglob,jglob,len,n,l

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

    !scattering to every processors
    do j = 1,nn
    do i = 1,mm
       if (myid.eq.0) then
	    do l=1,nprocs
	       xxx(l) = phi_out(i,j)
		enddo
       endif
       call MPI_Scatter(xxx,1,MPI_SP,&
            phi_out(i,j),1,MPI_SP,0,MPI_COMM_WORLD,ier)
    enddo
    enddo

END SUBROUTINE GATHER_GRID
# endif

!=====================================================================

! --- For AMR modeling !ykchoi
SUBROUTINE USE_NESTING_BC(istep,ratio,ng)
    USE GLOBAL
    IMPLICIT NONE
    INTEGER :: istep,ng
    INTEGER :: ratio

    REAL(SP) :: tmpA1, tmpA2, tmpB1, tmpB2
    REAL(SP) :: bdry1, bdry2

! time
      !tmp1=REAL(istep)/REAL(ratio)
      !tmp2=1.0_SP - tmp1;

	tmpA1=REAL(istep-1)/REAL(ratio)
	tmpA2=1.0_SP - tmpA1;

	tmpB1=REAL(istep)/REAL(ratio)
	tmpB2=1.0_SP - tmpB1;

! west boundary
# if defined(PARALLEL)
     IF (n_west .eq. MPI_PROC_NULL) THEN
# endif

      DO J=1,Nloc
      DO I=1,Nghost
        !ETA(I,J)=Z_NESTING_WEST_PAR(J,I,2)*tmp1&
        !        +Z_NESTING_WEST_PAR(J,I,1)*tmp2
        !U(I,J)=U_NESTING_WEST_PAR(J,I,2)*tmp1&
        !        +U_NESTING_WEST_PAR(J,I,1)*tmp2
        !V(I,J)=V_NESTING_WEST_PAR(J,I,2)*tmp1&
        !        +V_NESTING_WEST_PAR(J,I,1)*tmp2

	  !ETA(I,J)=( Z_NESTING_WEST_PAR(J,I,1) + Z_NESTING_WEST_PAR(J,I,2) )*0.5_SP
	  !U(I,J)=( U_NESTING_WEST_PAR(J,I,1) + U_NESTING_WEST_PAR(J,I,2) )*0.5_SP
	  !V(I,J)=( V_NESTING_WEST_PAR(J,I,1) + V_NESTING_WEST_PAR(J,I,2) )*0.5_SP

        bdry1 = Z_NESTING_WEST_PAR(J,I,2)*tmpA1  &
	         +Z_NESTING_WEST_PAR(J,I,1)*tmpA2
        !bdry2 = Z_NESTING_WEST_PAR(J,I,2)*tmpB1  &
	  !       +Z_NESTING_WEST_PAR(J,I,1)*tmpB2
        ETA(I,J) = bdry1
        !ETA(I,J) = 0.5_SP*( bdry1 + bdry2 )
	  
        bdry1 = U_NESTING_WEST_PAR(J,I,2)*tmpA1  &
	         +U_NESTING_WEST_PAR(J,I,1)*tmpA2
        !bdry2 = U_NESTING_WEST_PAR(J,I,2)*tmpB1  &
	  !       +U_NESTING_WEST_PAR(J,I,1)*tmpB2
        U(I,J) = bdry1
        !U(I,J) = 0.5_SP*( bdry1 + bdry2 )

        bdry1 = V_NESTING_WEST_PAR(J,I,2)*tmpA1  &
	         +V_NESTING_WEST_PAR(J,I,1)*tmpA2
        !bdry2 = V_NESTING_WEST_PAR(J,I,2)*tmpB1  &
	  !       +V_NESTING_WEST_PAR(J,I,1)*tmpB2
        V(I,J) = bdry1
        !V(I,J) = 0.5_SP*( bdry1 + bdry2 )

        HU(I,J)=Max( Depth(I,J)+ETA(I,J)*Gamma3,MinDepthFrc )*U(I,J)
        HV(I,J)=Max( Depth(I,J)+ETA(I,J)*Gamma3,MinDepthFrc )*V(I,J)
      ENDDO
      ENDDO

# if defined (PARALLEL)
     ENDIF  ! end in domain
# endif

! east boundary
# if defined(PARALLEL)
     IF (n_east .eq. MPI_PROC_NULL) THEN
# endif

      DO J=1,Nloc 
      DO I=Iend+1,Iend+Nghost
        !ETA(I,J)=Z_NESTING_EAST_PAR(J,I-Iend,2)*tmp1&
        !        +Z_NESTING_EAST_PAR(J,I-Iend,1)*tmp2
        !U(I,J)=U_NESTING_EAST_PAR(J,I-Iend,2)*tmp1&
        !        +U_NESTING_EAST_PAR(J,I-Iend,1)*tmp2
        !V(I,J)=V_NESTING_EAST_PAR(J,I-Iend,2)*tmp1&
        !        +V_NESTING_EAST_PAR(J,I-Iend,1)*tmp2

	  !ETA(I,J)=( Z_NESTING_EAST_PAR(J,I-Iend,1) + Z_NESTING_EAST_PAR(J,I-Iend,2) )*0.5_SP
	  !U(I,J)=( U_NESTING_EAST_PAR(J,I-Iend,1) + U_NESTING_EAST_PAR(J,I-Iend,2) )*0.5_SP
	  !V(I,J)=( V_NESTING_EAST_PAR(J,I-Iend,1) + V_NESTING_EAST_PAR(J,I-Iend,2) )*0.5_SP

        bdry1 = Z_NESTING_EAST_PAR(J,I-Iend,2)*tmpA1  &
	         +Z_NESTING_EAST_PAR(J,I-Iend,1)*tmpA2
        !bdry2 = Z_NESTING_EAST_PAR(J,I-Iend,2)*tmpB1  &
	  !       +Z_NESTING_EAST_PAR(J,I-Iend,1)*tmpB2
        ETA(I,J) = bdry1
        !ETA(I,J) = 0.5_SP*( bdry1 + bdry2 )
	  
        bdry1 = U_NESTING_EAST_PAR(J,I-Iend,2)*tmpA1  &
	         +U_NESTING_EAST_PAR(J,I-Iend,1)*tmpA2
        !bdry2 = U_NESTING_EAST_PAR(J,I-Iend,2)*tmpB1  &
	  !       +U_NESTING_EAST_PAR(J,I-Iend,1)*tmpB2
        U(I,J) = bdry1
        !U(I,J) = 0.5_SP*( bdry1 + bdry2 )

        bdry1 = V_NESTING_EAST_PAR(J,I-Iend,2)*tmpA1  &
	         +V_NESTING_EAST_PAR(J,I-Iend,1)*tmpA2
        !bdry2 = V_NESTING_EAST_PAR(J,I-Iend,2)*tmpB1  &
	  !       +V_NESTING_EAST_PAR(J,I-Iend,1)*tmpB2
        V(I,J) = bdry1
        !V(I,J) = 0.5_SP*( bdry1 + bdry2 )

        HU(I,J)=Max( Depth(I,J)+ETA(I,J)*Gamma3,MinDepthFrc )*U(I,J)
        HV(I,J)=Max( Depth(I,J)+ETA(I,J)*Gamma3,MinDepthFrc )*V(I,J)
      ENDDO
      ENDDO

# if defined (PARALLEL)
     ENDIF  ! end in domain
# endif

! south boundary
# if defined(PARALLEL)
     IF (n_suth .eq. MPI_PROC_NULL) THEN
# endif

      DO I=1,Mloc
      DO J=1,Nghost
        !ETA(I,J)=Z_NESTING_SOUTH_PAR(I,J,2)*tmp1&
        !        +Z_NESTING_SOUTH_PAR(I,J,1)*tmp2
        !U(I,J)=U_NESTING_SOUTH_PAR(I,J,2)*tmp1&
        !        +U_NESTING_SOUTH_PAR(I,J,1)*tmp2
        !V(I,J)=V_NESTING_SOUTH_PAR(I,J,2)*tmp1&
        !        +V_NESTING_SOUTH_PAR(I,J,1)*tmp2
	  !ETA(I,J)=( Z_NESTING_SOUTH_PAR(I,J,1) + Z_NESTING_SOUTH_PAR(I,J,2) )*0.5_SP
	  !U(I,J)=( U_NESTING_SOUTH_PAR(I,J,1) + U_NESTING_SOUTH_PAR(I,J,2) )*0.5_SP
	  !V(I,J)=( V_NESTING_SOUTH_PAR(I,J,1) + V_NESTING_SOUTH_PAR(I,J,2) )*0.5_SP

        bdry1 = Z_NESTING_SOUTH_PAR(I,J,2)*tmpA1  &
	         +Z_NESTING_SOUTH_PAR(I,J,1)*tmpA2
        !bdry2 = Z_NESTING_SOUTH_PAR(I,J,2)*tmpB1  &
	  !       +Z_NESTING_SOUTH_PAR(I,J,1)*tmpB2
        ETA(I,J) = bdry1
        !ETA(I,J) = 0.5_SP*( bdry1 + bdry2 )
	  
        bdry1 = U_NESTING_SOUTH_PAR(I,J,2)*tmpA1  &
	         +U_NESTING_SOUTH_PAR(I,J,1)*tmpA2
        !bdry2 = U_NESTING_SOUTH_PAR(I,J,2)*tmpB1  &
	  !       +U_NESTING_SOUTH_PAR(I,J,1)*tmpB2
        U(I,J) = bdry1
        !U(I,J) = 0.5_SP*( bdry1 + bdry2 )

        bdry1 = V_NESTING_SOUTH_PAR(I,J,2)*tmpA1  &
	         +V_NESTING_SOUTH_PAR(I,J,1)*tmpA2
        !bdry2 = V_NESTING_SOUTH_PAR(I,J,2)*tmpB1  &
	  !       +V_NESTING_SOUTH_PAR(I,J,1)*tmpB2
        V(I,J) = bdry1
        !V(I,J) = 0.5_SP*( bdry1 + bdry2 )

        HU(I,J)=Max( Depth(I,J)+ETA(I,J)*Gamma3,MinDepthFrc )*U(I,J)
        HV(I,J)=Max( Depth(I,J)+ETA(I,J)*Gamma3,MinDepthFrc )*V(I,J)
      ENDDO
      ENDDO 

# if defined (PARALLEL)
     ENDIF  ! end in domain
# endif

! north boundary
# if defined(PARALLEL)
     IF (n_nrth .eq. MPI_PROC_NULL) THEN
# endif

      DO I=1,Mloc
      DO J=Jend+1,Jend+Nghost
        !ETA(I,J)=Z_NESTING_NORTH_PAR(I,J-Jend,2)*tmp1&
        !        +Z_NESTING_NORTH_PAR(I,J-Jend,1)*tmp2
        !U(I,J)=U_NESTING_NORTH_PAR(I,J-Jend,2)*tmp1&
        !        +U_NESTING_NORTH_PAR(I,J-Jend,1)*tmp2
        !V(I,J)=V_NESTING_NORTH_PAR(I,J-Jend,2)*tmp1&
        !        +V_NESTING_NORTH_PAR(I,J-Jend,1)*tmp2
	  !ETA(I,J)=( Z_NESTING_NORTH_PAR(I,J-Jend,1) + Z_NESTING_NORTH_PAR(I,J-Jend,2) )*0.5_SP
	  !U(I,J)=( U_NESTING_NORTH_PAR(I,J-Jend,1) + U_NESTING_NORTH_PAR(I,J-Jend,2) )*0.5_SP
	  !V(I,J)=( V_NESTING_NORTH_PAR(I,J-Jend,1) + V_NESTING_NORTH_PAR(I,J-Jend,2) )*0.5_SP

        bdry1 = Z_NESTING_NORTH_PAR(I,J-Jend,2)*tmpA1  &
	         +Z_NESTING_NORTH_PAR(I,J-Jend,1)*tmpA2
        !bdry2 = Z_NESTING_NORTH_PAR(I,J-Jend,2)*tmpB1  &
	  !       +Z_NESTING_NORTH_PAR(I,J-Jend,1)*tmpB2
        ETA(I,J) = bdry1
        !ETA(I,J) = 0.5_SP*( bdry1 + bdry2 )
	  
        bdry1 = U_NESTING_NORTH_PAR(I,J-Jend,2)*tmpA1  &
	         +U_NESTING_NORTH_PAR(I,J-Jend,1)*tmpA2
        !bdry2 = U_NESTING_NORTH_PAR(I,J-Jend,2)*tmpB1  &
	  !       +U_NESTING_NORTH_PAR(I,J-Jend,1)*tmpB2
        U(I,J) = bdry1
        !U(I,J) = 0.5_SP*( bdry1 + bdry2 )

        bdry1 = V_NESTING_NORTH_PAR(I,J-Jend,2)*tmpA1  &
	         +V_NESTING_NORTH_PAR(I,J-Jend,1)*tmpA2
        !bdry2 = V_NESTING_NORTH_PAR(I,J-Jend,2)*tmpB1  &
	  !       +V_NESTING_NORTH_PAR(I,J-Jend,1)*tmpB2
        V(I,J) = bdry1
        !V(I,J) = 0.5_SP*( bdry1 + bdry2 )

        HU(I,J)=Max( Depth(I,J)+ETA(I,J)*Gamma3,MinDepthFrc )*U(I,J)
        HV(I,J)=Max( Depth(I,J)+ETA(I,J)*Gamma3,MinDepthFrc )*V(I,J)
      ENDDO
      ENDDO 

# if defined (PARALLEL)
     ENDIF  ! end in domain
# endif

END SUBROUTINE USE_NESTING_BC 