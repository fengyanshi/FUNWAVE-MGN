! --- For AMR modeling !ykchoi
MODULE LOCAL
       USE PARAM
       IMPLICIT NONE
       SAVE

! grid
# if defined(CARTESIAN)
         REAL(SP),DIMENSION(:),ALLOCATABLE :: DX_TILE,DY_TILE
# else
         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: DX_TILE,DY_TILE
# endif

! bathy
         REAL(SP),DIMENSION(:,:,:), ALLOCATABLE :: Depth_Tile, &
                 DepthX_Tile, DepthY_Tile
         INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: MASK_Tile,MASK9_Tile
# if defined(CARTESIAN)
!   dont need latitude
# else
         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: Lat_Theta_Tile
# endif

! variable
         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: U_Tile,V_Tile,Eta_Tile  
         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: Ubar_Tile,Vbar_Tile       

       CONTAINS

SUBROUTINE ALLOCATE_VAR_TILE
       USE PARAM
# if defined(PARALLEL)
       USE GLOBAL, ONLY: K,NumGrid,GridDimX,GridDimY,Nghost,ier,px,py,  &
                         MaxDimX,MaxDimX1,MaxDimY,MaxDimY1      
# else
       USE GLOBAL, ONLY: K,NumGrid,GridDimX,GridDimY,Nghost, &
                         MaxDimX,MaxDimX1,MaxDimY,MaxDimY1 
# endif
       IMPLICIT NONE

# if defined(PARALLEL)
       MaxDimX = MAXVAL(GridDimX)/px+2*Nghost
       MaxDimY = MAXVAL(GridDimY)/py+2*Nghost
# else
       MaxDimX = MAXVAL(GridDimX)+2*Nghost
       MaxDimY = MAXVAL(GridDimY)+2*Nghost
# endif
       MaxDimX1 = MaxDimX+1
       MaxDimY1 = MaxDimY+1
       
! static
# if defined(CARTESIAN)
       ALLOCATE ( DX_TILE(NumGrid) )
       ALLOCATE ( DY_TILE(NumGrid) )
# else
       ALLOCATE ( DX_TILE(MaxDimX,MaxDimY,NumGrid) )
       ALLOCATE ( DY_TILE(MaxDimX,MaxDimY,NumGrid) )
# endif
       ALLOCATE ( Depth_Tile(MaxDimX,MaxDimY,NumGrid) )
       ALLOCATE ( DepthX_Tile(MaxDimX1,MaxDimY,NumGrid) )
       ALLOCATE ( DepthY_Tile(MaxDimX,MaxDimY1,NumGrid) )
       ALLOCATE ( MASK_Tile(MaxDimX,MaxDimY,NumGrid) )
       ALLOCATE ( MASK9_Tile(MaxDimX,MaxDimY,NumGrid) )
# if defined(CARTESIAN)
!   no latitude for cartesian
# else
       ALLOCATE ( Lat_theta_Tile(MaxDimX,MaxDimY,NumGrid))
# endif
! dynamic
       ALLOCATE ( U_Tile(MaxDimX,MaxDimY,NumGrid))
       ALLOCATE ( V_Tile(MaxDimX,MaxDimY,NumGrid))
       ALLOCATE ( Ubar_Tile(MaxDimX,MaxDimY,NumGrid))
       ALLOCATE ( Vbar_Tile(MaxDimX,MaxDimY,NumGrid))

       ALLOCATE ( Eta_Tile(MaxDimX,MaxDimY,NumGrid)) 

! initialization
        MASK_TILE = 1
        MASK9_TILE = 1
        U_TILE = ZERO
        V_TILE = ZERO
        Ubar_TILE = ZERO
        Vbar_TILE = ZERO
        Eta_TILE = ZERO

END SUBROUTINE ALLOCATE_VAR_TILE

END MODULE LOCAL 