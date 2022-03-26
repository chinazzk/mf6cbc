    module variable_def
    use array_class
    implicit none
    
        type mf2k_cl            
		 character(len=50), pointer  :: file => null()
		 integer, pointer            :: nx => null()
		 integer, pointer            :: ny => null()
		 integer, pointer            :: nz => null()         
	     integer, pointer            :: preCBC  => null() 
	     integer, pointer            :: unitCBC => null()
	     character(len=12), pointer  :: formCBC => null()
    end type
    
    type advection_cl
      logical           :: action
	  type(array_cl)    :: qx,qy,qz,poro
	  integer           :: nt = 1              !total number of velocity snapshots
	  integer           :: it = 1              !timeshot interval associated with qx,qy,qz
	  integer           :: kper = 1            !stress period index
	  integer           :: kstp = 1            !timestep index (inside stress period)                 
	  real*8, pointer   :: time(:) => null()   !time discretization associated with velocity snapshots	  
      logical           :: switch = .TRUE.     !flag to switch velocity
    end type

    
    end module