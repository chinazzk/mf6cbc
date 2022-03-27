!  mod6_binary.f90 
!
!  FUNCTIONS:
!  mod6_binary - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: mod6_binary
!
!  PURPOSE:  Read cell-by-cell flow in MODFLOW6 binary budget file
!
!****************************************************************************

    program mod6_binary
    use variable_def ! mf6_cl
    use geometry_class
    use functions
    use advection_class
    implicit none
    
    ! Variables
    character(len=120)  :: fname_cbc, fname_dis ! binary budget file and binary grid file
    integer             :: iunit, ncol, nrow, nlay
    type(mf6_cl)       :: mf6
    type(geometry_cl)   :: geo
    type(advection_cl)  :: advection

    
    call allocate_mf6_ (mf6)
    fname_cbc = "./mf6.cbb"
    fname_dis = "./small_test_mf.dis.grb"
    
    ! read gird file, mf6 file is binary not unformated
    call open_fname (trim(fname_dis), iunit, 'unknown', 'sequential', 'binary')
    call read_mf6dis (iunit,ncol,nrow,nlay,mf6)
    ! binary cell by cell flow file
    open(iunit,file=trim(fname_cbc),status='old',form='binary',access='sequential',err=10) 
    
    ! print to screen
    write(*,*)
    write(*,14) nlay,nrow,ncol
    write(*,*)
14  format(1x,i10,' layers',i10,' rows',i10,' columns')
     
    mf6%nx = ncol
    mf6%ny = nrow
    mf6%nz = nlay
    mf6%unitCBC = iunit
    mf6%file    = fname_cbc

      
    call read_flux_from_mf6_ (mf6,advection,geo) ! read from cell-by-cell budget modflow
      
10    stop 'Could not open budget file from modflow' ! 
      
    end program mod6_binary
    
    

 
!******************************************************************
!     Read binary dis file for a MODFLOW6
!     only used for MODFLOW6
!    ja : a list of cell number (n) followed by its connecting cell numbers (m) for each of the m cells connected
!         to cell n. The number of values to provide for cell n is IAC(n). This list is sequentially provided
!         for the first to the last cell. The first value in the list must be cell n itself, and the remaining
!         cells must be listed in an increasing order (sorted from lowest number to highest). 
!    nja: When calculating the total number of connections,the connection between cell n and cell m is considered 
!          to be different from the connection between cell m and cell n. Thus, NJA is equal to the total number 
!           of connections, including n to m and m to n, and the total number of cells.
!    ia : is the number of connections (plus 1) for each cell. The sum of all the entries in IAC must be
!         equal to NJA.
!******************************************************************
    subroutine read_mf6dis (iu,ncol,nrow,nlay,mf6)
      use functions, only:  upper_case  
      use variable_def
      implicit none
      
      integer, intent(in)        :: iu
      integer, intent(out)       :: ncol,nrow,nlay
      type(mf6_cl),intent(inout) :: mf6
      character(len=50)          :: head1
      character(len=100)         :: txt
      integer                    :: nja,n,ncells
      real*8                     :: x_ori, y_ori, ang_rot
      real*8, pointer            :: delr(:),delc(:),top(:),botm(:)
      integer, pointer           :: idomain(:),icelltype(:)

!read 4 head of binary dis file
      do n=1,4
          read(iu,err=100,end=100) head1
      end do
!read 16 strings of binary dis file
      do n=1,16
          read(iu,err=100,end=100) txt
      end do
!read integer dis data      
      read(iu,err=100,end=100) ncells ! 读到文件尾部 跳到101
      read(iu,err=100,end=100) nlay
      read(iu,err=100,end=100) nrow
      read(iu,err=100,end=100) ncol
      read(iu,err=100,end=100) nja ! 4960
!read double dis data      
      read(iu,err=100,end=100) x_ori
      read(iu,err=100,end=100) y_ori
      read(iu,err=100,end=100) ang_rot
!read double array
      allocate(delr(ncol))
      allocate(delc(nrow))
      allocate(top(nrow*ncol))
      allocate(botm(ncells))
      read(iu,err=100,end=100) delr  ! row direction
      read(iu,err=100,end=100) delc
      read(iu,err=100,end=100) top
      read(iu,err=100,end=100) botm
!read integer array
      allocate(mf6%ia(ncells+1))
      allocate(mf6%ja(nja))
      allocate(idomain(ncells))
      allocate(icelltype(ncells))
      read(iu,err=100,end=100) mf6%ia   ! the number of connections (plus 1)for each cell
      read(iu,err=100,end=100) mf6%ja   ! a list of cell num  
      read(iu,err=100,end=100) idomain
      read(iu,err=100,end=100) icelltype
      
100   rewind(iu)  ! back to begin
      if(associated(delr))  deallocate (delr)
      if(associated(delc)) deallocate (delc)
      if(associated(top)) deallocate (top)
      if(associated(botm)) deallocate (botm)
      if(associated(icelltype)) deallocate (icelltype)
      if(associated(idomain)) deallocate (idomain)
      
      close(iu) ! close binary file
      return
      
    end subroutine
    
!*******************************************************    
!
!*******************************************************    
    subroutine allocate_mf6_ (mf6)
      use variable_def
 	  implicit none
      type(mf6_cl),    intent(inout) :: mf6
		 allocate (mf6%nx)
		 allocate (mf6%ny)
		 allocate (mf6%nz)          
	     !allocate (mf6%preCBC) 
	     allocate (mf6%unitCBC)
	     !allocate (mf6%formCBC)
	     allocate (mf6%file)
         !allocate (mf6%ncells)
		 mf6%nx = 0
		 mf6%ny = 0
		 mf6%nz = 0 
         !mf6%ncells = 0 ! add by zzk
	     !mf6%preCBC = 0 
	     mf6%unitCBC = 1
	     !mf6%formCBC ='UNFORMATTED'
	     mf6%file = ' '
    end subroutine   
    
    


