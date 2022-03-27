!  mod2005_binary.f90 
!
!  FUNCTIONS:
!  mod2005_binary - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: mod2005_binary
!
!  PURPOSE:  read binary file of cell-by-cell flow of mf2005 mf2k
!
!****************************************************************************

    program mod2005_binary
    use variable_def
    use geometry_class
    use functions
    use advection_class
    implicit none
    


    
    ! Variables
    character(len=120)  ::  fname     ! file name
    integer             ::  iprec, iunit, ncol, nrow, nlay, iufdgb
    type(mf2k_cl)       :: mf2k
    type(geometry_cl)   :: geo
    type(advection_cl)  :: advection
    
    call allocate_mf2k_ (mf2k)
    !fname = "./small_test_mf.cbb"
    fname = "./mf2005.cbc"
    ! check file 
    call open_fname (trim(fname), iunit, 'unknown', 'sequential', 'binary')
    call budget_precision (iprec,iunit,ncol,nrow,nlay)
    ! check file
    if(iprec.lt.1) then
         close (iunit)
         open(unit=iunit,file=trim(fname),status='old',form='unformatted',access='sequential',err=10)   !      
         call budget_precision(iprec,iunit,ncol,nrow,nlay)
    end if
    
    ! stop iprec !=1 or 2
    if(iprec.lt.1) then
        write(*,*) 'stopping because budget file is invalid'
        stop
    ! 1 single
    elseif(iprec.eq.1) then
        write(*,*) ' single precision budget file'
    ! 2 double
    else if(iprec.eq.2) then
        write(*,*) ' double precision budget file'
    end if
    
    ! print
    write(*,*)
    write(*,14) nlay,nrow,ncol
    !write(*,14) nlay,nrow,ncol
    write(*,*)
14  format(1x,i10,' layers',i10,' rows',i10,' columns')
     
    mf2k%preCBC = iprec
    mf2k%nx = ncol
    mf2k%ny = nrow
    mf2k%nz = nlay
    mf2k%unitCBC = iunit
    mf2k%file    = fname
      
    call read_flux_from_mf2k_ (mf2k,advection,geo) ! read from cell-by-cell budget modflow
      
10    stop 'Could not open budget file from modflow' ! 打开文件出错到此执行
      
      
    end program mod2005_binary
    
    

 
!     ******************************************************************
!     Determine single or double precision file type for a MODFLOW
!     budget file:  0=unrecognized, 1=single, 2=double.
!     ******************************************************************
    subroutine budget_precision (iprec,iu,ncol,nrow,nlay)
      use functions, only:  upper_case  ! 函数声明
      implicit none
      
      integer, intent(in)    :: iu
      integer, intent(out)   :: iprec,ncol,nrow,nlay
      real*8                 :: deltd,pertimd,totimd,vald ! double
      real                   :: delt,pertim,totim,val     ! single
      character(len=16)      :: text1,text2
      integer                :: kstp,kper,icode,nodes ! icode 0 1 2三种情况
      integer                :: nlst,n,icell,nc,nr,nl
      real, pointer          :: buff(:,:,:) ! single 1
      real*8, pointer        :: buffd(:,:,:) ! double precision 2
!
!  default is unrecognized file
     
      iprec=0  ! 0 1 2
!
!  single check
!
      ! KSTP time step number, interger 4bytes
      ! KPER stress period number, interger 4 bytes 
      ! DESC, a description of the array, 16 ANSI characters, 16 bytes.
      ! NCOL, the number of columns in the array, an integer, 4 bytes
      ! NROW, the number of rows in the array, an integer, 4 bytes.
      ! NLAY, the number of layers in the array, an integer, 4 bytes. 
      ! NLAY can be a negative number in which case the absolute value of NLAY 
      ! is the number of layers in the array. If NLAY is greater than zero, 
      ! a 3D array of real numbers follows NLAY.  The number of values is NCOL x NROW x NLAY. 
      ! To read it, you can use a loop over the layers that contains a loop over rows 
      ! that contains a loop over columns.
      ! 
      read(iu,err=100,end=100) kstp,kper,text1,ncol,nrow,nlay ! integer character
      
      text1 = upper_case (text1)
      
      icode=0
      
      ! If NLAY is less than zero, the compact format is being used.
      ! With the compact format, read the following:
      ! ITYPE, a value indicating how the data is stored, an integer, 4 bytes. icode here
      ! DELT, the length of the current time step, a real number, either 4 or 8 bytes.
      ! PERTIM: the time in the current stress period, a real number, either 4 or 8 bytes.
      ! TOTIM, the total elapsed time, a real number, either 4 or 8 bytes.
      if(nlay.lt.0) then
        nlay=-nlay ! change to a positvie value
        read(iu,err=50,end=50) icode,delt,pertim,totim  ! single
      end if
      
14    format(1x,i10,' layers',i10,' rows',i10,' columns')
      
      ! wrong data, rewind file and return
      if(ncol.lt.1 .or. nrow.lt.1 .or. nlay.lt.1) go to 100  ! wrong dimension
      if(ncol.gt.100000000 .or.nrow.gt.100000000 .or. nlay.gt.100000000) go to 100  ! wrong dimension
      if(ncol*nrow.gt.100000000 .or. ncol*nlay.gt.100000000 .or. nrow*nlay.gt.100000000) go to 100
      
      allocate (buff(ncol,nrow,nlay))
      allocate (buffd(ncol,nrow,nlay))
      nodes=ncol*nrow*nlay
!
!  read data depending on icode.  icode 0,1, or 2 are the only allowed
!  values because the first budget terms must be from the internal
!  flow package (bcf,lpf, or huf).
!
      ! If ITYPE = 0 or 1, read a 3D array of values. 
      ! Each value is either 4 or 8 bytes depending on whether the file 
      ! is saved with single- or double-precision data.
      if(icode.eq.0 .or. icode.eq.1) then
         read(iu,err=50,end=50) buff   ! single precision
         
      ! If ITYPE = 2 read NLIST which is the number of cells for which values will be stored.
      ! NLIST is an integer, 4 bytes.
      else if(icode.eq.2) then
         read(iu,err=50,end=50) nlst  ! NVAL or NLIST
         ! Next we read the values associated with the cells
         ! check whether NLIST is greater than zero
         if(nlst.lt.0) go to 50   ! no data goto 50 重新读取
         ! ! If it is greater than zero, compute NRC as NROW x NCOL.
         ! then in a loop from 1 to NLIST read first ICELL (an integer, 4 bytes) 
         ! and then NVAL values (real numbers either 4 or 8 bytes each).
         ! 
         if(nlst.gt.0) then       !  read a list of cells and their associated values
            do n=1,nlst
               read(iu,end=50,err=50) icell,val  ! single precision
               if(icell.le.0 .or. icell.gt.nodes) go to 50  ! start over with KSTP
            end do
         end if
      else
         go to 100
      end if
!
!  read 2nd header and check for valid type. second flow terms
!
      read(iu,err=50,end=50) kstp,kper,text2
      text2 = upper_case (text2)
      if(text1.eq.'         STORAGE' .and. text2.eq.'   CONSTANT HEAD') then
           iprec=1  ! single precision 
           go to 100
      else if(text1.eq.'   CONSTANT HEAD' .and. text2.eq.'FLOW RIGHT FACE ') then
           iprec=1
           go to 100
      end if
!
!  double check  双精度 检查以 双精度 读取 再检查text是否正确
!
50    rewind(iu)
      read(iu,err=100,end=100) kstp,kper,text1,nc,nr,nl
      icode=0
      if(nl.lt.0) then
        nl=-nl
        read(iu,err=100,end=100) icode,deltd,pertimd,totimd
      end if
!
!  read data depending on icode.  icode 0,1, or 2 are the only allowed
!  values because the first budget terms must be from the internal
!  flow package (bcf,lpf, or huf).
!
      if(icode.eq.0 .or. icode.eq.1) then
         read(iu,err=100,end=100) buffd
      else if(icode.eq.2) then
         read(iu,err=100,end=100) nlst
         if(nlst.lt.0) go to 100
         if(nlst.gt.0) then
            do n=1,nlst
               read(iu,end=100,err=100) icell,vald
               if(icell.le.0 .or. icell.gt.nodes) go to 100
            end do
         end if
      else
         go to 100
      end if
!
!  read 2nd header and check for valid type.
!
      read(iu,err=100,end=100) kstp,kper,text2
      if(text1.eq.'         storage' .and.  text2.eq.'   constant head') then
           iprec=2
      else if(text1.eq.'   constant head' .and. text2.eq.'flow right face ') then
           iprec=2
      end if
!
100   rewind(iu)  ! 回到文件头部
      if(associated(buff))  deallocate (buff)
      if(associated(buffd)) deallocate (buffd)
      return
      
    end subroutine
    
    
    
    subroutine allocate_mf2k_ (mf2k)
      use variable_def
 	  implicit none
      type(mf2k_cl),    intent(inout) :: mf2k
		 allocate (mf2k%nx)
		 allocate (mf2k%ny)
		 allocate (mf2k%nz)          
	     allocate (mf2k%preCBC) 
	     allocate (mf2k%unitCBC)
	     allocate (mf2k%formCBC)
	     allocate (mf2k%file)
		 mf2k%nx = 0
		 mf2k%ny = 0
		 mf2k%nz = 0          
	     mf2k%preCBC = 0 
	     mf2k%unitCBC = 1
	     mf2k%formCBC ='UNFORMATTED'
	     mf2k%file = ' '
    end subroutine   
    
    


