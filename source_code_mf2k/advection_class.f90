    module advection_class
    
    implicit none
    
    
    contains
    
    
    
!********************************************************************************************
!   SUBROUTINE TO READ FACE FLUXES FROM MODFLOW
!********************************************************************************************
! This program reads the output binary file from modflow cell-by-cell flow terms
! and calculates the darcy velocity field.
!
! - qx,qy,qz are darcy velocities of grid interfaces
!
! - watch-out: axes from modflow are not equal to axes in RW3D
!
! - to go from (i,j,k) to modflow coordinates (jmod,imod,kmod):
!
!     jmod = i
!     imod = nrow-j+1
!     kmod = nlay-k+1
!
!*********************************************************************************************
  subroutine read_flux_from_mf2k_ (mf2k,advection,geo)
	 use array_class
	 use geometry_class
	 use functions,only: generate_unit,open_fname, upper_case
	 !use global_variables, only: fdbg
	 !use loops_particles, only: nmove
	 use variable_def
	 implicit none
	 type(mf2k_cl),      intent(in)    :: mf2k
	 type(geometry_cl),  intent(in)    :: geo
	 type(advection_cl), intent(inout) :: advection
     character(len=16)                 :: text,ctmp
	 real*8                            :: dx,dy,dz
     integer*4                         :: kstp,kper,ncol,nrow,nlay
     integer                           :: iunit,ierror,ioerr,imod,jmod,kmod,i,j,k,iufdbg,alloc_err
	 logical                           :: exists
	 real*8                            :: Qsto
	 real*8                            :: dzero
	 integer                           :: nc,nr,nl
	 integer                           :: inbud,iprec
     real*8                            :: deltd,pertimd,totimd,vald(20)
     real                              :: delt,pertim,totim,val(20)
     character(len=16)                 :: text1,text2
     integer                           :: icode,nodes
     integer                           :: nlist,n,icell
     integer                           :: itype,nval
     integer, save, pointer            :: ibuff(:,:,:)
     real,    save, pointer            :: buff(:,:,:)
     real*8,  save, pointer            :: buffd(:,:,:)
     logical                           :: KeepReading,assignflow

!....initialize

     dzero = 0.d0

     kper = advection%kper
     kstp = advection%kstp
     
     ncol = mf2k%nx
     nrow = mf2k%ny
     nlay = mf2k%nz
     
     if (.not.associated(buff))  allocate (buff(ncol,nrow,nlay))
     if (.not.associated(buffd)) allocate (buffd(ncol,nrow,nlay))
     if (.not.associated(ibuff)) allocate (ibuff(ncol,nrow,nlay))
     if (.not.associated(advection%qx%values)) advection%qx = make_array_ (0.d0,ncol+1,nrow,  nlay  )
	 if (.not.associated(advection%qy%values)) advection%qy = make_array_ (0.d0,ncol,  nrow+1,nlay  )
	 if (.not.associated(advection%qz%values)) advection%qz = make_array_ (0.d0,ncol,  nrow,  nlay+1)
     
     inbud = mf2k%unitCBC
     iprec = mf2k%preCBC

     assignflow = .FALSE.

!....read one time step


    loop: do

      read(inbud,end=1000,err=1000) kstp,kper,text,nc,nr,nl
      
      ! check if right stress period and time step is correct
      ! if still not there, exit and continue with the same velocity field

      if (kper > advection%kper  .or. (kper == advection%kper .and. kstp > advection%kstp) ) then
          backspace(inbud) ! 返回上一行
          exit
      end if

      if (kper == advection%kper  .and. kstp == advection%kstp ) assignflow =.TRUE. ! 

      ! start reading one internal flow term
      
      itype=0
      if(nl.lt.0) then  ! Nlay < 0 compact format 
         if(iprec.eq.1) then
           read(inbud) itype,delt,pertim,totim  ! single 4
         else
           read(inbud) itype,deltd,pertimd,totimd ! double 8
         end if
         nval=1 ! default value
         if(itype.eq.5) then
            read(inbud) nval
            if(nval.gt.1) then  ! If NVAL > 1, read (NVAL - 1) copies of CTMP. 
               do n=2,nval
                 read(inbud) ctmp ! CTMP is 16 ANSI characters, 16 bytes
               end do
            end if
         end if
         ! f ITYPE = 2 or 5 read NLIST which is the number of cells for which values will be stored
         if(itype.eq. 2 .or. itype.eq.5) read(inbud) nlist ! integer, 4 bytes.
      end if
      ! below for nlay > 0 
!------------------------------------------------------------------------------------------------------------------------------------------------------
!   read the budget term data under the following conditions:
!------------------------------------------------------------------------------------------------------------------------------------------------------
! if itype = 0 or 1, read a 3d array of values.
! if itype = 2 or 5, read a list of cells and their associated values.
! if itype = 3, read a 2d layer indicator array followed by a 2d array of ! values to be assigned to the layer indicated by the layer indicator array.
! if itype = 4 read a 2d array of values associated with layer 1.
!------------------------------------------------------------------------------------------------------------------------------------------------------
      if(itype.eq.0 .or. itype.eq.1) then ! If ITYPE = 0 or 1, read a 3D array of values.  
!  full 3-d array
         if(iprec.eq.1) then
           read(inbud) buff  ! single
           !buffd=buff
         else
           read(inbud) buffd ! double
         end if
      else if(itype.eq.3) then ! If ITYPE = 3, read a 2D layer indicator array followed by a 2D array of value
!  1-layer array with layer indicator array
         buffd=dzero
         read(inbud) ((ibuff(j,i,1),j=1,ncol),i=1,nrow) !  first read NROW*NCOL integer values (4 bytes each) row major order
         if(iprec.eq.1) then ! single
           read(inbud) ((buff(j,i,1),j=1,ncol),i=1,nrow) ! Then read NROW*NCOL real-number values (either 4 or 8 bytes each)
           do 265 i=1,nrow
           do 265 j=1,ncol
           buffd(j,i,1)=buff(j,i,1) ! 将该值赋值给 buffd
265        continue
         else ! double
           read(inbud) ((buffd(j,i,1),j=1,ncol),i=1,nrow)
         end if
         do 270 i=1,nrow
         do 270 j=1,ncol
         if(ibuff(j,i,1).ne.1) then  ! ????
            buffd(j,i,ibuff(j,i,1))=buffd(j,i,1)
            buffd(j,i,1)=dzero
         end if
270      continue
      else if(itype.eq.4) then ! read NROW*NCOL real-number values (either 4 or 8 bytes each) in row major order.
!  1-layer array that defines layer 1
         if(iprec.eq.1) then ! single
           read(inbud) ((buff(j,i,1),j=1,ncol),i=1,nrow)
           do 275 i=1,nrow
           do 275 j=1,ncol
           buffd(j,i,1)=buff(j,i,1)
275        continue
         else ! double
           read(inbud) ((buffd(j,i,1),j=1,ncol),i=1,nrow)
         end if
         if(nlay.gt.1) then ! 其他层赋值0？？？
            do 280 k=2,nlay
            do 280 i=1,nrow
            do 280 j=1,ncol
            buffd(j,i,k)=dzero
280         continue
         end if
      else if(nlist.gt.0) then  ! itype= 2 or 5 If ITYPE = 2 or 5 check whether NLIST is greater than zero.
!  list -- read only if the values need to be skipped.
         do 300 n=1,nlist !  then NVAL values (real numbers either 4 or 8 bytes each).
         if(iprec.eq.1) then ! ead NVAL which is the number of values associated with each cell.
           read(inbud) icell,(val(i),i=1,nval) ! in a loop from 1 to NLIST read first ICELL 
!            k= (icell-1)/nrc + 1
!            i= ( (icell - (k-1)*nrc)-1 )/ncol +1
!            j= icell - (k-1)*nrc - (i-1)*ncol
         else
           read(inbud) icell,(vald(i),i=1,nval)
!            k= (icell-1)/nrc + 1
!            i= ( (icell - (k-1)*nrc)-1 )/ncol +1
!            j= icell - (k-1)*nrc - (i-1)*ncol
         end if
300      continue
      end if
!
!-----------------------------------------------------------------------
!     associate internal flow terms to darcy fluxes
!-----------------------------------------------------------------------      
      if (assignflow) then
      
      text = upper_case (text)
      
      if(text.eq.'   CONSTANT HEAD') then
      
         ! do nothing
      
      elseif(text.eq.'FLOW RIGHT FACE ') then
            
            ! assign buffer to qx

            advection%qx%values = 0.d0

            do kmod=1,nlay      ! 4层rw3d
	         do imod=1,nrow     ! 10行rw3d
	           do jmod=1,ncol-1 ! 20列rw3d 最后一列没有 qx？qx=0
                   i=jmod        !    即jmod 1-19的buff赋值给 jmod+1的qx   2-20 对应其第一列没有为2-20 qx[1,:,:]和qx[21,:,:]为零任何情况下？
	               j=nrow-imod+1 !buff, 1-10, qx 10-1
                   k=nlay-kmod+1 !z方向相反       4-1
	               !dy = value_array_ (geo%dy,1,j,1)
	               !dz = value_array_ (geo%dz,1,1,k)   ! buff[cols20,rows10,layers4]
                   dy=1
                   dz=1
	               if (iprec == 1) advection%qx%values(i+1,j,k) = buff(jmod,imod,kmod)/dy/dz
	               if (iprec /= 1) advection%qx%values(i+1,j,k) = buffd(jmod,imod,kmod)/dy/dz
	           end do
	         end do
	        end do     
      
      elseif(text.eq.'FLOW FRONT FACE ') then

            ! assign buffer to qy
     
            advection%qy%values = 0.d0
                        
            do kmod=1,nlay     ! 4
	         do imod=1,nrow    ! 10     1-10
	          do jmod=1,ncol   ! 20
                   i=jmod        ! 1-20
	               j=nrow-imod   ! 9-0,   10-1     qy(:,11,:)=0     
                   k=nlay-kmod+1  ! 4-1
	               !dx = value_array_ (geo%dx,i,1,1)
	               !dz = value_array_ (geo%dz,1,1,k)
                   dx = 1
                   dz = 1
	               if (iprec == 1) advection%qy%values(i,j+1,k) = - buff(jmod,imod,kmod)/dx/dz
	               if (iprec /= 1) advection%qy%values(i,j+1,k) = - buffd(jmod,imod,kmod)/dx/dz
	          end do
	         end do
	        end do

      
      elseif(text.eq.'FLOW LOWER FACE ') then

            ! assign buffer to qz

            advection%qz%values = 0.d0

            do kmod=1,nlay-1    ! 1-3   最后一层buff忽略 对应qz的第一层忽略
	         do imod=1,nrow     ! 10
	          do jmod=1,ncol    ! 20
                   i=jmod       ! 1-20
	               j=nrow-imod+1  ! 10-1
                   k=nlay-kmod    !4-2
	               !dx = value_array_ (geo%dx,i,1,1)
	               !dy = value_array_ (geo%dy,1,j,1)
                   dx = 1
                   dy = 1
	               if (iprec == 1) advection%qz%values(i,j,k+1) = - buff(jmod,imod,kmod)/dx/dy
	               if (iprec /= 1) advection%qz%values(i,j,k+1) = - buffd(jmod,imod,kmod)/dx/dy
	          end do
	         end do
	        end do

      end if

      end if

        
 end do loop    

 1000 return


  end subroutine
  
  
  end module