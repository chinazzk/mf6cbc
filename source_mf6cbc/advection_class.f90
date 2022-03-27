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
  subroutine read_flux_from_mf6_ (mf6,advection,geo)
	 use array_class
	 use geometry_class
	 use functions,only: generate_unit,open_fname, upper_case
	 use variable_def ! mf6_cl
	 implicit none
	 type(mf6_cl),      intent(in)     :: mf6
	 type(geometry_cl),  intent(in)    :: geo
	 type(advection_cl), intent(inout) :: advection
     character(len=16)                 :: text,text2
     character(len=16), pointer        :: auxtxt(:)
	 real*8                            :: dx,dy,dz
     integer*4                         :: kstp,kper,ncol,nrow,nlay
     integer                           :: iunit,imod,jmod,kmod,i,j,k
	 integer                           :: inbud
     real*8                            :: deltd,pertimd,totimd,vald(20)
     integer                           :: ndim1,ndim2,ndim3
     integer                           :: n,m, ncon,ipos,ncells,id,idx1,idx2,idy1,idy2,idz1,idz2
     real                              :: q,face1,face2
     integer                           :: itype,ndat,nlist
     real*8,  save, pointer            :: buff2d(:,:), buff(:), id1(:),id2(:),buffq(:,:)
     integer,pointer                   :: ia(:), ja(:)
     integer,pointer                   :: icell(:), buffi(:,:,:)
     logical                           :: assignflow

!....initialize

     kper = advection%kper
     kstp = advection%kstp
     
     ncol = mf6%nx
     nrow = mf6%ny
     nlay = mf6%nz
     ncells = ncol*nrow*nlay
     allocate(icell(ncells))
     allocate(buffi(ncol,nrow,nlay))
     allocate(buffq(ncells,ncells))
     buffq=0.d0 
     
     icell = (/(i,i=1,ncol*nrow*nlay)/)
     buffi = reshape(source=icell,shape=[ncol,nrow,nlay])
    
     inbud = mf6%unitCBC

     assignflow = .FALSE.

!....read one time step

     loop: do
         read(inbud,end=1000,err=1000) kstp,kper,text,ndim1,ndim2,ndim3
         if (ndim3 .lt. 0) then
             ndim3 = -ndim3
         end if
      ! check if right stress period and time step is correct
      ! if still not there, exit and continue with the same velocity field
         if (kper > advection%kper  .or. (kper == advection%kper .and. kstp > advection%kstp) ) then
             backspace(inbud) ! 返回上一行
             exit
         end if
         if (kper == advection%kper  .and. kstp == advection%kstp ) assignflow =.TRUE. ! 
! start reading one internal flow term
         itype=0
         read(inbud) itype,deltd,pertimd,totimd ! double 8
         if (itype.eq.1) then
             j= ndim1*ndim2*ndim3
             if (.not.associated(buff)) allocate(buff(j))
             read(inbud) buff
         end if
         if (itype.eq.6) then
              do n=1,4
                  read(inbud) text2
              end do
              read(inbud) ndat
              if (.not.associated(auxtxt)) allocate(auxtxt(ndat-1))
              read(inbud) (auxtxt(n),n=1,ndat-1)
              read(inbud) nlist
              if (.not.associated(id1)) allocate(id1(nlist))
              if (.not.associated(id2)) allocate(id2(nlist))
              if (.not.associated(buff2d)) allocate(buff2d(ndat,nlist))
              read(inbud) ((id1(n),id2(n),(buff2d(i,n),i=1,ndat)),n=1,nlist)
          end if
!-----------------------------------------------------------------------
!     associate internal flow terms to darcy fluxes
!-----------------------------------------------------------------------      
          if (assignflow) then
              text = upper_case (text)
              if(text.eq.'   CONSTANT HEAD') then
              ! do nothing
              elseif(text.eq.'    FLOW-JA-FACE') then
                  iunit=100
                  open(iunit,file='cell_face_flow.txt') ! 创建文件
                  call open_fname('cell_face_flow.txt', iunit)
                  do n=1, 800 ! 按 cell num 循环
                      print *, 'This is cell: ', n
                      ncon = mf6%ia(n+1) - mf6%ia(n) - 1 ! 该cell 的 connections 数量 减去了自身
                      if (ncon .lt. 0) ncon=0
                      print *, 'number of connected cells is: ', ncon
                      do ipos = mf6%ia(n)+1, mf6%ia(n+1)-1 ! 该cell的connection 序号 包括本身
                          m = mf6%ja(ipos) ! cell number
                          q = buff(ipos) ! inter flow 该cell #flows are positive for the cell in question.  轮到哪个cell 她所在的流速为正
                          buffq(n,m) = q ! 赋值 构成 待选数组
                          print *,  n,m, q 
                          write(iunit, '(i8.3,i8.3, g21.14E2)')  n,m, q ! 'cell n, cell m, flow : '
                      end do
                  end do        
! assign buffer to qx,qy,qz
! Initialize velocity field
                  if (.not.associated(advection%qx%values)) advection%qx = make_array_ (0.d0,mf6%nx+1,mf6%ny,  mf6%nz ) ! cols+1 rows layers
                  if (.not.associated(advection%qy%values)) advection%qy = make_array_ (0.d0,mf6%nx,  mf6%ny+1,mf6%nz ) ! cols+1 rows+1 layers
                  if (.not.associated(advection%qz%values)) advection%qz = make_array_ (0.d0,mf6%nx,  mf6%ny,  mf6%nz+1) ! cols+1 rows layers+1
                  advection%qx%values = 0.d0
                  advection%qy%values = 0.d0
                  advection%qz%values = 0.d0
                  !open(111,file='check_face_flow.txt') ! 创建文件记录输出结果用以验证
                  do kmod=1,nlay      ! 层
                    do imod=1,nrow     ! 行
                      do jmod=1,ncol ! 列
                          id = buffi(jmod,imod,kmod)
!qx.........................................................
                          i=jmod        ! 相同 rw3d的行 等于 modflow的列？？
                          j=nrow-imod+1 !坐标系方向相反 y
                          k=nlay-kmod+1 !z方向相反 
                          !face1=jmod-1! 下标
                          face2=jmod+1!
                          if (face2 .le. ncol) then
                              idx2 = buffi(face2,imod,kmod)
                              advection%qx%values(i+1,j,k) = - buffq(id,idx2)
                          end if
                          !!both face exist in buffq
                          !if (face1 .ge. 1 .and. face2 .le. ncol) then
                          !    idx1 = buffi(face1,imod,kmod) ! face1
                          !    advection%qx%values(i,j,k) = buffq(id,idx1)
                          !    write(111,'(i8.3,i8.3, g21.14E2)' ) id,idx1,buffq(id,idx1)
                          !    idx2 = buffi(face2,imod,kmod) ! face2 
                          !    advection%qx%values(i+1,j,k) = buffq(id,idx2)
                          !    write(111,'(i8.3,i8.3, g21.14E2)' ) id,idx2,buffq(id,idx2)
                          !!only face2 exist in buffq    
                          !elseif (face1 .lt.1 .and. face2 .le. ncol) then
                          !    idx2 = buffi(face2,imod,kmod) ! face2 
                          !    advection%qx%values(i+1,j,k) = buffq(id,idx2)
                          !    write(111,'(i8.3,i8.3, g21.14E2)' ) id,idx2,buffq(id,idx2)
                          !!only face1 exist in buffq
                          !elseif (face1 .ge. 1 .and. face2 .gt. ncol) then
                          !    idx1 = buffi(face1,imod,kmod) ! face1
                          !    advection%qx%values(i,j,k) = buffq(id,idx1)
                          !    write(111,'(i8.3,i8.3, g21.14E2)' ) id,idx1,buffq(id,idx1)
                          !end if
!qy.........................................................    
                          j=nrow-imod
                          face2=imod+1
                          if (face2 .le. nrow) then
                              idy2 = buffi(jmod,face2,kmod)
                              advection%qy%values(i,j+1,k) = buffq(id,idy2)
                          end if
                          !face1=imod-1
                          !face2=imod+1
                          !if (face1 .ge. 1 .and. face2 .le. nrow ) then !限制下标范围不能超出
                          !    idy1 = buffi(jmod,face1,kmod)
                          !    advection%qy%values(i,j,k) = buffq(id,idy1)
                          !    write(111,'(i8.3,i8.3, g21.14E2)' ) id,idy1,buffq(id,idy1)
                          !    idy2 = buffi(jmod,face2,kmod)
                          !    advection%qy%values(i,j+1,k) = buffq(id,idy2)
                          !    write(111,'(i8.3,i8.3, g21.14E2)' ) id,idy2,buffq(id,idy2)                                
                          !elseif (face1 .lt. 1 .and. face2 .le. nrow) then ! 不大于
                          !    idy2 = buffi(jmod,face2,kmod)
                          !    advection%qy%values(i,j+1,k) = buffq(id,idy2)
                          !    write(111,'(i8.3,i8.3, g21.14E2)' ) id,idy2,buffq(id,idy2)                      
                          !elseif (face1 .ge. 1 .and. face2 .gt. nrow .and.j>0) then
                          !    idy1 = buffi(jmod,face1,kmod)
                          !    advection%qy%values(i,j,k) = buffq(id,idy1)
                          !    write(111,'(i8.3,i8.3, g21.14E2)' ) id,idy1,buffq(id,idy1)                              
                          !end if

!qz.........................................................   
                          j=nrow-imod+1
                          k=nlay-kmod
                          face2=kmod+1
                          if (face2 .le. nlay) then
                              idz2 = buffi(jmod,imod,face2)
                              advection%qz%values(i,j,k+1) = buffq(id,idz2)
                          end if          
                          !face1=kmod-1
                          !face2=kmod+1
                          !j=nrow-imod+1
                          !k=nlay-kmod
                          !if (face1 .ge. 1 .and. face2 .le. nlay ) then
                          !    idz1 = buffi(jmod,imod,face1)
                          !    advection%qz%values(i,j,k) = buffq(id,idz1)
                          !    write(111,'(i8.3,i8.3, g21.14E2)' ) id,idz1,buffq(id,idz1)
                          !    idz2 = buffi(jmod,imod,face2)
                          !    advection%qz%values(i,j,k+1) = buffq(id,idz2)
                          !    write(111,'(i8.3,i8.3, g21.14E2)') id,idz2,buffq(id,idz2)
                          !elseif (face1 .lt. 1 .and. face2 .le. nlay) then
                          !    idz2 = buffi(jmod,imod,face2)                   
                          !    advection%qz%values(i,j,k+1) = buffq(id,idz2)
                          !    write(111,'(i8.3,i8.3, g21.14E2)') id,idz2,buffq(id,idz2)
                          !elseif (face1 .ge. 1 .and. face2 .gt. nlay .and. k>0) then
                          !    idz1 = buffi(jmod,imod,face1)
                          !    advection%qz%values(i,j,k) = buffq(id,idz1)
                          !    write(111,'(i8.3,i8.3, g21.14E2)' ) id,idz1,buffq(id,idz1)
                          !end if
                      end do
                    end do
                  end do      
              end if
          end if       
     end do loop    

 1000 return


  end subroutine
  
  
  end module