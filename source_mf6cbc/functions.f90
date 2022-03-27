MODULE functions
interface open_fname
	     module procedure open_fname_1,open_fname_2;	end interface
CONTAINS
!*********************************************************************
!          function to change string character to upper case
!*********************************************************************
    function  upper_case (string)  result (new_string) ! like C
!     -------------------------------------------------------------
!        Convert a string or character to upper case
!         (valid for ASCII or EBCDIC processors)
!     -------------------------------------------------------------
    implicit none
    character (len = *), intent(in) :: string     ! unknown length
    character (len = len(string))   :: new_string ! same length
        character (len = 26), parameter ::               &
            UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ',  &
            lower = 'abcdefghijklmnopqrstuvwxyz'
        integer :: k    ! loop counter
        integer :: loc  ! position in alphabet
        new_string = string             ! copy everything
        do k = 1, len_trim (string)     ! to change letters
           loc = index ( lower, string (k:k) )              ! locate
        if (loc /= 0 ) new_string (k:k) = UPPER(loc:loc) ! convert
        end do ! over string characters
    end function upper_case
    
    
!    文件可用
! *********************************************
    function generate_unit (i) result (value)
       implicit none
	   integer, intent(in) :: i
	   integer             :: value
	   logical             :: connected
          value = i
		  do
		    inquire(unit=value,opened=connected)
			if (.not.connected) exit
			value = value + 1
		  end do
    end function generate_unit
    
!***********************************************************************
!       SUBROUTINE TO OPEN A FILE
!***********************************************************************

   subroutine open_fname_1 (fname,iunit) !the file is positioned at the end-of-file 
      implicit none
	  character(len=*),optional,  intent(in)  :: fname
	  integer,                    intent(out) :: iunit
      logical                                 :: connected,exists

       if(.not.present(fname)) then
           iunit=6
       else
	       inquire(file=fname,exist=exists,opened=connected)
           if(.not.exists) then 
                write(*,*) '>> file not found: ',trim(fname)
                stop
           end if  	              
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
		       iunit = generate_unit(100)
	           open(iunit,file=fname,status='unknown',access='append')
           end if       
	   end if
   end subroutine

   ! 
   subroutine open_fname_2 (fname,iunit,estatus,acceso,formato) !the file is positioned at the end-of-file 
      implicit none
	  character(len=*),           intent(in)  :: fname
	  integer,                    intent(out) :: iunit
	  character(len=*),           intent(in)  :: estatus,acceso,formato
      logical                                 :: connected,exists

	       inquire(file=fname,exist=exists,opened=connected) ! check the file status zzk 0 1
           if(.not.exists) then 
                write(*,*) '>> file not found: ',trim(fname)
                stop
           end if  	              
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
		       iunit = generate_unit(100) 
	           open(iunit,file=fname,status=estatus,access=acceso,form=formato)
           end if       

   end subroutine
  
   
   
   
   !***********************************************************************
!         LOCATE VALUE WITHIN ARRAY (LOCATE)
!***********************************************************************
!
! Given an array "xx" of length "n", and given a value "x", this routine
! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
! must be monotonic单调, either increasing or decreasing.  j=is-1 or j=ie is
! returned to indicate that x is out of range.
!
! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
!-----------------------------------------------------------------------
      subroutine locate (xx,n,is,ie,x,j)

      implicit real*8 (a-h,o-z)

      real*8 xx(n),x

!---- modification from initial subroutine (DANI)

      if (xx(n) == x ) then
	      j = n-1
		  return
	  end if

      if (x <= 0.d0 ) then
	      j = 1
		  return
	  end if

!--------------------------------- end modification
!
! Initialize lower and upper methods:
!
      if(is.le.0) is = 1
      jl = is-1
      ju = ie
      if(xx(n).le.x) then
            j = ie
            return
      end if
!
! If we are not done then compute a midpoint:
!
 10   if(ju-jl.gt.1) then
            jm = (ju+jl)/2
!
! Replace the lower or upper limit with the midpoint:
!
            if((xx(ie).gt.xx(is)).eqv.			 &
               (x.gt.xx(jm))) then
                  jl = jm
            else
                  ju = jm
            endif
            go to 10
      endif
!
! Return with the array index:
!
      j = jl
      
      end subroutine locate
   
      
      
!*********************************************************************************************
!      trilinear interpolation
!*********************************************************************************************
      function interpolation_trilinear (q,xx,yy,zz) result (value)
	  implicit none
	  real*8, intent(in)  :: xx,yy,zz
	  real*8, intent(in)  :: q(2,2,2)
	  real*8              :: value

      value = (1.d0-xx) * (1.d0-yy) * (1.d0-zz) * q(1,1,1) +			   &
                    xx  * (1.d0-yy) * (1.d0-zz) * q(2,1,1) +			   &
                    xx  *       yy  * (1.d0-zz) * q(2,2,1) +    		   &
              (1.d0-xx) *       yy  * (1.d0-zz) * q(1,2,1) +			   &
              (1.d0-xx) * (1.d0-yy) *       zz  * q(1,1,2) +			   &
                    xx  * (1.d0-yy) *       zz  * q(2,1,2) +			   &
                    xx  *       yy  *       zz  * q(2,2,2) +			   &
              (1.d0-xx) *       yy  *       zz  * q(1,2,2)
 
      end function
      
    subroutine open_fname_normal (fname,iunit) !the file is positioned at the beginning-of-file 
      implicit none
	  character(len=*), optional, intent(in)  :: fname
	  integer,                    intent(out) :: iunit
      logical                                 :: connected,exists

       if(.not.present(fname)) then
           iunit=6
       else
	       inquire(file=fname,exist=exists,opened=connected)
           if(.not.exists) then 
                write(*,*) '>> file not found: ',trim(fname)
                stop
           end if   
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
		       iunit = generate_unit(100)
	           open(iunit,file=fname,status='unknown',access='sequential')
           end if       
	   end if
   end subroutine

END MODULE