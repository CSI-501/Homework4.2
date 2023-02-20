program CramersRule

    ! System of equations. 2x2, 3x3
    ! The main program is written for you. Read through the comments and
    ! see how the main program works.
    ! 2 Special Notes!!!!!
    ! 1: Take note of how the logial variable 'Success' will either write
    !    the solution or 'No Solution' to the output file.
    ! 2: Take note of how inside the do loop, allocating and deallocating
    !    memory for the arrays Matrix1, b, and x are done so the amount of
    !    memory allocated changes for each system. You cannot allocate more
    !    memory for an array until currently allocated memory is deallocated.
    
    implicit none
    
    ! Declare varialble
    integer :: n, row, col
    real, allocatable :: Matrix1(:,:), b(:), x(:)
    logical :: Success
    
    ! Open the input and output files.
    open(42,file='Data2.txt')
    open(43,file='Data2Out.txt')
    
    ! Solve each system in the input files.
    do
       ! Read in size of first system.
       read(42,*) n
       if (n .eq. 0) exit  ! Quit if zero.
       
       ! Allocate memory for system, right hand side, and solution vector.
       allocate(Matrix1(n,n), b(n), x(n))
       
       ! Read in the system. Ask if you do not understand how this works!
       do row = 1, n
          read(42,*) (Matrix1(row, col), col = 1, n), b(row)
       enddo
       
       ! Use cramers rule to get solution.
       call Cramer(Matrix1, b, n, x, Success)
       
       if (Success) then
          ! Write solution to file
          do row = 1, n
             write(43,*) x(row)
          enddo
          write(43,*)
       else ! This happens when there is no unique solution.
          write(43,*) 'No Solution'
          write(43,*)
       endif
       
       ! clean up memory and go back up to top for next system.
       deallocate(Matrix1, b, x)
       
    enddo
    
    ! close files
    close(42)
    close(43)
 
 end program CramersRule
 
 subroutine Cramer(M, b, n, x, Success)

   ! This subroutine does Cramer's Rule
   ! Declare and initialize your variables first.
   implicit none
   
   real :: M(n,n), b(n), x(n), detM, detW, determinant
   real, allocatable :: workM(:,:)
   integer :: n, i
   logical :: Success
       
   ! Set logic flags.
   Success = .true.

   ! Find the determinant of M first. print it to screen.
   detM = determinant(M, n)
   print*, 'the determinant is ', detM

   ! If it is zero, set the Success logical variable and quit.
   if ( detM == 0.0 ) then
      Success = .false.
      RETURN
   end if 

   ! Allocate memory for a working matrix for column substituion. 
   allocate(workM(n,n))
   
   ! Then, for each column, i, substitute column i with vector b and get 
   ! that determinant. Compute the ith solution.

   do i = 1, n
      workM = M
      call ColumnInsert(workM, b, n, i, workM)
      detW = determinant(workM, n)
      x(i) = detW / detM
   end do
       
   ! deallocate memory for the working matrix.
   deallocate(workM)
       
 end subroutine Cramer
 
 subroutine ColumnInsert(M, b, n, col, MatOut)

   ! This subroutine takes vector b and inserts in into matrix M at column col.
   ! Don't forget to set MatOut = M before you substitute the column in.

   real :: M(n, n), MatOut(n, n), b(n)
   integer :: i, n, col

   MatOut = M
   do i = 1, n
      MatOut(i, col) = b(i)
   enddo

       
 
 end subroutine ColumnInsert
 
 function determinant(M, n) result(Det)
 
   ! Clear the memory for the variables
   implicit none
 
   ! Define the variable types
   real :: Det
   integer :: n
   real :: M(n,n)

   ! Perform determinant for 2x2 matrix.
   if ( n == 2 ) then
      Det = (M(1,1) * M(2,2)) - (M(1,2) * M(2,1))
   
   ! Perform determinant for 3x3 matrix.
   else if ( n == 3 ) then
      Det = (M(1,1) * M(2,2) * M(3,3)) + (M(1,2) * M(2,3) * M(3,1)) + (M(1,3) * M(2,1) * M(3,2)) &
      - (M(1,3) * M(2,2) * M(3,1)) - (M(1,2) * M(2,1) * M(3,3)) - (M(1,1) * M(2,3) * M(3,2))
   end if
       
 
 end function determinant