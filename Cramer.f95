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
    integer :: n, row, col, i
    real, allocatable :: Matrix1(:,:), b(:), x(:)
    real :: detA, detM, determinant
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
 
       
    ! Find the determinant of M first. print it to screen.
    ! If it is zero, set the Success logical variable and quit.
 
       
    ! Allocate memory for a working matrix for column substituion. Then, for each
    ! column, i, substitute column i with vector b and get that determinant. 
    ! Compute the ith solution.
 
       
    ! deallocate memory for the working matrix.
    
       
 end subroutine Cramer
 
 subroutine ColumnInsert(M, b, n, col, MatOut)
 
    ! This subroutine takes vector b and inserts in into matrix M at column col.
    ! Don't forget to set MatOut = M before you substitute the column in.
       
 
 end subroutine ColumnInsert
 
 function Determinant(M, n) result(Det)
 
    ! This function computes the determinant of matices of size 2 or 3. This
    ! should just be an if...else....endif. One is the formula for a 2x2
    ! and the other is the formula for the 3x3 determinant.
       
 
 end function Determinant