!
! A serial map program
!
! Map the square function sqeezing all of this data through 1 CPU.


subroutine nonCudaMap(in,out,size) 

real, dimension(1:size), intent(in) :: in
real, dimension(1:size), intent(out) :: out
integer, intent(in) :: size

        do i = 1, size
                out(i) = in(i)*in(i)
        end do

end subroutine nonCudaMap

program SerialMap

        implicit none
        integer i, size, iargc
        character(len=32) :: arg
        real, dimension(:), allocatable :: h_in, h_out

        if (iargc() .ne. 1) then
                print *,'Usage: map-non-cuda #'
                stop
        end if

        call getarg(1,arg)
        read(arg,*) size

        allocate(h_in(1:size));
        allocate(h_out(1:size));

        do i = 1, size
                h_in(i) = i
        end do

        call nonCudaMap(h_in, h_out,size)

        do  i = 1, size
                print *,i, h_in(i), h_out(i)
        end do

        deallocate(h_in);
        deallocate(h_out);

end program SerialMap
