! A very simple CUDA implementation of Map

module map_module
contains 

attributes(global) subroutine map(in,out,size)
        implicit none
        real, dimension(:), intent(in) :: in
        real, dimension(:), intent(out) :: out

        integer, value :: size
        integer i

        i = threadIdx%x
        out(i) = in(i)*in(i)

end subroutine map

end module map_module

program Map_Cuda_Small

        use cudafor
        use map_module

        implicit none
        integer size, i, iargc, ierr
        character (len=32) arg
        type(cudaDeviceProp) :: prop

        real, dimension(:), allocatable :: h_in, h_out
        real, dimension(:), allocatable, device :: d_in, d_out

        if (iargc() .ne. 1) then
                print *, 'Usage map-cuda-small #'
                stop
        end if

        call getarg(1,arg)
        read(arg,*) size

!
! make sure size is <= 1024 for the small case
!

        ierr = cudaGetDeviceProperties(prop,0) 
        if (prop%maxThreadsPerBlock .lt.size) then
                print *, 'Size must be <= ',prop%maxThreadsPerBlock
                stop
        endif 

        allocate(h_in(1:size))
        allocate(h_out(1:size))
        allocate(d_in(1:size))
        allocate(d_out(1:size))

        do i = 1, size
                h_in(i) = i
        end do

        ! transfer data to device

        d_in = h_in

        ! Launch the kernel

        call map<<<1,size>>>(d_in, d_out, size)

        ! transfer data back to host

        h_out = d_out

        do i = 1, size 
                print *, i, h_in(i), h_out(i)
        end do

        deallocate(h_in, h_out, d_in, d_out)

end program Map_Cuda_Small
