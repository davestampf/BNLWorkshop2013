program deviceQuery

        use cudafor
        implicit none

        type(cudaDeviceProp) :: prop
        integer :: nDevices, i, ierr

        ! how many cuda devices on this node?

        ierr = cudaGetDeviceCount(nDevices)

        if (ierr .ne. cudaSuccess) then
                print *, 'Failed to get # of devices'
                stop
        else if (nDevices .eq. 0) then
                print *, 'No Cuda Devices'
                stop
        else
                print *, 'Found ',nDevices,'.'
        end if

        do i = 0, nDevices-1

                ierr = cudaGetDeviceProperties(prop,i)
                if (ierr .ne. cudaSuccess) then
                        print *, 'Failed to get device properties'
                        stop
                endif
                print *, "Device # ",i
		print *, " name = ",trim(prop%name)
		print *, " version = ",prop%major,".",prop%minor
		print *, " total global memory = ",prop%totalGlobalMem
		print *, " shared Memory/Block = ",prop%sharedMemPerBlock
		print *, " registers/block = ",prop%regsPerBlock
		print *, " warp size = ",prop%warpSize
		print *, " Max threads/block = ",prop%maxThreadsPerBlock
		print *, " Max Threads Dim = ",prop%maxThreadsDim(1) ," x ", & 
                         prop%maxThreadsDim(2)," x ",prop%maxThreadsDim(3)
		print *, " Max Grid Size = ",prop%maxGridSize(1)," x ", &
                        prop%maxGridSize(2)," x ",prop%maxGridSize(3)
		print *, " Multi-processor count = ",prop%multiProcessorCount
		print *, " Max Threads/multiprocessor = ",prop%maxThreadsPerMultiProcessor

        end do

end program deviceQuery

