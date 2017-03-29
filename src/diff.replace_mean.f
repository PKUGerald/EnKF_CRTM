14d13
<   USE mpi_module
21,22c20,21
<   integer            :: ii, jj, kk, m, n, fid, iargc, nmi,mstart, mend, ie
<   character (len=10) :: wrf_file, filec, wrf_ifile
---
>   integer            :: ii, jj, kk, m, n, fid, iargc
>   character (len=10) :: wrf_file, filec
24,25d22
<   integer, parameter :: n_cdf = 1000
<   real, dimension(n_cdf-1,2) :: erfinv
36,38c33
<   real, allocatable, dimension(:,:,:)   :: dat3d, work,xq_n,xq_p
<   real, allocatable, dimension(:,:)     :: cen_loc
<   real, allocatable, dimension(:,:,:,:) :: xq_nsend,xq_psend
---
>   real, allocatable, dimension(:,:,:)   :: dat3d, work
40,41d34
< 
<   call parallel_start()
44c37
<   varnum = 23
---
>   varnum = 17
58,59c51,52
<   varname = (/'U_2       ', 'V_2       ', 'W_2       ', 'PH_2      ', 'P         ', &
<               'T_2       ', 'MU_2      ', 'MUB       ', 'PHB       ', 'PB        ', &
---
>   varname = (/'U_2       ', 'V_2       ', 'W_2       ', 'PH_2      ', &
>               'T_2       ', 'MU_2      ', 'MUB       ',  &
61,62c54,55
<               'U10       ', 'V10       ', 'QVAPOR    ', 'QCLOUD    ', 'QRAIN     ', &
<               'QICE      ', 'QSNOW     ', 'QGRAUP    ', 'TSK       '/)
---
>               'U10       ', 'V10       ', 'QVAPOR    ', 'QCLOUD    ', 'QRAIN', &
>               'TSK       '/)
71,81c64,65
<   if(mod(numbers_en,nprocs).eq.0) then
<     nmi=numbers_en/nprocs
<   else
<     nmi=numbers_en/nprocs+1
<   endif
<   mstart=my_proc_id*nmi+1
<   mend=(my_proc_id+1)*nmi
< 
<  do_wrf_var  : do m = 1, varnum
<  var = varname(m)
<  if(my_proc_id==0) write(*,*)var
---
>   do_wrf_var  : do m = 1, varnum
> 
83,92c67,74
<    call wrf_var_dimension ( var, ix, jx, kx, ii, jj, kk )
<    allocate( xa    ( ii, jj, kk, nmi ) )
<    allocate( xm    ( ii, jj, kk ) )
<    allocate( gfs   ( ii, jj, kk ) )
<    allocate( dat3d ( ii, jj, kk ) )
<    allocate( work  ( ii, jj, kk ) )
<    allocate( xq_p  ( ii, jj, kk ) )
<    allocate( xq_n  ( ii, jj, kk ) )
<    allocate( xq_psend    ( ii, jj, kk, nmi ) )
<    allocate( xq_nsend    ( ii, jj, kk, nmi ) )
---
>       var = varname(m)
>       write(*,*)var
>       call wrf_var_dimension ( var, ix, jx, kx, ii, jj, kk )
>       allocate( xa    ( ii, jj, kk, numbers_en ) )
>       allocate( xm    ( ii, jj, kk ) )
>       allocate( gfs   ( ii, jj, kk ) )
>       allocate( dat3d ( ii, jj, kk ) )
>       allocate( work  ( ii, jj, kk ) )
95,109d76
<    if ( kk > 1 ) then
<      call get_variable3d( 'fort.70010', var, ii, jj, kk, 1, gfs )
<    else if ( kk == 1 ) then
<      call get_variable2d( 'fort.70010', var, ii, jj, 1, gfs )
<    endif
< 
<   work = 0.
<   xa = 0.
<   do_ensemble_member : do n=1,nmi   !mstart,mend
<   ie = (n-1)*nprocs+my_proc_id+1
<   !ie = nmi*my_proc_id+n 
<   if (ie <= numbers_en) then
< !.... get ensemble and calculate average
<       write(wrf_file,'(a5,i5.5)')'fort.',i_unit+ie
< !....... get data and sum
111c78
<         call get_variable3d( wrf_file, var, ii, jj, kk, 1, dat3d )
---
>          call get_variable3d( 'fort.70010', var, ii, jj, kk, 1, gfs )
113c80
<         call get_variable2d( wrf_file, var, ii, jj, 1, dat3d )
---
>          call get_variable2d( 'fort.70010', var, ii, jj, 1, gfs )
115,122c82,96
<       xa(:,:,:,n)=dat3d(:,:,:)
<   endif
<   end do do_ensemble_member
< 
<   CALL MPI_Allreduce(sum(xa,4),work,ii*jj*kk,MPI_REAL,MPI_SUM,comm,ierr)
<   xm = work/float(numbers_en)
<   if ( my_proc_id == 0 ) write(*,*)'nmi=',nmi,'xq-mean ',minval(xm),'~',maxval(xm)
<   work = 0.
---
> 
> !.... get ensemble and calculate average
>       work = 0.
>       do_ensemble_member : do n = 1, numbers_en
>          write(wrf_file,'(a5,i5.5)')'fort.',i_unit+n
> !....... get data and sum
>          if ( kk > 1 ) then
>             call get_variable3d( wrf_file, var, ii, jj, kk, 1, dat3d )
>          else if ( kk == 1 ) then
>             call get_variable2d( wrf_file, var, ii, jj, 1, dat3d )
>          endif
>          xa(:,:,:,n)=dat3d(:,:,:)
>          work = work + dat3d
>       end do do_ensemble_member
>       xm = work/float(numbers_en)
126,188c100,128
<   do n=1,nmi !mstart,mend
<   ie = (n-1)*nprocs+my_proc_id+1
<   if (ie <= numbers_en) then
<     work=xa(:,:,:,n)-xm+gfs
<   endif
<   end do
< 
< !!! Removing negative Q-value by Minamide 2015.5.26 
<   if (var=='QCLOUD    ' .or. var=='QRAIN     ' .or. var=='QICE      ' .or. &
<       var=='QGRAUP    ' .or. var=='QSNOW     ') then 
<    xq_p = 0.
<    xq_n = 0.
<    xq_psend = 0.
<    xq_nsend = 0.
<    do n=1,nmi !mstart,mend
<     ie = (n-1)*nprocs+my_proc_id+1
<     if (ie <= numbers_en) then
<      where(work >= 0.) xq_psend(:,:,:,n) = work
<      where(work < 0.) xq_nsend(:,:,:,n) = work
<     endif
<    enddo
<    call MPI_Allreduce(sum(xq_psend,4),xq_p,ii*jj*kk,MPI_REAL,MPI_SUM,comm,ierr)
<    call MPI_Allreduce(sum(xq_nsend,4),xq_n,ii*jj*kk,MPI_REAL,MPI_SUM,comm,ierr)
<    if( my_proc_id == 0 ) write(*,*)'original xq value',minval(work),'~',maxval(work)
<    do n=1,nmi
<     ie = (n-1)*nprocs+my_proc_id+1
<     if(ie<=numbers_en) then
<       where(work < 0.) work = 0.
<       where(xq_p >= abs(xq_n).and.xq_p > 0.) work = work*(xq_p+xq_n)/xq_p
<       where(xq_p < abs(xq_n).or. xq_p == 0.)  work = 0.
<     endif
<    enddo
<    if( my_proc_id == 0 ) write(*,*)'non-negative xq value',minval(work),'~',maxval(work)
<   endif
< 
<   do n=1,nmi !mstart,mend
<   ie = (n-1)*nprocs+my_proc_id+1
<   if (ie <= numbers_en) then
<     write(wrf_file,'(a5,i5.5)')'fort.',o_unit+ie
<     if(my_proc_id==0) write(*,*)'output to ', wrf_file
<     call open_file( wrf_file, nf_write, fid )
< 
<     if ( kk > 1 ) then
<        call write_variable3d(fid, var, ii, jj, kk, 1, work )
<     else if ( kk == 1 ) then
<        call write_variable2d(fid, var, ii, jj,  1, work )
<     endif
<     call close_file( fid )
<   endif
<   end do 
<  
<   deallocate( xa    )
<   deallocate( xm    )
<   deallocate( gfs   )
<   deallocate( dat3d )
<   deallocate( work  )
<   deallocate( xq_p  )
<   deallocate( xq_n  )
<   deallocate( xq_psend  )
<   deallocate( xq_nsend  )
<   call mpi_barrier(comm,ierr)
< 
<  end do do_wrf_var
---
>       do n = 1, numbers_en
>          write(wrf_file,'(a5,i5.5)')'fort.',o_unit+n
>          write(*,*)'output to ', wrf_file
>          call open_file( wrf_file, nf_write, fid )
>           work(:,:,:)=xa(:,:,:,n)-xm(:,:,:)+gfs(:,:,:)
> 
>          if ( kk > 1 ) then
>             call write_variable3d(fid, var, ii, jj, kk, 1, work )
>          else if ( kk == 1 ) then
>             call write_variable2d(fid, var, ii, jj,  1, work )
>          endif
>          call close_file( fid )
>       enddo 
> 
> !      write(wrf_file,'(a5,i5.5)')'fort.',o_unit+numbers_en+1
> !      write(*,*)'output to ', wrf_file
> !      call open_file( wrf_file, nf_write, fid )
> !      if ( kk > 1 ) then
> !         call write_variable3d(fid, var, ii, jj, kk, 1, gfs)
> !      else if ( kk == 1 ) then
> !         call write_variable2d(fid, var, ii, jj,  1, gfs)
> !      endif
> !      call close_file( fid )
> 
>       deallocate( xa    )
>       deallocate( xm    )
>       deallocate( gfs   )
>       deallocate( dat3d )
>       deallocate( work  )
190c130
<   call parallel_finish()
---
>   end do do_wrf_var
192c132
<   if(my_proc_id==0) write(*,*)'!!! Successful completion of replace_mean.exe !!!'
---
>   write(*,*)'!!! Successful completion of replace_mean.exe !!!'
247a188
> !==============================================================================
