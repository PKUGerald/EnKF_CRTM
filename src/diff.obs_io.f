36d35
< !*******************************************************************************
38d36
< !*******************************************************************************
48d45
< !*******************************************************************************
50d46
< !*******************************************************************************
88,95d83
< 
< !---edited by Minamide 2014.11.18
< ! Radiance data (satellite observation)
< if ( use_radiance ) then
<   if ( .not. use_ideal_obs  ) call get_radiance ( ix, jx, kx, proj, times )
< endif
< 
< 
100,101d87
< 
< !*******************************************************************************
103d88
< !*******************************************************************************
141,144d125
< !---edited by Minamide 2014.11.18
< !. calculate satellite radiance records
< if ( use_radiance) iobs = iobs + raw%radiance%num
< 
155,157c136
< allocate( obs%roi      ( iobs, 3 ) )
< allocate( obs%sat      ( iobs    ) )
< allocate( obs%ch       ( iobs    ) )
---
> allocate( obs%roi      ( iobs, 2 ) )
165,166d143
< obs%sat        = '            '
< obs%ch         = -888888
170d146
< !*******************************************************************************
173d148
< !*******************************************************************************
185d159
< 
306,312d279
< !---edited by Minamide 2014.11.18
< !....... Get Satellite radiance
<    if ( use_radiance ) then
<       call sort_radiance_data( wrf_file, ix, jx, kx, proj, 'radiance', &
<               datathin_radiance, hroi_radiance, vroi_radiance, grid_id )
<    endif
< 
968,1070d934
< 
< !---edited by Minamide 2014.11.18 (for satellite data)
< !=======================================================================================
<   subroutine get_radiance ( ix, jx, kx, proj, times )
<   use constants
<   use namelist_define
<   use mapinfo_define
<   use obs_define
<   use map_utils
<   use netcdf
<   use wrf_tools
<   use radar
<   implicit none
<   type(proj_info), intent(in)         :: proj
<   integer, intent(in)                 :: ix, jx, kx
<   real, dimension(ix, jx      )       :: xlong
<   real, dimension(ix, jx, kx+1)       :: ph
<   real, dimension(ix, jx, kx+1)       :: phb
<   character (len=80)                  :: radiance_file
<   character (len=80)                  :: times
<   character (len=12)                  :: so_time
<   character (len=12)                  :: sat_id
<   integer                             :: i, n, iost, num
<   integer                             :: ch_info,hroi_rad,hroi_drad
<   real                                :: lat, lon, tb, err
<   real                                :: s, h, rx, ry, ir1, jr1, is, js
< 
<   if ( my_proc_id == 0 ) then
<     write(6, *)'   '
<     write(6, *)'---------------------------------------------------'
<     write(6, *)'.... Getting Satellite Radiance Data  ....'
<   endif
< ! get data
<   if ( use_radiance ) then
< 
<   radiance_file = ""
<   radiance_file = 'radiance_'//times(1:4)//times(6:7)//         &
<                times(9:10)//times(12:13)//times(15:16)//'_so'
<   i = len_trim( radiance_file )
<   if ( my_proc_id == 0 ) write(6, *)'.... ', radiance_file(1:i),'  ....'
<   open (10, file = radiance_file(1:i), status = 'old', form = 'formatted', iostat =iost )
<   if( iost .ne. 0 ) then
<       if ( my_proc_id == 0 )write(*,*)'============================================'
<       if ( my_proc_id == 0 )write(*,*)radiance_file(1:i),' does not exist, please check it.'
<       if ( my_proc_id == 0 )write(*,*)'============================================'
< !         stop 'get_radiance'
<   endif
< 
< !...... get the data number
<      num = 0
<      do_get_raw_data_loop_read : do
<         read(10, '(2a12,i12,3f12.3)', iostat = iost ) so_time, sat_id, ch_info, lat, lon, tb
<         if( iost .ne. 0 ) exit
<         num = num + 1
<      end do do_get_raw_data_loop_read
< 
< !...... allocate
<      allocate( raw%radiance%lat( num ) )
<      allocate( raw%radiance%lon( num ) )
<      allocate( raw%radiance%platform( num ) )
<      allocate( raw%radiance%ch( num ) )
<      allocate( raw%radiance%ii( num ) )
<      allocate( raw%radiance%jj( num ) )
<      allocate( raw%radiance%tb( num ) )
<      allocate( raw%radiance%hroi( num ) )
<      allocate( raw%radiance%hroi_d( num ) )
<      allocate( raw%radiance%err( num ) )
< 
< !...... get data
<      rewind(10)
<      num = 0
<      do_get_raw_data_loop : do
< 
< !......... Enkf with odd data, and verify with even data
<         read(10, '(2a12,i12,3f12.3,2i12,f12.3)', iostat = iost ) so_time, sat_id, ch_info, lat, lon, tb, hroi_rad,hroi_drad,err
<         if( iost .ne. 0 ) exit
< !......... calculate radar center's position according to wrf domain grid
<         call latlon_to_ij( proj, lat, lon, is, js )
< !......... evaluate
<            if ( is > 2 .and. is < proj%nx-1 .and. js > 2 .and. js <proj%ny-1) then
<               num = num + 1
<               raw%radiance%lat(num) = lat
<               raw%radiance%lon(num) = lon
<               raw%radiance%platform(num) = sat_id
<               raw%radiance%ch(num) = ch_info
<               raw%radiance%ii(num) = is
<               raw%radiance%jj(num) = js
<               raw%radiance%tb(num) = tb
<               raw%radiance%hroi(num) = (hroi_rad*1000)/proj%dx
<               raw%radiance%hroi_d(num) = (hroi_drad*1000)/proj%dx
<               raw%radiance%err(num) = err
<            else
<            endif
< 
<      end do do_get_raw_data_loop
<      raw%radiance%num = num
<      close (10)
< 
<   endif   !if ( use_radiance ) then
< end subroutine get_radiance
< 
< 
< !=======================================================================================
2113c1977
< !......... T
---
> !......... T:
2403,2493d2266
< 
< 
< !---edited by Minamide 2014.11.18 (for satellite data)
< !=======================================================================================
< subroutine sort_radiance_data( wrf_file, ix, jx, kx, proj, instrument, datathin,hroi, vroi, grid_id )
< 
< use constants
< use namelist_define
< use mpi_module
< use obs_define
< use netcdf
< use map_utils
< !----------------------------------------------------------------------------
< implicit none
< 
< type(proj_info)                      :: proj
< character (len=10), intent(in)       :: wrf_file
< integer, intent(in)                  :: ix, jx, kx
< integer, intent(in)                  :: datathin, hroi, vroi, grid_id
< character (len=10), intent(in)       :: instrument
< 
< integer                              :: i, j, k, nr, n, ista, sta
< 
< real                                 :: x, y, u, v
< real, allocatable, dimension(:,:)    :: data
< character (len=12), allocatable, dimension(:) :: data_sat
< integer, allocatable, dimension(:)   :: data_ch, data_hroi, data_hroi_d
< 
< 
< integer                              :: start_data, inter_data, iroi, ngxn
< 
< !----------------------------------------------------------------------------
< !. get data
< 
< sta = 0
< do nr = 1, raw%radiance%num
<    sta = sta + 1
< enddo
< allocate(data(sta,4))
< allocate(data_sat(sta))
< allocate(data_ch(sta))
< allocate(data_hroi(sta))
< allocate(data_hroi_d(sta))
< 
< ista = 0
< do nr = 1, raw%radiance%num
<     ista = ista + 1
<     data(ista,1) = raw%radiance%tb(nr)
<     data(ista,2) = raw%radiance%err(nr)
<     data(ista,3) = raw%radiance%ii(nr)
<     data(ista,4) = raw%radiance%jj(nr)
<     data_sat(ista) = raw%radiance%platform(nr)
<     data_ch(ista)  = raw%radiance%ch(nr)
<     data_hroi(ista)  = raw%radiance%hroi(nr)
<     data_hroi_d(ista)  = raw%radiance%hroi_d(nr)
< enddo
< !----------------------------------------------------------------------------
< !. data thinning
< start_data = 1
< inter_data = 1
< if ( datathin .lt. -1 ) start_data = abs(datathin) - 1
< if ( abs(datathin) .gt. 1 ) inter_data = abs(datathin)
< 
< iroi = 0
< !inter_data does not work  2014.11.20 Minamide
< do_reports : do n = start_data, ista,inter_data
<   iroi = iroi + 1
<   call cal_hroi ( instrument, grid_id, iroi, ngxn )
<   x = data(n,3)
<   y = data(n,4)
< 
<   if(x.le.1. .and. x.ge.real(ix) .and. y.le.1. .and. y.ge.real(jx) ) cycle do_reports
<   obs%num                 = obs%num + 1
<   obs%dat( obs%num )      = data(n,1)
<   obs%type(obs%num)       = 'Radiance  '
<   obs%err(obs%num)        = data(n,2)
<   obs%position(obs%num,1) = data(n,3)
<   obs%position(obs%num,2) = data(n,4)
<   obs%sat(obs%num)        = data_sat(n)
<   obs%ch(obs%num)         = data_ch(n)
<   obs%roi(obs%num,1)      = data_hroi(n)*ngxn
<   obs%roi(obs%num,3)      = data_hroi_d(n)*ngxn
< !  obs%roi(obs%num,1)      = 1
<   obs%roi(obs%num,2)      = vroi
< end do do_reports
< 
<   return
< 
< end subroutine sort_radiance_data
< !========================================================================================
< 
2523,2524d2295
<    allocate( obs%sat      ( obs%num    ) )
<    allocate( obs%ch       ( obs%num    ) )
2532,2533d2302
<    obs%sat        = '            '
<    obs%ch         = -888888
