131c131
< real      :: fac,d,alpha,var,cov,y_hxm,corr_coef,d_ogn
---
> real      :: fac,d,alpha,var,cov,y_hxm,corr_coef
163,168d162
< ! for satellite radiance
< integer :: iob_radmin,iob_radmax
< real, dimension(obs%num) :: yasend_tb, ym_radiance
< real, dimension(ni,nj,nk)     :: xq_n,xq_p
< real, dimension(ni,nj,nk,nm)  :: xq_nsend,xq_psend
< 
336,349d329
< 
< !--ensemble loop for satellite radiance
<  !---every grid is calculated in subroutine xb_to_radiance
< 
< if(raw%radiance%num.ne.0) then
<   yasend_tb=0.
<   do ie = 1, numbers_en+1
<     yasend_tb = 0.0
<     write( filename, '(a5,i5.5)') wrf_file(1:5), iunit+ie-1
<     call xb_to_radiance(filename,proj,ix,jx,kx,xlong,xlat,xland,iob_radmin,iob_radmax,yasend_tb)
<     yasend(iob_radmin:iob_radmax,ie) = yasend_tb(iob_radmin:iob_radmax)
<   enddo
< endif
< 
351,357d330
< 
< !calcurate mean of ya (radmean) by Minamide 2015.9.25
< ym_radiance = 0
< do ie = 1, numbers_en
<   ym_radiance = ym_radiance + ya(:,ie)/float(numbers_en)
< enddo
< 
395c368
<    if ( obstype=='RadarRef  ' .or. obstype=='RadarHydro') then
---
>    if ( obstype=='RadarRef  ' .or. obstype=='RadarHydro') then 
453,502d425
< !
< !!---relaxation for Radiance assimilation by Minamide 2015.3.14
<    d    = fac * var + error * error
<    alpha = 1.0/(1.0+sqrt(error*error/d))
< ! --- for Successive Covariance Localization
<    ngx = obs%roi(iob,1)
<    if (obstype=='Radiance  ') then
<      if(varname=='QCLOUD    ' .or. varname=='QRAIN     ' .or. varname=='QICE      ' .or. & !varname=='QVAPOR    ' .or. 
<                varname=='QGRAUP    ' .or. varname=='QSNOW     ') then
<        if(obs%roi(iob,1) == 0) then
<          update_flag = 0
<        else
<          ngx = obs%roi(iob,1)
<        endif
<      else
<        if(obs%roi(iob,3) == 0) then
<          update_flag = 0
<        else
<          ngx = obs%roi(iob,3)
<          !d = max(fac * var + error * error, y_hxm * y_hxm)
<          !alpha = 1.0/(1.0+sqrt((d-fac * var)/d))
<          !if ( my_proc_id == 0 .and. sqrt(d-fac * var) > error .and. varname=='T         ')&
<          !     write(*,*) 'observation-error inflated to ',sqrt(d-fac * var)
<        endif
<      endif
<      d = max(fac * var + error * error, y_hxm * y_hxm)
<      alpha = 1.0/(1.0+sqrt((d-fac * var)/d))
<      if ( my_proc_id == 0 .and. sqrt(d-fac * var) > error .and. varname=='T         ')&
<           write(*,*) 'observation-error inflated to ',sqrt(d-fac * var)
<      !if  ( my_proc_id == 0 .and.  varname=='T         ') write(*,*)'hroi for dynamic field = ',ngx
<    endif
< !!
< ! --- for Observation Error Inflation
< !   if (obstype=='Radiance  ') then
< !      d = max(fac * var + error * error, y_hxm * y_hxm)
< !      alpha = 1.0/(1.0+sqrt((d-fac * var)/d))
< !      if ( my_proc_id == 0 .and. sqrt(d-fac * var) > error .and. varname=='T         ')&
< !           write(*,*) 'observation-error inflated to ',sqrt(d-fac * var)
< !   endif
< !!  excluding pressure fields
< !   if ((obstype=='Radiance  ') .and.  &
< !       (varname=='PH        ' .or. varname=='MU        ' .or. varname=='PSFC      ' .or. varname=='P         ' .or. &
< !        varname=='U         ' .or. varname=='V         ' .or. varname=='U10       ' .or. varname=='V10       ')) then
< !     update_flag = 0
< !   else if ((obstype /= 'Radiance  ') .and.  &
< !       (varname=='QCLOUD    ' .or. varname=='QRAIN     ' .or. varname=='QICE      ' .or. &
< !        varname=='QGRAUP    ' .or. varname=='QSNOW     ' )) then
< !     update_flag = 0
< !   endif
< !!---relaxation end
504d426
< 
771,772c693
< !     if(varname(1:1).eq.'Q') then
<      if(varname.eq.'QVAPOR    ') then
---
>      if(varname(1:1).eq.'Q') then
802d722
<    ngx = obs%roi(iob,1)
814,819d733
< 
<       !if (obs%roi(iob,1) == 0) then 
<       !  call corr(real(obs%position(iiob,1)-obs%position(iob,1)), real(obs%position(iiob,2)-obs%position(iob,2)), &
<       !          real(obs%position(iiob,3)-obs%position(iob,3)), obs%roi(iob,3), obs%roi(iob,2), corr_coef)
<       !else
<       ! original
822d735
<       !endif
835,851d747
< !!! relaxation by quality contoling
< !   if ((obstype=='Radiance  ') .and. (abs(y_hxm)>max(error*3.,sqrt(fac *var))))then
< !! error = oma_omb method
< !    d_ogn = d
< !     d    = fac * var + max(abs(y_hxm-fac*var*y_hxm/d)*abs(y_hxm),error**2)
< !     alpha =1.0/(1.0+sqrt(max(abs(y_hxm-fac*var*y_hxm/d_ogn)*abs(y_hxm),error**2)/d))
< !!error = y_hxm
< !     d    = fac * var + abs(y_hxm) * abs(y_hxm)
< !     alpha = 1.0/(1.0+sqrt(abs(y_hxm)*abs(y_hxm)/d))
< !! 
< 
< !!!---- OEI
<    if (obs%type(iob)=='Radiance  ') then
<       d = max(fac * var + error * error, y_hxm * y_hxm)
<       alpha = 1.0/(1.0+sqrt((d-fac * var)/d))
<    endif
< !!! relaxation end
980,1015d875
< !!! Removing negative Q-value by Minamide 2015.5.26 
< if ( my_proc_id==0 ) write(*,*) 'updating negative values'
< if(raw%radiance%num.ne.0) then
<  do m=1,nv
<    varname=enkfvar(m)
<    xq_p = 0.
<    xq_n = 0.
<    xq_psend = 0.
<    xq_nsend = 0.
<    if (varname=='QCLOUD    ' .or. varname=='QRAIN     ' .or. varname=='QICE      ' .or. &
<        varname=='QGRAUP    ' .or. varname=='QSNOW     ') then
<     do n=1,nm
<      ie=(n-1)*nmcpu+gid+1
<      if(ie==numbers_en+1) write(*,*)'original xq value',minval(x(:,:,:,m,n)),'~',maxval(x(:,:,:,m,n))
<      if(ie<=numbers_en) then
<        where(x(:,:,:,m,n) >= 0.) xq_psend(:,:,:,n) = x(:,:,:,m,n)
<        where(x(:,:,:,m,n) < 0.) xq_nsend(:,:,:,n) = x(:,:,:,m,n)
<      endif
<     enddo
<     call MPI_Allreduce(sum(xq_psend,4),xq_p,ni*nj*nk,MPI_REAL,MPI_SUM,comm,ierr)
<     call MPI_Allreduce(sum(xq_nsend,4),xq_n,ni*nj*nk,MPI_REAL,MPI_SUM,comm,ierr)
<     if ( my_proc_id==0 ) write(*,*) 'xq_p',minval(xq_p),'~',maxval(xq_p)
<     if ( my_proc_id==0 ) write(*,*) 'xq_n',minval(xq_n),'~',maxval(xq_n)
<     do n=1,nm
<      ie=(n-1)*nmcpu+gid+1
<      if(ie<=numbers_en) then
<        where(x(:,:,:,m,n) < 0.) x(:,:,:,m,n) = 0.
<        where(xq_p >= abs(xq_n).and.xq_p > 0.) x(:,:,:,m,n) = x(:,:,:,m,n)*(xq_p+xq_n)/xq_p
<        where(xq_p < abs(xq_n).or. xq_p == 0.) x(:,:,:,m,n) = 0.
<      endif
<      if(ie<=numbers_en+1) where(xm(:,:,:,m) < 0.) x(:,:,:,m,n) = 0.
<      if(ie==numbers_en+1) write(*,*)'non-negative xq value',minval(x(:,:,:,m,n)),'~',maxval(x(:,:,:,m,n))
<     enddo
<    endif
<  enddo
< endif
