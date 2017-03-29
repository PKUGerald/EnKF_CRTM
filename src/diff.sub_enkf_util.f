148,154c148
<   
< !-- use_radiance
<    use_radiance   = .false.
<    datathin_radiance = 999
<    hroi_radiance     = 999
<    vroi_radiance     = 999
<  
---
>    
250,254d243
<    read ( unit = namelist_unit, nml = radiance, iostat = iost )
<    if( iost .ne. 0 ) then
<        write(*,*)'radiance, please check it.'
<        stop 'read_namelist radiance'
<    endif
