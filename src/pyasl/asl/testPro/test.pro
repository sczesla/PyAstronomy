pro test_airtovac
  openw, fw, 'airvac2.test', /get_lun, WIDTH=1000
  printf, fw, "#wvl, airtovac, vactoair"
  for i=0,700 do begin
    wvl = 3000.0d + double(i)*10.0d
    vactoair, wvl, vta
    airtovac, wvl, atv
    printf, fw, wvl, atv, vta, FORMAT="(F20.10, F20.10, F20.10)"
  endfor
  close, fw
end


pro test_mphase
  openw, fw, 'mphase.test', /get_lun, WIDTH=1000
  printf, fw, "#jd, illumfrac"
  for i=0,1000 do begin
    jd = randomu(seed) * 10000.0d + 2456024.0d
    MPHASE, jd, f
    printf, fw, jd, f, FORMAT="(F20.10, F20.10)"
  endfor
  close, fw
end


pro test_moonpos
  openw, fw, 'moonpos.test', /get_lun, WIDTH=1000
  printf, fw, "#jd, ra, dec, dis, geolong, geolat"
  for i=0,1000 do begin
    jd = randomu(seed) * 10000.0d + 2456024.0d
    MOONPOS, jd, ra, dec, dis, geolong, geolat
    printf, fw, jd, ra, dec, dis, geolong, geolat, FORMAT="(F20.10, F20.10, F20.10, F20.10, F20.10, F20.10)"
  endfor
  close, fw
end


pro test_hadec2altaz
  openw, fw, 'hadec2altaz.test', /get_lun, WIDTH=1000
  printf, fw, "#ha, dec, lat, alt, azt"
  for i=0,1000 do begin
    ha = randomu(seed) * 360.0d - 180.0d0
    dec = randomu(seed)*180.0d0-90.0d0
    lat = randomu(seed)*180.0d0-90.0d0
    hadec2altaz, ha, dec, lat, alt, az
    printf, fw, ha, dec, lat, alt, az, FORMAT="(F20.10, F20.10, F20.10, F20.10, F20.10)"
  endfor
  close, fw
end



pro test_co_refract
  openw, fw, 'co_refract.test', /get_lun, WIDTH=1000
  printf, fw, "#ja, obsalt, P, T, aout"
  for i=0,1000 do begin
    a = randomu(seed) * 90.0d
    obsalt = randomu(seed) * 5000.0d
    P = 800.0 + randomu(seed)*400.0d0
    T = -10.0d0 + randomu(seed)*30.0d0 + 273.15d0
    aout = co_refract(a, altitude=obsalt, pressure=P, temperature=T, /to_observed)
    printf, fw, a, obsalt, P, T, aout, FORMAT="(F20.10, F20.10, F20.10, F20.10, F20.10)"
  endfor
  close, fw
end


pro test_co_refract_forward
  openw, fw, 'co_refract_forward.test', /get_lun, WIDTH=1000
  printf, fw, "#a, P, T, R"
  for i=0,1000 do begin
    a = randomu(seed) * 90.0d
    P = 800.0 + randomu(seed)*400.0d0
    T = -10.0d0 + randomu(seed)*30.0d0 + 273.15d0
    R = co_refract_forward(a, P=P, T=T)
    printf, fw, a, P, T, R, FORMAT="(F20.10, F20.10, F20.10, F20.10)"
  endfor
  close, fw
end


pro test_co_aberration
  openw, fw, 'co_aberration.test', /get_lun, WIDTH=1000
  printf, fw, "#jd,  ra,  dec, dra, ddec, eps"
  for i=0,1000 do begin
    jd = randomu(seed) * 10000.0d + 2456024.0d
    ra = randomu(seed)*360.0d0
    dec = randomu(seed)*180.0d0-90.0d0
    co_aberration, jd, ra, dec, d_ra, d_dec
    printf, fw, jd, ra, dec, d_ra, d_dec, FORMAT="(F20.10, F20.10, F20.10, F20.10, F20.10)"
  endfor
  close, fw
end


pro test_co_nutate
  openw, fw, 'co_nutate.test', /get_lun, WIDTH=1000
  printf, fw, "#jd,  ra,  dec, dra, ddec, eps, dpsi, deps"
  for i=0,1000 do begin
    jd = randomu(seed) * 10000.0d + 2456024.0d
    ra = randomu(seed)*360.0
    dec = randomu(seed)*180.0-90.0
    co_nutate, jd, ra, dec, d_ra, d_dec, eps=eps, d_psi=dpsi, d_eps=deps
    printf, fw, jd, ra, dec, d_ra, d_dec, eps, dpsi, deps, FORMAT="(F20.10, F20.10, F20.10, F20.10, F20.10, F20.10, F20.10, F20.10)"
  endfor
  close, fw
end 

pro test_nutate
  openw, fw, 'nutate.test', /get_lun, WIDTH=1000
  printf, fw, "#ja,  longitude,  obliquity"
  for i=0,1000 do begin
    jd = randomu(seed) * 10000.0d + 2456024.0d
    nutate, jd, long, obliq
    printf, fw, jd, long, obliq, FORMAT="(F20.10, F20.10, F20.10)"
  endfor
  close, fw
end 

pro test_sunpos
  openw, fw, 'sunpos.test', /get_lun, WIDTH=1000
  for i=0,1000 do begin
    jd = randomu(seed) * 10000.0d + 2456024.0d
    sunpos, jd, ra, dec, elong, obliquity
    printf, fw, jd, ra, dec, elong, obliquity, FORMAT="(F20.10, F20.10, F20.10, F20.10, F20.10)"
  endfor
  close, fw
end

pro test_eq2hor
  openw, fw, 'eq2hor.test', /get_lun, WIDTH=1000
  printf, fw, "#ra,  dec,  jd,  lon,  lat,  altitude,  alt, az"
  for i=0,1000 do begin
    ra = randomu(seed)*360.0
    de = randomu(seed)*180.0-90.0
    jd = randomu(seed) * 10000.0d + 2456024.0d
    lon = randomu(seed)*360.0
    lat = randomu(seed)*180.0-90.0
    altitude = randomu(seed) * 8000.0d
    eq2hor, ra, de, jd, alt, az, lon=lon, lat=lat, ALTITUDE=altitude
    printf, fw, ra, de, jd, lon, lat, altitude, alt, az, FORMAT="(F20.10, F20.10, F20.10, F20.10, F20.10, F20.10,F20.10,F20.10)"
  endfor
  close, fw
end

pro test_bprecess
  openw, fw, 'bprecess.test', /get_lun, WIDTH=1000
  printf, fw, "#ra,  dec,  ra50,  dec50,  parallax,  mu_radec,  rad_vel"
  for i=0,1000 do begin
    ra = randomu(seed)*360.0
    de = randomu(seed)*180.0-90.0
    r5 = 0.0
    d5 = 0.0
    bprecess, ra, de, r5, d5
    printf, fw, ra, de, r5, d5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, FORMAT="(F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,  F20.10,  F20.10,  F20.10,  F20.10,  F20.10, F20.10)"
  endfor
  for i=0,1000 do begin
    ra = randomu(seed)*360.0
    de = randomu(seed)*180.0-90.0
    mrd = [randomn(seed)*20.0, randomn(seed)*20.0]
    mrd1 = mrd
    r5 = 0.0
    d5 = 0.0
    bprecess, ra, de, r5, d5, MU_RADEC=mrd
    printf, fw, ra, de, r5, d5, 0.0, 0.0, mrd1[0], mrd1[1], mrd[0], mrd[1], 0.0, 0.0, FORMAT="(F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,  F20.10,  F20.10,  F20.10,  F20.10,  F20.10, F20.10)"
  endfor
  for i=0,1000 do begin
    ra = randomu(seed)*360.0
    de = randomu(seed)*180.0-90.0
    mrd = [randomn(seed)*20.0, randomn(seed)*20.0]
    mrd1 = mrd
    par = randomu(seed)
    par1 = par
    r5 = 0.0
    d5 = 0.0
    bprecess, ra, de, r5, d5, MU_RADEC=mrd, PARALLAX=par
    printf, fw, ra, de, r5, d5, par1, par, mrd1[0], mrd1[1], mrd[0], mrd[1], 0.0, 0.0, FORMAT="(F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,  F20.10,  F20.10,  F20.10,  F20.10,  F20.10, F20.10)"
  endfor
  for i=0,1000 do begin
    ra = randomu(seed)*360.0
    de = randomu(seed)*180.0-90.0
    mrd = [randomn(seed)*20.0, randomn(seed)*20.0]
    mrd1 = mrd
    par = randomu(seed)
    par1 = par
    rv = randomu(seed)*100.0
    rv1 = rv
    r5 = 0.0
    d5 = 0.0
    bprecess, ra, de, r5, d5, MU_RADEC=mrd, PARALLAX=par
    printf, fw, ra, de, r5, d5, 0.0, 0.0, mrd1[0], mrd1[1], mrd[0], mrd[1], rv1, rv, FORMAT="(F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,  F20.10,  F20.10,  F20.10,  F20.10,  F20.10, F20.10)"
  endfor
  close, fw
end

pro test_premat
  openw, fw, 'premat.test', /get_lun, WIDTH=1000
  for i=0, 1000 do begin
    e1 = randomu(seed)*2000.0
    e2 = e1 + randomu(seed)*1000.0
    r = premat(e1, e2)
    printf, fw, e1, e2, r[0,0], r[0,1], r[0,2], r[1,0], r[1,1], r[1,2], r[2,0], r[2,1], r[2,2], 0, FORMAT="(F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)"
    r = premat(e1, e2, /FK4)
    printf, fw, e1, e2, r[0,0], r[0,1], r[0,2], r[1,0], r[1,1], r[1,2], r[2,0], r[2,1], r[2,2], 1, FORMAT="(F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)"
  endfor
  close, fw
end

pro test_precess
  openw, fw, 'precess.test', /get_lun, WIDTH=1000
  for i=0, 1000 do begin
    e1 = randomu(seed)*2000.0
    e2 = e1 + randomu(seed)*1000.0
    ra = randomu(seed)*360.0
    de = randomu(seed)*180.0-90.0
    ra1 = ra
    de1 = de
    precess, ra, de, e1, e2
    printf, fw, e1,e2, ra1, de1, ra, de, 0, 0, FORMAT="(F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10)"
  endfor
  for i=0, 1000 do begin
    e1 = randomu(seed)*2000.0
    e2 = e1 + randomu(seed)*1000.0
    ra = randomu(seed)*360.0
    de = randomu(seed)*180.0-90.0
    ra1 = ra
    de1 = de
    precess, ra, de, e1, e2, /FK4
    printf, fw, e1,e2, ra1, de1, ra, de, 1, 0, FORMAT="(F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10)"
  endfor
  for i=0, 1000 do begin
    e1 = randomu(seed)*2000.0
    e2 = e1 + randomu(seed)*1000.0
    ra = randomu(seed)*360.0/ 180.0* 3.14159265359
    de = randomu(seed)*180.0-90.0/ 180.0* 3.14159265359
    ra1 = ra
    de1 = de
    precess, ra, de, e1, e2, /FK4, /RADIAN
    printf, fw, e1,e2, ra1, de1, ra, de, 1, 1, FORMAT="(F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10)"
  endfor
  close, fw
end

pro test_precess_xyz
  openw, fw, 'precess_xyz.test', /get_lun, WIDTH=1000
  for i=0, 1000 do begin
    e1 = randomu(seed,/double)*2000.0
    e2 = e1 + randomu(seed,/double)*1000.0
    x = randomu(seed,/double)*2e9-1e9
    y = randomu(seed,/double)*2e9-1e9
    z = randomu(seed,/double)*2e9-1e9
    x1 = x
    y1 = y
    z1 = z
    precess_xyz, x,y,z, e1, e2
    printf, fw, e1,e2, x1, y1, z1, x, y, z, FORMAT="(F50.25,F50.25,F50.25,F50.25,F50.25,F50.25,F50.25,F50.25)"
  endfor
  close, fw
end

pro test_xyz
  openw, fw, 'xyz.test', /get_lun, WIDTH=1000
  for i=0, 1000 do begin
    e1 = randomu(seed,/double)*200.0+1900.0
    d = randomu(seed,/double)*10000.0+50000.0
    x=0
    y=0
    z=0
    xv=0
    yv=0
    zv=0
    xyz, d, x, y, z, xv, yv, zv, equinox=e1
    printf, fw, d, x,y,z,xv,yv,zv, e1, FORMAT="(F30.15,F30.15,F30.15,F30.15,F30.15,F30.15,F30.15,F30.15)"
  endfor
  close, fw
end

pro test_helio_jd
  openw, fw, 'helio_jd.test', /get_lun, WIDTH=1000
  for i=0, 1000 do begin
    d = randomu(seed)*10000.0+50000.0
    ra = randomu(seed)*360.0
    de = randomu(seed)*180.0-90.0
    hjd = helio_jd(d ,ra , de)
    printf, fw, hjd, d, ra, de, FORMAT="(F30.15,F30.15,F30.15,F30.15,F30.15)"
  endfor
  close, fw
end


pro test_daycnv
  openw, fw, 'daycnv.test', /get_lun, WIDTH=1000
  for i=0, 1000 do begin
     xjd = randomu(seed,/double)*3e7
     yr=0
     mn=0
     day=0
     hr=0
     DAYCNV, xjd, yr, mn, day, hr
     printf, fw, xjd, yr, mn, day, hr, FORMAT="(F30.15,F30.15,F30.15,F30.15,F30.15)"
  endfor
  close, fw
end

pro test_aitoff
  openw, fw, 'aitoff.test', /get_lun, WIDTH=1000
  x = 0
  y = 0
  for i=0, 1000 do begin
    l = 360.0d0 * randomu(seed,/double) - 180.0d0
    b = 180.0d0 * randomu(seed,/double) - 90.0d0
    aitoff, l, b, x, y
    printf, fw, l, b, x, y, FORMAT="(F30.15,F30.15,F30.15,F30.15)"
  endfor
  close, fw
end

pro test_baryvel
  openw, fw, 'baryvel.test', /get_lun, WIDTH=1000
  for i=0, 10000 do begin
    jd = 2.4d6 + double(i)*10.37681d0
    baryvel, jd, 0, a, b
    printf, fw, a[0], a[1], a[2], b[0], b[1], b[2], FORMAT="(F30.15,F30.15,F30.15,F30.15,F30.15,F30.15)"
  endfor
  close, fw

  openw, fw, 'baryvel2.test', /get_lun, WIDTH=1000
  for i=0, 100 do begin
    jd = 2.4d6 + double(i)*1000.37681d0
    for k=0, 10 do begin
      deq = jd + double(k-5) * 10000.0d
      baryvel, jd, deq, a, b
      printf, fw, a[0], a[1], a[2], b[0], b[1], b[2], jd, deq, FORMAT="(F30.15,F30.15,F30.15,F30.15,F30.15,F30.15,F30.15,F30.15)"
    endfor
  endfor
  close, fw
end


pro test_posangle
  openw, fw, 'posangle.test', /get_lun, WIDTH=1000
  for i=0, 1000 do begin
    ra1 = 360.0d0 * randomu(seed,/double) /15.0d0
    de1 = 180.0d0 * randomu(seed,/double) - 90.0d0
    ra2 = 360.0d0 * randomu(seed,/double) /15.0d0
    de2 = 180.0d0 * randomu(seed,/double) - 90.0d0
    posang, 1, ra1, de1, ra2, de2, a
    printf, fw, ra1, de1, ra2, de2, a, FORMAT="(F30.15,F30.15,F30.15,F30.15, F30.15)"
  endfor
  close, fw
end


pro test_helcorr
; pro helcorr,obs_long,obs_lat,obs_alt,ra2000,dec2000,jd $
;            ,corr,hjd,DEBUG=debug
  openw, fw, 'helcorr.test', /get_lun, WIDTH=1000
  corr = 0.0d0
  hjd = 0.0d0
  for i=0, 1000 do begin
    obs_long = 720.0d0 * randomu(seed,/double) - 360.0d0
    obs_lat = 180.0d0 * randomu(seed,/double) - 90.0d0
    obs_alt = 8000.0d0 * randomu(seed,/double)
    ; Convert into hours!
    ra2000 = 360.0d0 * randomu(seed,/double) / 15.0d0
    dec2000 = 180.0d0 * randomu(seed,/double) - 90.0d0
    jd = 2444240.0d0 + randomu(seed,/double) * (100.0d0*365.0d0)
    rjd = jd - 2.4d6
    helcorr, obs_long, obs_lat, obs_alt, ra2000, dec2000, rjd, corr, hjd
    printf, fw, obs_long, obs_lat, obs_alt, ra2000, dec2000, jd, corr, hjd, FORMAT="(F30.15,F30.15,F30.15,F30.15,F30.15,F30.15,F30.15,F30.15)"
  endfor
  close, fw
end


pro create_test_data
 ; note that co_refract has to be compiled first...
 TEST_BPRECESS
 TEST_PREMAT
 TEST_PRECESS
 TEST_PRECESS_XYZ
 TEST_XYZ
 TEST_HELIO_JD
 TEST_DAYCNV
 TEST_AITOFF
 TEST_BARYVEL
 TEST_EQ2HOR
 test_sunpos
 test_nutate
 test_co_nutate
 test_co_aberration
 test_co_refract_forward
 test_co_refract
 test_hadec2altaz
 test_moonpos
 test_mphase
 test_posangle
 test_helcorr
 test_airtovac
end