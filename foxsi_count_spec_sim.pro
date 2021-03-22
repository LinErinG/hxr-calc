Pro foxsi_count_spec_sim, eph, photon_flux, int_time=int_time, counting_stat=counting_stat,foxsi_spacecraft=foxsi_spacecraft

  ; This script is to generate simulated FOXSI (sounding rocket or spacecraft) count spectrum
  ; from an input photon spectrum.
  ;
  ; Inputs:
  ;   eph: array of photon energies (keV). MUST BE EVENLY SPACED.
  ;   photon_flux: photon fluxes that correspond to eph (photons s^-1 keV^-1 cm^-2)
  ; 
  ; Keywords:
  ;   int_time: integration time, default is 1 second.
  ;   counting_stat: if set, possion noises will be added to the count spectrum.
  ;   foxsi_spacecraft: if set, the effective area of FOXSI spacecraft will be used;
  ;                     default is to use the effective area of FOXSI-4 sounding rocket.
  ;

  Default, int_time, 1. ;second
  Default, counting_stat, 0
  Default, foxsi_spacecraft, 0
  
  ; Read effective area
  If foxsi_spacecraft eq 0 then file = 'foxsi_effarea_data/foxsi4_eff_area_10shell_cdte_al380.csv' $
    Else file = 'foxsi_effarea_data/foxsi_spacecraft_effective_area_per_module.csv'
  data = read_csv(file)
  energy = data.field1   ;keV
  eff_area = data.field2    ;cm2
  eff_area = interpol(eff_area, energy, eph)

  ; FOXSI count spectrum (counts s^-1 keV^-1)
  count_flux =  photon_flux*eff_area

  ; Convolve with a Gaussian to include detector energy resolution
  energy_resolution = 0.8 ;keV
  bin_size = mean(get_edges(eph, /width)) 
  sigma = energy_resolution/bin_size/(2.*sqrt(2*alog(2)))
  count_flux = GAUSS_SMOOTH(count_flux, sigma, kernel=kernel, /edge_truncate)

  ; Add Poisson noises if needed
  If counting_stat eq 1 then begin
    new_flux = count_flux*int_time
    FOR k=0, n_elements(count_flux)-1 DO BEGIN
      IF count_flux[k] GT 0 THEN Begin
        new_flux[k] = randomu(seed, 1, poisson=count_flux[k]*int_time,/double)
        seed = !NULL
      ENDIF ELSE  new_flux[k] = count_flux[k]*int_time
    ENDFOR
    count_flux = new_flux/int_time
  Endif

  If int_time eq 1 then ytitle_str='count flux (counts/s/keV)' $
    else ytitle_str='counts/keV (in '+strtrim(string(int_time),2)+'s)'
  If foxsi_spacecraft eq 1 then begin
    title_str = 'FOXSI spacecraft simulated count spectrum (per module)'
    xrange = [1,100]
  Endif else begin
    title_str = 'FOXSI-4 sounding rocket simulated count spectrum'
    xrange = [4,30]
  Endelse

  set_line_color
  plot,eph,count_flux*int_time,/ylog,/xlog,/xsty,xr=xrange,color=0,background=1,chars=1.5,$
    xtitle='Energy (keV)',ytitle=ytitle_str,title=title_str,psym=10

  ;stop

End
