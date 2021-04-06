Pro rhessi_count_spec_sim, eph, photon_flux, int_time=int_time, counting_stat=counting_stat,$
        radial_offset_=radil_offset_, atten_state=atten_state, detectors=detectors

  ; This script is to generate simulated RHESSI count spectrum from an input photon spectrum.
  ; 
  ; -----------------------------------------------------------------------------------------------------
  ; DISCLAIMER:
  ; This routine only takes into account photon absorbtion for RHESSI instrument response (i.e. no compton
  ; scattering etc), therefore it works only for energies below ~50 keV. Background is not considered 
  ; either. Users should always use the full RHESSI response if they want to look at a specific event. 
  ; ------------------------------------------------------------------------------------------------------
  ;       
  ; Inputs:
  ;   eph: array of photon energies (keV). MUST BE EVENLY SPACED.
  ;   photon_flux: photon fluxes that correspond to eph (photons s^-1 keV^-1 cm^-2)
  ;
  ; Keywords:
  ;   int_time: integration time, default is 1 second.
  ;   counting_stat: if set, possion noises will be added to the count spectrum.
  ;   radial_offset_: radial offsets [degrees] of center of the image relative to the imaging axis.
  ;                   (See headers of hessi_grm.pro & hessi_drm_4image.pro for more information.)
  ;   atten_state: scalar, 0,1,2,3 ==> none, thin, thick, both attenuators in.
  ;   detectors: array of detector indices (0,...,17 correspond to 1t,2t,5...,9t,1r,2r,...,). 
  ;              default: [1]
  ;   

  Default, int_time, 1. ;second
  Default, counting_stat, 0
  Default, radial_offset_, 0.25
  Default, atten_state, 0
  Default, detectors, [1]

  ; Get RHESSI effective area (note: only diagonal elements here) 
  bin_size = mean(get_edges(eph, /width))
  e_edges = [eph[0]-0.5*bin_size,eph+0.5*bin_size]
  ind = where(e_edges ge 1. and e_edges le 20000.)
  e_edges = e_edges(ind)
  eff_area = rhessi_eff_area2(e_edges,radial_offset_, atten_state, detectors=detectors)
  
  ; FOXSI count spectrum (counts s^-1 keV^-1)
  eph = eph(ind[0:-2])
  photon_flux = photon_flux(ind[0:-2])
  count_flux =  photon_flux*eff_area

  ; Convolve with a Gaussian to include detector energy resolution
  energy_resolution = 1. ;keV
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
  title_str = 'RHESSI simulated count spectrum'
  xrange = [1,500]
  
  set_line_color
  plot,eph,count_flux*int_time,/ylog,/xlog,/xsty,xr=xrange,color=0,background=1,chars=1.5,$
    xtitle='Energy (keV)',ytitle=ytitle_str,title=title_str,psym=10
  al_legend,['attenuator state: '+strtrim(string(round(atten_state)),2), $
             'detectors: '+ strjoin(strtrim(string(detectors),2),' ')], /right,chars=1.2
             
  ;stop
  
End
