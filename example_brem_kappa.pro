Pro example_brem_kappa, photon_spec

  ; Example script to use brem_nontherm.pro for bremsstrahlung calculation:
  ; calculate x-ray spectrum from a kappa electron distribution and compare with
  ; the results calculated by F_THIN_KAPPA.
  ;
  ; Outputs:
  ;   photon_spec: structure containing photon energy array (in keV) and the corresponding
  ;                photon fluxes (photons s^-1 keV^-1 cm^-2)
  ;
  ; Calling Sequence:
  ;   example_brem_kappa, photon_spec
  ;

  eph = findgen(1000)         ; Photon energies (keV)
  emin = 0.001               ; Min electron energy (keV)
  emax = 2.e3                ; Max electron energy
  bin= 0.01                  ; Electron energy bin (in log space)
  nbin = fix((alog(emax)-alog(emin))/bin)
  eel = findgen(nbin+1)*bin + alog(emin)
  eel = get_edges( exp(eel), /mean )    ;Electron energies (keV)
  kappa = 3.5           ; Kappa index 
  T = 1.                ; Temperature (keV)
  norm_e = 1.e5         ; Normalization factor for electron density distribution
                        ; (equivalent to total electron density, cm^-3)
  n_0 = 1.e10                ; Ambient density (cm^-3)
  vol = (1.d8)^3             ; volume of source (cm^3)

  ; Compute electron density distribution (cm^-3/keV)
  keV2K = 1.16d7        ;conversion factor from keV to K
  e_dist = EKAPPA(eel, [T*keV2K,kappa])   
  n_e = norm_e*e_dist

  ; Photon spectrum (photons s^-1 keV^-1 cm^-2)
  flux = brem_nontherm( eel, n_e, eph, n_0, vol )

  ; Here's the RHESSI kappa thin-target calculation for comparison:
  a0 = n_0*norm_e*vol*1d-49            ; Normalization factor
  flux_compare = F_THIN_KAPPA( eph, [a0,T,kappa,emax] )
 
  loadct, 0
  hsi_linecolors
  plot, eph, flux, /xlo, /ylo, xra=[1.,1000.], thick=2, color=0, background=1,charsize=1.5, $
    xtitle='Photon energy [keV]', ytitle='X-ray flux [phot s!U-1!N keV!U-1!N cm!U-2!N', $
    title = 'Bremsstrahlung X-ray spectrum (kappa)'
  oplot, eph, flux_compare, color=6, thick=2
  al_legend, ['brem_nontherm','f_thin_kappa'], lines=0, thick=2, col=[0,6], box=0, $
    charsize=1.5, linsi=0.3, /right

  photon_spec = create_struct("energy_keV", eph, "photon_flux", flux)

  stop

End
