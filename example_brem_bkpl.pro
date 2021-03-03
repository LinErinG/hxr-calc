Pro example_brem_bkpl
    
    ; Example script to use brem_nontherm.pro for bremsstrahlung calculation: 
    ; calculate x-ray spectrum from a broken-power-law electron distribution   
  
    eph = findgen(1000)        ; Photon energies
    emin = 10.                 ; Min electron energy
    emax = 2.e3                ; Max electron energy
    ebrk = 100.                ; break energy
    bin= 0.01                  ; Electron energy bin (in log space)
    nbin = fix((alog(emax)-alog(emin))/bin)
    eel = findgen(nbin+1)*bin + alog(emin)
    eel = get_edges( exp(eel), /mean )		;Electron energies
    n_0 = 1.e10                ; Ambient density
    delta1 = 3.5               ; Electron spectral index below ebrk
    delta2 = 8.                ; Electron spectral index above ebrk
    
    ; Compute electron density distribution
    brm2_distrn, eel, delta1, delta2, emin, ebrk, emax, e_dist
    n_e = 1.e5*e_dist
    
    gamma = eel/510.98d00+1       ; relativistic gamma (kinetic E / rest E + 1)
    c = 3e10                  ; speed of light in cm/s
    v = c*sqrt(1-1./gamma^2)      ; relativistic velocity
    vol = (1.d8)^3            ; volume of source

    ; Photon spectrum
    flux = brem_nontherm( eel, n_e, eph, n_0, vol )

    ; Here's the RHESSI power-law thin-target calculation for comparison:
    el_f_tot = tsum( eel, n_e*v )    ; Total flux density of nonthermal electrons 
                                     ; (n_e*v integrated over energy)
    a0 = n_0*el_f_tot*vol
    flux_compare = a0*brm2_thintarget( eph, [1.,delta1-0.5,ebrk,delta2-0.5,emin,emax] )
            ; Brm2_ThinTarget uses electron FLUX density distribution, so the spectral indices 
            ; to be put in here is different from the spectral indices for the electron density
            ; distribution. The relation is delta_f=delta-0.5 for sub-relativistic cases.

    loadct, 0
    hsi_linecolors
    popen, 'brem_example_bkpl', xsi=7, ysi=5
    plot, eph, flux, /xlo, /ylo, xra=[1.,1000.], thick=4, $
        xtitle='Photon energy [keV]', ytitle='X-ray flux [phot s!U-1!N keV!U-1!N cm!U-2!N', $
        title = 'Bremsstrahlung X-ray spectrum (broken power law)'
    oplot, eph, flux_compare, color=6, thick=4
    al_legend, ['brem_nontherm','brm2_thintarget'], lines=0, thick=4, col=[0,6], box=0, /right
    pclose
    spawn, 'open brem_example_bkpl.ps'
	
    ;stop
	
End
