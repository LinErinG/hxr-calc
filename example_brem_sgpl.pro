Pro example_brem_sgpl
    
	; Example script to use brem_nontherm.pro for bremsstrahlung calculation: 
	; calculate photon spectrum from a single-power-law electron distribution   
  
    eph = findgen(101)        ; Photon energies
    emin = 0.1                ; Min electron energy
    emax = 1.e4          	  ; Max electron energy
    bin= 0.001        		  ; Electron energy bin (in log space)
    nbin = fix((alog(emax)-alog(emin))/bin)
    eel = findgen(nbin+1)*bin + alog(emin)
    eel = get_edges( exp(eel), /mean )		;Electron energies
    n_0 = 1.e10               ; Ambient density
    delta = 3.                ; Electron spectral index
    n_e = 1.e5*eel^(-delta)         ; Electron density distribution that follows a single power law 
    gamma = eel/510.98d00+1     	; relativistic gamma (kinetic E / rest E + 1)
    c = 3e10         		  ; speed of light in cm/s
    v = c*sqrt(1-1./gamma^2)		; relativistic velocity
    vol = (1.d8)^3			  ; volume of source

	; Photon spectrum
	flux = brem_nontherm( eel, n_e, eph, n_0, vol )

    ; Here's the RHESSI power-law thin-target calculation for comparison:
    el_f_tot = tsum( eel, n_e*v )		; Total flux density of nonthermal electrons 
										; (n_e*v integrated over energy)
 	a0 = n_0*el_f_tot*vol
	flux_compare = a0*brm2_thintarget( eph, [1.,delta-0.5,1.e5,delta-0.5,0.1,1e4] )
        	; Brm2_ThinTarget uses electron FLUX density distribution, so the spectral indices 
			; to be put in here is different from the spectral indices for the electron density
			; distribution. The relation is delta_f=delta-0.5 for sub-relativistic cases.

	loadct, 0
	hsi_linecolors
	popen, 'brem_example_sgpl', xsi=7, ysi=5
	plot, eph, flux, /xlo, /ylo, xra=[1.,100.], thick=4, $
		xtitle='Photon energy [keV]', ytitle='X-ray flux [phot s!U-1!N keV!U-1!N cm!U-2!N', $
		title = 'Bremsstrahlung X-ray spectrum (single power law)'
	oplot, eph, flux_compare, color=6, thick=4
	al_legend, ['brem_nontherm','brm2_thintarget'], lines=0, thick=4, col=[0,6], box=0, /right
	pclose
	spawn, 'open brem_example_sgpl.ps'
	
	;stop
	
End
	
