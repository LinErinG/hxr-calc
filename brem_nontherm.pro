FUNCTION	brem_nontherm, eel, n_e, eph, n_0, vol, stop=stop

	; SSW has a bremsstrahlung cross-section calculator.  I think this works for 
	; relativistic cases, *but I should check that.*

	; Here's the general equation to work with.
	; See https://hesperia.gsfc.nasa.gov/ssw/packages/xray/doc/brm_thin_doc.pdf or 
	; 13.2.11 from Aschwanden's text or other sources.
	; I = 1./(4*!pi*r^2) * n_0 * integral[ cross * v * n_e ]	
	; integrated from photon energy to infinity
	; If using the OSPEX calculators, multiply above equation by extra term of 1/(mc^2)^2 
	; to get the dimensions right.  Then also multiply by n_e_tot

	flux = fltarr(n_elements(eph))	; create flux array
	r = 1.496e13	; 1 AU in cm
	gamma = eel/511+1			; relativistic gamma
	c = 3.e10					; speed of light
	v = c*sqrt(1-1./gamma^2)	; relativistic velocity
	z = 1.2						; avg atomic number of targets
	mc2 = 510.98d+00
	
	if keyword_set( STOP ) then stop
	
	; Total number density of nonthermal electrons (n_e integrated over energy)
	;n_e_tot = tsum( eel, n_e )

	for i=1, n_elements(eph)-1 do begin
		photon_energy = eph[i]
	
		; lower limit on integration is the photon energy.
		; Upper limit is top of EEL array, so be careful not to use this if the electron 
		; energies go past 1 MeV.
	
		brm_bremcross, eel, photon_energy+fltarr(n_elements(eel)), z, cross
		cross = cross/mc2         ; The cross section calculated by brm_bremcross is normalized and 
		                          ; in units of cm^2. To get a cross section in cm^2/keV, we need to
		                          ; to divide it by mc2(eletron rest energy).
		
		ind = where(eel ge photon_energy and finite(cross))
		integral = tsum( eel[ind], cross[ind]*v[ind]*n_e[ind] )    ; use trapezoidal rule for integration
		;; Removed the factor of v - this is assuming the input is electron flux density.
		;integral = int_tabulated( eel[ind], cross[ind]*n_e[ind] )
		flux[i] = 1./(4*!pi*r^2) * n_0 * integral * vol;;; * n_e_tot
	endfor
	
	return, flux; / 300.
	
	; plot, eph, flux, /xlo, /ylo, xra=[1.,100.]

END
