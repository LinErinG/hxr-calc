FUNCTION	brem_nontherm, eel, n_e, eph, n_0, vol

	; PURPOSE:
	; Calculate the thin-target bremsstrahlung emission from an arbitrary electron distribution
	; 
	; DESCRIPTION:
	; The general equation for bremsstrahlung calculation is:
	; I = 1./(4*pi*r^2) * n_0 * integral[ cross_section * v * n_e ] 
	; (integrated from photon energy to infinity).
	; See e.g. 13.2.11 from Aschwanden's text (Physics of the Solar Corona) 
	; or Holman et al 2018.
	; The cross section is calculated by SSW brm_bremcross.pro, which uses Eqn (4) from Haug, 1997 
	; (closely follows Formula 3BN in Kock & Motz 1959) and includes Elwert, 1939 factor (see routine 
	; header for more info).
	; 	
	; CALLING SEQUENCE:
	; Result = brem_nontherm( eel, n_e, eph, n_0, vol )
	; 
	; INPUTS:
	; eel: array of electron energies (keV)
	; n_e: electron (differential) density (cm^-3/keV) corresponding to eel 
	; eph: array of photon energies (keV)
	; n_0: local ion density (cm^-3)
	; vol: volume of source (cm^3)
	;
	; OUTPUTS:
	; Photon fluxes (photons s^-1 keV^-1 cm^-2) that correspond to the input eph array.
	;
	; HISTORY:
	; March 2021   Yixian Zhang & Lindsay Glesener 
	;
	
	
	flux = fltarr(n_elements(eph))		
	r = 1.496e13	                	; 1 AU in cm
	mc2 = 510.98d+00                        ; electron rest energy
	gamma = eel/mc2+1	  		; relativistic gamma
	c = 3.e10				; speed of light
	v = c*sqrt(1-1./gamma^2)		; relativistic velocity
	z = 1.2					; avg atomic number of targets
	
	for i=1, n_elements(eph)-1 do begin
		photon_energy = eph[i]
		brm_bremcross, eel, photon_energy+fltarr(n_elements(eel)), z, cross
		cross = cross/mc2         ; The cross section calculated by brm_bremcross is normalized and 
		                          ; in units of cm^2. To get a cross section in cm^2/keV, we need to
		                          ; divide it by mc2(electron rest energy).			  
		ind = where(eel ge photon_energy and finite(cross))
		integral = tsum( eel[ind], cross[ind]*v[ind]*n_e[ind] )    ; use trapezoidal rule for integration
			; Note: Lower limit on integration is the photon energy. Upper limit is top of EEL array. 
		flux[i] = 1./(4*!pi*r^2) * n_0 * integral * vol
	endfor

	return, flux

END
