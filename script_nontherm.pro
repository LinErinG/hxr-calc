;
; Scripts for practicing the non-thermal bremsstrahlung calculation.
;

; SSW has a bremsstrahlung cross-section calculator: brm_bremcross.pro
; Uses Eqn (4) from Haug, 1997 (closely follows Formula 3BN in Kock & Motz 1959)
; Includes Elwert, 1939 factor (see routine header for more info).

; Here's the general equation to work with (13.2.11 from Aschwanden's text and other sources)
; I = 1./(4*!pi*r^2) * n_0 * integral[ cross * v * n_e ]	; integrated from photon energy to infinity
; Integrate this over a cell volume to get the emission from that volume. (???)

eph = findgen(100)			; observed photon energies
flux = fltarr(n_elements(eph))	; create empty flux array
eel = findgen(1.e4)/10.		; Electron energies (note it doesn't go to infinity!)

r = 1.496e13				; 1 AU in cm
n_0 = 1.e10					; local ion density (chosen)
n_e = 1.e5*eel^(-4)			; electron density (chosen; this is our input distribution)
vol = (1.d8)^3				; test volume (chosen)
gamma = eel/511+1			; relativistic gamma (kinetic E / rest E + 1)
c = 3.e10					; speed of light in cm/s
v = c*sqrt(1-1./gamma^2)	; relativistic velocity
z = 1.2				        ; avg atomic number of targets

.run
for i=1, n_elements(eph)-1 do begin
	photon_energy = eph[i]

	; lower limit on integration is the photon energy.
	; Upper limit is top of EEL array, so be careful not to use this if the electron 
	; energies go past 1 MeV.

	brm_bremcross, eel, photon_energy+fltarr(n_elements(eel)), z, cross
	ind = where(eel ge photon_energy and finite(cross))
	integral = tsum( eel[ind], cross[ind]/511*v[ind]*n_e[ind] )  ;use trapezoidal rule for integration
	         ; The cross section calculated by brm_bremcross is normalized and in units of cm^2. 
	         ; To get a cross section in cm^2/keV, we need to divide it by mc2(eletron rest energy).
	flux[i] = 1./(4*!pi*r^2) * n_0 * integral * vol
endfor
end

plot, eph, flux, /xlo, /ylo, xra=[1.,100.]

;; This looks like a nice clean line up to 100 keV.

;; This is now coded into a function brem_nontherm.pro
; Below is practice using that routine.  This prints out a plot comparing the output of 
; my routine to that of F_POW.

eph = findgen(101)				; Photon energies
emin = 0.1						; Min electron energy
emax = 1e4						; Max electron energy
bin= 0.001						; Elec energy bin (in log space)
nbin = fix((alog(emax)-alog(emin))/bin)
eel = findgen(nbin+1)*bin + alog(emin)
eel = get_edges( exp(eel), /mean)
n_0 = 1.e10						; Ambient density
delta = 3.						; Electron spectral index
n_e = 1.e5*eel^(-delta)			; Electron density distribution, single power law
gamma = eel/510.98d00+1     ; relativistic gamma (kinetic E / rest E + 1)
c = 3e10         ; speed of light in cm/s
v = c*sqrt(1-1./gamma^2)

; Total flux density of nonthermal electrons (n_e*v integrated over energy)
el_f_tot = tsum( eel, n_e*v )
vol = (1.d8)^3

flux = brem_nontherm( eel, n_e, eph, n_0, vol );, /stop )

; Here's the RHESSI power-law thin-target calculation for comparison:
;;;brm2_distrn, eel, delta, delta, min(eel), 20., max(eel), fcn	; Normalized to integrate to 1.
;;;fcn *= int_tabulated(eel,n_e)
a0 = n_0*el_f_tot*vol
flux_compare = a0*brm2_thintarget( eph, [1.,delta-0.5,1.e5,delta-0.5,0.1,1e4] )
        ; Brm2_ThinTarget uses electron FLUX density distribution, so the spectral indices 
	; to be put in here is different from the spectral indices for our electron density
	; distribution. The relation is delta_f=delta-0.5 for sub-relativistic cases.

loadct, 0
hsi_linecolors
popen, 'brem_test', xsi=7, ysi=5
plot, eph, flux, /xlo, /ylo, xra=[1.,100.], thick=4, $
	xtitle='Photon energy [keV]', ytitle='X-ray flux [phot s!U-1!N keV!U-1!N cm!U-2!N', $
	title = 'Bremsstrahlung X-ray spectrum from electrons with alpha='+strtrim(delta,2)
oplot, eph, flux_compare, color=6, thick=4
al_legend, ['brem_nontherm','brm2_thintarget'], lines=0, thick=4, col=[0,6], box=0, /right
pclose
spawn, 'open brem_test.ps'

print, alog10(flux[10]/flux[99])
print, median( flux_compare / flux )
print, median( flux / flux_compare )
plot, flux_compare / flux, /xlog, /ylog, xrange=[1.,100], yrange=[0.1,10]
;print, n_e_tot

;;; scratch work, just reminding myself of some conversions.
v = findgen(1000)/1000.*1.e11
c = 3.e10
m = 9.1d-28
erg2keV = 6.24e8
mc2 = 511
gamma = 1./sqrt(1-(v/c)^2)
ke = (gamma-1)*mc2
ke2 = 0.5*m*v^2*erg2kev
plot, v, ke;, yra=[0,100]
oplot, v, ke2

; Practice using the cross section calculator.  This computes the emission at a 
; given photon energy for a variety of initial electron energies.
; Consider electron energies 1-1000 keV.  This means this should not be used if the 
; electron energies go up past 1 MeV.
photon_energy = 10.
eel = findgen(1.e3)/10.
eph = fltarr(1.e3) + photon_energy
z = 1.
brm_bremcross, eel, eph, z, cross
plot, eel, cross, /ylog, charsi=1.5, xtit='Electron energy [keV]', $
	ytit='Bremsstrahlung cross section [cm!U2!N?]', $
	title = 'Cross section for '+strtrim(photon_energy,2)+' keV photons'

