# hxr-calc

This directory contains the IDL code to calculate bremsstrahlung x-ray spectrum from electron distributions.

* brem_nontherm.pro: main bremsstrahlung calculation code
* example_brem_* .pro: example scripts to compute X-ray spectrum from various electron distributions
* foxsi_count_spec_sim.pro: routine to simulate FOXSI (rocket or spacecraft) count spectrum from an input photon spectrum
* rhessi_count_spec_sim.pro: routine to simulate RHESSI count spectrum from an input photon spectrum
* foxsi_effarea_data: folder containing FOXSI effective area data, used in foxsi_count_spec_sim.pro
* rhessi_eff_area2.pro: function to compute RHESSI effective area (diagonal elements only), used in rhessi_count_spec_sim.pro

### Examples of Usage
To calculate X-ray spectrum from a single power-law electron distribution:
```
example_brem_sgpl, photon_spec
```
To further simulate the nonthermal count spectrum for FOXSI spacecraft (60s integration time, without noise):
```
foxsi_count_spec_sim, photon_spec.energy_keV,photon_spec.photon_flux, foxsi_spacecraft=1, int_time=60., counting_stat=0
```
To simulate count spectrum for RHESSI (detectors 1&3, thin attenuator)
```
rhessi_count_spec_sim, photon_spec.energy_keV,photon_spec.photon_flux,detectors=[1,3],atten_state=1
```
