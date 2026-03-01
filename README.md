# Multi-Group Cross Section Generation

This project is a Julia based tool that uses continuous-energy nuclear data to generate multigroup cross sections. The code assesses basic reactor physics phenomena such as Westcott g-factor estimates, resonance self-shielding using the Narrow Resonance (NR) approximation, Doppler broadening, and continuous weighting spectrum creation. The existing implementation is designed to assess **Plutonium-239 (Pu-239)** cross sections over a 12-group energy structure, using a **Uranium-235 (U-235)** Watt fission spectrum in the fast region as the basis for a continuous weighting spectrum.

<p align="center">
  <img src="https://github.com/user-attachments/assets/d054f922-059b-4634-8a05-9bab87cfb6e9" width="49%" />
  <img src="https://github.com/user-attachments/assets/84b550c0-27cd-4fae-adfe-bfcc41a5507e" width="50%" />
</p>


## Methodology

### 1. Continuous Weighting Spectrum

To accurately collapse continuous-energy cross sections into multigroup constants, an appropriate intra-group neutron flux spectrum $\phi(E)$ is required. The script synthesizes a continuous 3-region neutron flux spectrum:

* **Thermal Region ($E < 0.2$ eV):** Modeled using a Maxwell-Boltzmann distribution.

$$\phi_{th}(E) = \frac{E}{(k_B T_n)^2} \exp\left(-\frac{E}{k_B T_n}\right)$$



where $T_n$ is the neutron temperature and $k_B$ is the Boltzmann constant ($8.617 \times 10^{-5}$ eV/K).
* **Epithermal Region ($0.2 \text{ eV} \le E < 100 \text{ keV}$):** Modeled using a standard $1/E$ (Fermi) slowing-down spectrum.

$$\phi_{epi}(E) = \frac{C}{E}$$



where $C$ is a continuity constant calculated at the thermal-epithermal boundary.
* **Fast Region ($E \ge 100 \text{ keV}$):** Modeled using an isotope-dependent Watt Fission Spectrum. For this implementation, the **U-235** ENDF/B evaluated parameters are utilized: $a = 0.988 \text{ MeV}^{-1}$ and $b = 2.249 \text{ MeV}^{-1}$.

$$\phi_{fast}(E) = D \cdot \exp(-a E_{MeV}) \sinh(\sqrt{b E_{MeV}})$$



where $D$ is a matching constant computed at the epithermal-fast boundary to ensure spectrum continuity.

### 2. Doppler Broadening

For epithermal energies, resonances must be broadened to account for the thermal motion of the target nuclei. The code implements an optimized SIGMA1 kernel method.
The Doppler width is given by:


$$\Delta = \sqrt{\frac{4 k_B T E}{A}}$$


where $A = 239.0$ for Pu-239.

To optimize performance, the numerical integration is strictly bounded to a local energy window surrounding the target energy ($\pm 8 \Delta \sqrt{E}$). The broadened cross section $\sigma(E, T)$ is evaluated as:


$$\sigma(E, T) = \frac{1}{\Delta \sqrt{\pi}} \int_{E_{min}}^{E_{max}} \sigma(E') \frac{E'}{E} \left[ e^{-\frac{(\sqrt{E} - \sqrt{E'})^2}{\Delta^2}} - e^{-\frac{(\sqrt{E} + \sqrt{E'})^2}{\Delta^2}} \right] dE'$$


Integration is performed numerically using the trapezoidal rule.

### 3. Resonance Self-Shielding

To correct for spatial and energetic flux depression in the vicinity of large resonances, the **Narrow Resonance (NR) Approximation** is applied. A background cross-section $\sigma_0$ (defaulted to 50.0 barns) is assumed.
The self-shielded flux $\phi_{sh}(E)$ is defined as:


$$\phi_{sh}(E) = \phi_{macro}(E) \left( \frac{\sigma_0}{\sigma_{tot}(E, T) + \sigma_0} \right)$$




where $\sigma_{tot}(E, T)$ is the Doppler-broadened total cross section.

### 4. Multigroup Collapsing

The continuous cross sections are collapsed into a standard 12-group structure. Two variants of group constants are calculated:

* **Infinite Dilution:**

$$\sigma_{g, \infty} = \frac{\int_{g} \sigma(E) \phi_{macro}(E) dE}{\int_{g} \phi_{macro}(E) dE}$$


* **Self-Shielded:**

$$\sigma_{g, sh} = \frac{\int_{g} \sigma(E) \phi_{sh}(E) dE}{\int_{g} \phi_{sh}(E) dE}$$



*(Note: The integrals are evaluated using the discrete trapezoidal rule `trapz` provided in the data utilities.)*

### 5. Westcott g-factor

The non-1/v behavior of the absorption and fission cross sections in the thermal region ($10^{-5}$ to $1.0$ eV) is quantified using the Westcott g-factor.


$$g = \frac{\bar{\sigma}_{Maxwell}}{\sigma(E_0)} \frac{2}{\sqrt{\pi}} \sqrt{\frac{T_{mod}}{T_0}}$$


where:
- $$\bar{\sigma}_{\text{Maxwell}}$$ is the cross-section averaged over a pure Maxwellian flux profile at the moderator temperature.
- $$\sigma(E_0)$$ is the reference cross-section evaluated at $$E_0 = 0.0253 \text{ eV}$$ (thermal reference).
- $$T_0 = 293.6 \text{ K}$$ (Reference temperature).

## Code Requirement

```julia
import Pkg
Pkg.add(["CSV", "DataFrames", "Plots", "LaTeXStrings"])
```


