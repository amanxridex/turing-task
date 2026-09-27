# Task Metadata
- **Project:** Nvidia_STEM_SciCode_3k
- **Domain:** Physics / Condensed Matter Physics
- **Subdomain:** 2D Materials & Semiclassical Carrier Transport
- **Title:** Carrier Transport Relaxation and Pauli-Blocking Modulated Mobility in Bilayer Graphene
- **Paper Reference:** Amit Varshney and SSZ Ashraf, "Transport Relaxation Mechanisms in Bilayer Graphene: Effects of Pauli Blocking", arXiv:2609.19643v1 [cond-mat.mes-hall] (2026). https://arxiv.org/abs/2609.19643

---

# Subproblem 1: Elastic Disorder Scattering Rates in Bilayer Graphene

## Prompt
In Bernal-stacked bilayer graphene (BLG), low-energy electrons behave as massive chiral quasiparticles described by a two-band effective Hamiltonian with an approximately parabolic energy dispersion \( E_k = \frac{\hbar^2 k^2}{2 m^*} \), where \( m^* \approx 0.033 \, m_e \) is the effective electron mass and \( k \) is the magnitude of the 2D wavevector. Due to the two-sublattice pseudospin structure of BLG, the low-energy spinor eigenstates have winding number 2, resulting in a chiral wave-function overlap factor between initial wavevector \(\mathbf{k}\) and scattered wavevector \(\mathbf{k}'\) given by \( u_{kk'} = \frac{1 + \cos(2\theta)}{2} = \cos^2\theta \), where \(\theta\) is the scattering angle between \(\mathbf{k}\) and \(\mathbf{k}'\).

Carrier momentum relaxation by weak, finite-range neutral impurities (NI) is modeled within the first Born approximation using a Gaussian impurity potential \( V(r) = V_0 \exp\left(-\frac{r^2}{2 a^2}\right) \), whose 2D Fourier transform is \( V(q) = 2\pi a^2 V_0 \exp\left(-\frac{q^2 a^2}{2}\right) \), where \( a \) is the Gaussian correlation length and \( V_0 \) is the potential amplitude. For elastic scattering at carrier energy \( E_k \), energy conservation imposes \( k' = k = k_F = \frac{\sqrt{2 m^* E_k}}{\hbar} \), and the momentum transfer is \( q = 2 k_F \sin(\theta/2) \). The Boltzmann transport relaxation rate includes the transport weighting factor \( (1 - \cos\theta) \). Integrating over the scattering angle yields the closed-form expression for the neutral impurity scattering rate:
\[
\frac{1}{\tau_{NI}} = \frac{2\pi n_{NI}}{\hbar} (2\pi a^2 V_0)^2 D(E_F) R_{BLG}(x_{NI})
\]
where \( n_{NI} \) is the 2D neutral impurity density, \( D(E_F) = \frac{g_s g_v m^*}{2\pi \hbar^2} = \frac{2 m^*}{\pi \hbar^2} \) is the total electronic density of states of BLG (with spin degeneracy \( g_s = 2 \) and valley degeneracy \( g_v = 2 \)), \( x_{NI} = 2 k_F^2 a^2 \), and \( R_{BLG}(x) \) is the dimensionless BLG transport kernel:
\[
R_{BLG}(x) = e^{-x} \left[ \frac{1}{2} I_0(x) - \frac{3}{4} I_1(x) + \frac{1}{2} I_2(x) - \frac{1}{4} I_3(x) \right]
\]
Here \( I_n(x) \) is the modified Bessel function of the first kind of order \( n \). In the short-range contact limit where the correlation length \( a \to 0 \) while the integrated potential strength \( U_0 = 2\pi a^2 V_0 \) is held constant, \( x_{NI} \to 0 \), \( R_{BLG}(0) = \frac{1}{2} \), yielding the contact-limit relaxation rate \( \frac{1}{\tau_{NI}^{(0)}} = \frac{2 n_{NI} m^* U_0^2}{\hbar^3} \). For comparison, a conventional two-dimensional electron gas (2DEG) without chiral wavefunction overlap has the transport kernel \( R_{2DEG}(x) = e^{-x}[I_0(x) - I_1(x)] \).

Write a Python function `blg_elastic_disorder_rates` that calculates the Fermi wavevector \( k_F \), the dimensionless argument \( x_{NI} \), the BLG transport kernel \( R_{BLG}(x_{NI}) \), the conventional 2DEG kernel \( R_{2DEG}(x_{NI}) \), the neutral impurity scattering rate \( \tau_{NI}^{-1} \), and optionally the contact-limit rate \( \tau_{NI}^{(0)-1} \).

```python
def blg_elastic_disorder_rates(Ek, m_star=3.0060966e-32, n_NI=1.0e15, a_NI=0.5e-9, V0_NI=1.602176634e-19, U0_contact=None):
    """
    Compute elastic neutral impurity scattering rates and transport kernels in bilayer graphene.

    Parameters
    ----------
    Ek : float
        Carrier kinetic energy in Joules (J). Must be a positive finite float.
    m_star : float, optional
        Effective electron mass in kilograms (kg). Default is 0.033 * m_e = 3.0060966e-32 kg.
    n_NI : float, optional
        Sheet density of neutral impurities in m^-2. Default is 1.0e15 m^-2 (10^11 cm^-2).
    a_NI : float, optional
        Gaussian correlation length of the impurity potential in meters (m). Default is 0.5e-9 m (0.5 nm).
    V0_NI : float, optional
        Amplitude of the neutral impurity potential in Joules (J). Default is 1.0 eV = 1.602176634e-19 J.
    U0_contact : float or None, optional
        Integrated contact potential strength in J*m^2. If provided, evaluates contact-limit rate. Default is None.

    Returns
    -------
    dict
        A dictionary containing:
        - 'rate_NI' : float, Gaussian neutral impurity scattering rate in s^-1.
        - 'rate_NI_contact' : float, contact-limit scattering rate in s^-1 (0.0 if U0_contact is None).
        - 'kF' : float, Fermi wavevector in m^-1.
        - 'x_NI' : float, dimensionless argument 2*(kF*a_NI)^2.
        - 'R_BLG' : float, dimensionless BLG chiral transport kernel value.
        - 'R_2DEG' : float, dimensionless conventional 2DEG transport kernel value.

    Raises
    ------
    ValueError
        If any input parameter is None, negative, zero, or non-finite.
    """
    pass
```

## Background
Momentum relaxation in bilayer graphene devices is fundamentally governed by the interplay between the electronic band structure and disorder potentials. Unlike monolayer graphene where carriers are massless Dirac fermions with linear dispersion, Bernal-stacked bilayer graphene features an approximately parabolic band dispersion near the charge neutrality point. Consequently, the low-energy density of states is nearly constant with respect to carrier energy, mirroring a conventional two-dimensional electron gas.

However, the pseudospin chirality of bilayer graphene introduces an angular overlap factor that strongly suppresses small-angle scattering and modifies the backscattering weighting compared to conventional semiconductor quantum wells. When carriers scatter elastically from neutral disorder, such as point defects, vacancies, or neutral adsorbates modeled by a Gaussian potential, the finite range of the potential cuts off large momentum transfers. Transforming the angular collision integral over the chiral wavefunction overlap and transport factor yields a characteristic combination of modified Bessel functions of orders zero through three. In the extreme short-range limit, this kernel approaches exactly one half, demonstrating that chiral pseudospin suppresses the momentum relaxation rate by a factor of two relative to an ordinary two-dimensional electron gas with identical effective mass and scattering potential.

## Testing Template
```python
def test_case_1():
    import numpy as np
    out = blg_elastic_disorder_rates(0.1 * 1.602176634e-19)
    assert np.isclose(out['rate_NI'], 3.010429408887473e+12, rtol=1e-4)
    assert np.isclose(out['kF'], 294303550.5772216, rtol=1e-4)
    assert np.isclose(out['x_NI'], 0.043307289941179605, rtol=1e-4)
    assert np.isclose(out['R_BLG'], 0.46358936248662097, rtol=1e-4)

def test_case_2():
    import numpy as np
    out = blg_elastic_disorder_rates(0.1 * 1.602176634e-19)
    assert not np.isclose(out['R_BLG'], out['R_2DEG'], rtol=1e-2)
    assert out['R_BLG'] < out['R_2DEG']
    assert np.isclose(out['R_2DEG'], 0.9373253772625255, rtol=1e-4)

def test_case_3():
    import numpy as np
    U0_val = 1.57e-37
    out = blg_elastic_disorder_rates(0.05 * 1.602176634e-19, a_NI=1e-12, U0_contact=U0_val)
    expected_contact = 2.0 * 1.0e15 * 3.0060966e-32 * (U0_val**2) / (1.054571817e-34**3)
    assert np.isclose(out['rate_NI_contact'], expected_contact, rtol=1e-5)
    assert np.isclose(out['R_BLG'], 0.5, atol=1e-3)

def test_case_4():
    import numpy as np
    out_low = blg_elastic_disorder_rates(0.02 * 1.602176634e-19)
    out_high = blg_elastic_disorder_rates(0.30 * 1.602176634e-19)
    assert out_low['rate_NI'] > out_high['rate_NI']
    assert out_low['x_NI'] < out_high['x_NI']

def test_case_5():
    import numpy as np
    out = blg_elastic_disorder_rates(0.1 * 1.602176634e-19, n_NI=2.0e15)
    out_ref = blg_elastic_disorder_rates(0.1 * 1.602176634e-19, n_NI=1.0e15)
    assert np.isclose(out['rate_NI'], 2.0 * out_ref['rate_NI'], rtol=1e-6)

def test_case_6():
    import numpy as np
    out = blg_elastic_disorder_rates(0.1 * 1.602176634e-19)
    for key in ['rate_NI', 'rate_NI_contact', 'kF', 'x_NI', 'R_BLG', 'R_2DEG']:
        assert key in out
        assert isinstance(out[key], float)
    assert out['rate_NI_contact'] == 0.0

def test_case_7():
    import numpy as np
    bad_inputs = [None, -1.0, 0.0, float('nan'), float('inf')]
    for bad in bad_inputs:
        try:
            blg_elastic_disorder_rates(bad)
            assert False, f"Expected ValueError for bad Ek={bad}"
        except (ValueError, TypeError):
            pass
```

## Solution
```python
def blg_elastic_disorder_rates(Ek, m_star=3.0060966e-32, n_NI=1.0e15, a_NI=0.5e-9, V0_NI=1.602176634e-19, U0_contact=None):
    """
    Compute elastic neutral impurity scattering rates and transport kernels in bilayer graphene.

    Parameters
    ----------
    Ek : float
        Carrier kinetic energy in Joules (J). Must be a positive finite float.
    m_star : float, optional
        Effective electron mass in kilograms (kg). Default is 0.033 * m_e = 3.0060966e-32 kg.
    n_NI : float, optional
        Sheet density of neutral impurities in m^-2. Default is 1.0e15 m^-2 (10^11 cm^-2).
    a_NI : float, optional
        Gaussian correlation length of the impurity potential in meters (m). Default is 0.5e-9 m (0.5 nm).
    V0_NI : float, optional
        Amplitude of the neutral impurity potential in Joules (J). Default is 1.0 eV = 1.602176634e-19 J.
    U0_contact : float or None, optional
        Integrated contact potential strength in J*m^2. If provided, evaluates contact-limit rate. Default is None.

    Returns
    -------
    dict
        A dictionary containing:
        - 'rate_NI' : float, Gaussian neutral impurity scattering rate in s^-1.
        - 'rate_NI_contact' : float, contact-limit scattering rate in s^-1 (0.0 if U0_contact is None).
        - 'kF' : float, Fermi wavevector in m^-1.
        - 'x_NI' : float, dimensionless argument 2*(kF*a_NI)^2.
        - 'R_BLG' : float, dimensionless BLG chiral transport kernel value.
        - 'R_2DEG' : float, dimensionless conventional 2DEG transport kernel value.

    Raises
    ------
    ValueError
        If any input parameter is None, negative, zero, or non-finite.
    """
    import numpy as np
    from scipy.special import iv

    HBAR = 1.054571817e-34

    if Ek is None or not np.isfinite(Ek) or Ek <= 0:
        raise ValueError("Carrier energy Ek must be a positive finite float.")
    if m_star is None or not np.isfinite(m_star) or m_star <= 0:
        raise ValueError("Effective mass m_star must be a positive finite float.")
    if n_NI is None or not np.isfinite(n_NI) or n_NI < 0:
        raise ValueError("Neutral impurity density n_NI must be a non-negative finite float.")
    if a_NI is None or not np.isfinite(a_NI) or a_NI < 0:
        raise ValueError("Correlation length a_NI must be a non-negative finite float.")
    if V0_NI is None or not np.isfinite(V0_NI):
        raise ValueError("Potential amplitude V0_NI must be a finite float.")

    kF = float(np.sqrt(2.0 * m_star * Ek) / HBAR)
    x_NI = float(2.0 * (kF * a_NI)**2)

    term_blg = 0.5 * iv(0, x_NI) - 0.75 * iv(1, x_NI) + 0.5 * iv(2, x_NI) - 0.25 * iv(3, x_NI)
    R_BLG_val = float(np.exp(-x_NI) * term_blg)

    term_2deg = iv(0, x_NI) - iv(1, x_NI)
    R_2DEG_val = float(np.exp(-x_NI) * term_2deg)

    D_EF = 2.0 * m_star / (np.pi * (HBAR**2))
    U0 = 2.0 * np.pi * (a_NI**2) * V0_NI
    rate_NI = float((2.0 * np.pi * n_NI / HBAR) * (U0**2) * D_EF * R_BLG_val)

    rate_NI_contact = 0.0
    if U0_contact is not None:
        if not np.isfinite(U0_contact) or U0_contact < 0:
            raise ValueError("U0_contact must be a non-negative finite float.")
        rate_NI_contact = float(2.0 * n_NI * m_star * (U0_contact**2) / (HBAR**3))

    return {
        'rate_NI': float(rate_NI),
        'rate_NI_contact': float(rate_NI_contact),
        'kF': float(kF),
        'x_NI': float(x_NI),
        'R_BLG': float(R_BLG_val),
        'R_2DEG': float(R_2DEG_val)
    }
```

---

# Subproblem 2: Inelastic Phonon Scattering Rates with Pauli Blocking

## Prompt
Carrier momentum relaxation in bilayer graphene at finite lattice temperatures is mediated by intrinsic acoustic phonons (AP) and optical phonons (OP). Because phonon absorption and emission change the carrier energy between initial state \( E_k \) and final state \( E_{k'} = E_k \pm \hbar\omega_q \), the transition probability is constrained by the Fermi–Dirac occupation of the final state via the Pauli blocking (PB) ratio \( \frac{1 - f(E_{k'})}{1 - f(E_k)} \), where \( f(E) = \frac{1}{1 + \exp\left(\frac{E - \mu}{k_B T}\right)} \), \( \mu \) is the chemical potential, and \( T \) is the absolute temperature.

For longitudinal acoustic phonons with linear dispersion \( \omega_q = v_{lb} q \), the characteristic acoustic frequency evaluated at the characteristic momentum transfer in BLG is:
\[
\omega_{ba} = \frac{4 \sqrt{2 m^* E_k} v_{lb}}{\pi \hbar}
\]
The equilibrium acoustic phonon occupation follows the Bose–Einstein distribution \( n_{ph} = \frac{1}{\exp\left(\frac{\hbar \omega_{ba}}{k_B T}\right) - 1} \). Integrating over the BLG chiral angular factor yields the acoustic phonon scattering rate with Pauli blocking:
\[
\frac{1}{\tau_{AP}^{PB}} = \frac{V_a^2 m^* \omega_{ba}}{4 \rho_b v_{lb}^2 \hbar^2 (1 - f(E_k))} \left[ n_{ph}(1 - f(E_k + \hbar\omega_{ba})) + (n_{ph} + 1)(1 - f(E_k - \hbar\omega_{ba})) \Theta(E_k - \hbar\omega_{ba}) \right]
\]
where \( V_a \) is the acoustic deformation potential in Joules, \( \rho_b \) is the areal mass density of BLG, \( v_{lb} \) is the longitudinal sound velocity, and \( \Theta(x) \) is the Heaviside step function enforcing the kinematic threshold for phonon emission. In the absence of Pauli blocking, this expression reduces to:
\[
\frac{1}{\tau_{AP}^{noPB}} = \frac{V_a^2 m^* \omega_{ba}}{4 \rho_b v_{lb}^2 \hbar^2} \left[ n_{ph} + (n_{ph} + 1) \Theta(E_k - \hbar\omega_{ba}) \right]
\]

For intrinsic optical phonons, the dispersion is treated as dispersionless with fixed optical phonon energy \( \hbar\omega_o \approx 200 \text{ meV} \), with phonon occupation \( N_o = \frac{1}{\exp\left(\frac{\hbar \omega_o}{k_B T}\right) - 1} \). The optical deformation potential coupling is \( V_o \) (in J/m). The optical phonon scattering rates with and without Pauli blocking are:
\[
\frac{1}{\tau_{OP}^{PB}} = \frac{V_o^2 m^*}{2 \rho_b \omega_o \hbar^2} \left[ N_o \frac{1 - f(E_k + \hbar\omega_o)}{1 - f(E_k)} + (N_o + 1)\Theta(E_k - \hbar\omega_o) \frac{1 - f(E_k - \hbar\omega_o)}{1 - f(E_k)} \right]
\]
\[
\frac{1}{\tau_{OP}^{noPB}} = \frac{V_o^2 m^*}{2 \rho_b \omega_o \hbar^2} \left[ N_o + (N_o + 1)\Theta(E_k - \hbar\omega_o) \right]
\]

Write a Python function `blg_inelastic_phonon_rates` that computes the acoustic and optical phonon scattering rates both with and without Pauli blocking, as well as the characteristic frequencies and occupation numbers.

```python
def blg_inelastic_phonon_rates(Ek, mu, T, m_star=3.0060966e-32, Va=3.0441356e-18, rho_b=15.2e-7, v_lb=2.12e4, hw_o=3.204353268e-20, Vo=4.005441585e-9):
    """
    Compute inelastic acoustic and optical phonon scattering rates with and without Pauli blocking.

    Parameters
    ----------
    Ek : float
        Carrier kinetic energy in Joules (J). Must be a positive finite float.
    mu : float
        Chemical potential in Joules (J). Must be a finite float.
    T : float
        Absolute temperature in Kelvin (K). Must be a positive finite float.
    m_star : float, optional
        Effective electron mass in kg. Default is 3.0060966e-32 kg.
    Va : float, optional
        Acoustic deformation potential in Joules (J). Default is 19.0 eV = 3.0441356e-18 J.
    rho_b : float, optional
        Areal mass density of BLG in kg/m^2. Default is 15.2e-7 kg/m^2.
    v_lb : float, optional
        Longitudinal acoustic phonon velocity in m/s. Default is 2.12e4 m/s.
    hw_o : float, optional
        Optical phonon energy in Joules (J). Default is 200 meV = 3.204353268e-20 J.
    Vo : float, optional
        Optical deformation potential in J/m. Default is 2.5 eV/Angstrom = 4.005441585e-9 J/m.

    Returns
    -------
    dict
        A dictionary containing:
        - 'rate_AP_PB' : float, acoustic phonon scattering rate with Pauli blocking in s^-1.
        - 'rate_AP_noPB' : float, acoustic phonon scattering rate without Pauli blocking in s^-1.
        - 'rate_OP_PB' : float, optical phonon scattering rate with Pauli blocking in s^-1.
        - 'rate_OP_noPB' : float, optical phonon scattering rate without Pauli blocking in s^-1.
        - 'omega_ba' : float, characteristic acoustic phonon frequency in rad/s.
        - 'n_ph_AP' : float, acoustic phonon Bose-Einstein occupation number.
        - 'N_o_OP' : float, optical phonon Bose-Einstein occupation number.

    Raises
    ------
    ValueError
        If any parameter is None, negative, zero, or non-finite.
    """
    pass
```

## Background
In addition to static disorder, charge carriers in bilayer graphene experience inelastic momentum relaxation from lattice vibrations. At low to moderate temperatures, longitudinal acoustic phonons dominate intrinsic electron-phonon scattering. Unlike simplified treatments that assume purely elastic equipartition scattering, retaining the finite acoustic phonon energy exposes a kinematic threshold for phonon emission: an electron must have sufficient initial kinetic energy to emit a phonon of energy \(\hbar \omega_{ba}\).

At higher carrier energies and elevated temperatures, optical phonons become activated. Because the optical phonon branch has a substantial energy scale of roughly 200 meV, optical phonon emission is strictly forbidden for cold carriers below this threshold. When the carrier energy exceeds this threshold, the opening of the emission channel triggers a sharp, multi-order-of-magnitude surge in the scattering rate. Furthermore, in degenerate bilayer graphene where the Fermi level resides inside the conduction band, the Pauli exclusion principle significantly reduces the phase space available for scattering into already-occupied electronic states. Incorporating the ratio of unoccupied final to initial states directly reveals how degenerate screening and thermal broadening dynamically regulate carrier cooling and hot-electron transport.

## Testing Template
```python
def test_case_1():
    import numpy as np
    Ek = 0.1 * 1.602176634e-19
    mu = 0.05 * 1.602176634e-19
    T = 150.0
    out = blg_inelastic_phonon_rates(Ek, mu, T)
    assert np.isclose(out['rate_AP_PB'], 3.637106721906641e+11, rtol=1e-4)
    assert np.isclose(out['rate_AP_noPB'], 3.6492034234896735e+11, rtol=1e-4)
    assert np.isclose(out['omega_ba'], 7944041077518.73, rtol=1e-4)
    assert np.isclose(out['n_ph_AP'], 2.005666662085697, rtol=1e-4)

def test_case_2():
    import numpy as np
    # Optical phonon emission threshold check (threshold at 200 meV = 0.20 eV)
    mu = 0.05 * 1.602176634e-19
    T = 200.0
    out_below = blg_inelastic_phonon_rates(0.15 * 1.602176634e-19, mu, T)
    out_above = blg_inelastic_phonon_rates(0.25 * 1.602176634e-19, mu, T)
    assert out_below['rate_OP_PB'] < 1e7
    assert out_above['rate_OP_PB'] > 1e10
    assert not np.isclose(out_below['rate_OP_PB'], out_above['rate_OP_PB'], rtol=1e-1)

def test_case_3():
    import numpy as np
    # Temperature dependence of phonon occupations and rates
    Ek = 0.1 * 1.602176634e-19
    mu = 0.05 * 1.602176634e-19
    out_cold = blg_inelastic_phonon_rates(Ek, mu, 100.0)
    out_hot = blg_inelastic_phonon_rates(Ek, mu, 300.0)
    assert out_hot['n_ph_AP'] > out_cold['n_ph_AP']
    assert out_hot['rate_AP_PB'] > out_cold['rate_AP_PB']
    assert out_hot['N_o_OP'] > out_cold['N_o_OP']

def test_case_4():
    import numpy as np
    # Pauli blocking effect: rate_AP_PB differs from rate_AP_noPB
    Ek = 0.06 * 1.602176634e-19
    mu = 0.05 * 1.602176634e-19
    T = 150.0
    out = blg_inelastic_phonon_rates(Ek, mu, T)
    assert out['rate_AP_PB'] != out['rate_AP_noPB']
    assert np.isclose(out['rate_AP_PB'], out['rate_AP_noPB'], rtol=0.1)

def test_case_5():
    import numpy as np
    # Emission threshold strictly off when Ek < hw_ba
    Ek = 1e-24 # Extremely low energy where Ek < hw_ba
    mu = 0.0
    T = 150.0
    out = blg_inelastic_phonon_rates(Ek, mu, T)
    assert out['rate_AP_noPB'] > 0.0
    assert out['rate_AP_PB'] > 0.0

def test_case_6():
    import numpy as np
    Ek = 0.1 * 1.602176634e-19
    mu = 0.05 * 1.602176634e-19
    T = 150.0
    out = blg_inelastic_phonon_rates(Ek, mu, T)
    expected_keys = ['rate_AP_PB', 'rate_AP_noPB', 'rate_OP_PB', 'rate_OP_noPB', 'omega_ba', 'n_ph_AP', 'N_o_OP']
    for k in expected_keys:
        assert k in out
        assert isinstance(out[k], float)

def test_case_7():
    import numpy as np
    bad_tuples = [
        (None, 0.05, 150.0),
        (-0.1, 0.05, 150.0),
        (0.1, None, 150.0),
        (0.1, 0.05, -50.0),
        (0.1, 0.05, 0.0)
    ]
    for ek, mu, t in bad_tuples:
        try:
            blg_inelastic_phonon_rates(ek, mu, t)
            assert False, "Should raise for invalid inputs"
        except (ValueError, TypeError):
            pass
```

## Solution
```python
def blg_inelastic_phonon_rates(Ek, mu, T, m_star=3.0060966e-32, Va=3.0441356e-18, rho_b=15.2e-7, v_lb=2.12e4, hw_o=3.204353268e-20, Vo=4.005441585e-9):
    """
    Compute inelastic acoustic and optical phonon scattering rates with and without Pauli blocking.

    Parameters
    ----------
    Ek : float
        Carrier kinetic energy in Joules (J). Must be a positive finite float.
    mu : float
        Chemical potential in Joules (J). Must be a finite float.
    T : float
        Absolute temperature in Kelvin (K). Must be a positive finite float.
    m_star : float, optional
        Effective electron mass in kg. Default is 3.0060966e-32 kg.
    Va : float, optional
        Acoustic deformation potential in Joules (J). Default is 19.0 eV = 3.0441356e-18 J.
    rho_b : float, optional
        Areal mass density of BLG in kg/m^2. Default is 15.2e-7 kg/m^2.
    v_lb : float, optional
        Longitudinal acoustic phonon velocity in m/s. Default is 2.12e4 m/s.
    hw_o : float, optional
        Optical phonon energy in Joules (J). Default is 200 meV = 3.204353268e-20 J.
    Vo : float, optional
        Optical deformation potential in J/m. Default is 2.5 eV/Angstrom = 4.005441585e-9 J/m.

    Returns
    -------
    dict
        A dictionary containing:
        - 'rate_AP_PB' : float, acoustic phonon scattering rate with Pauli blocking in s^-1.
        - 'rate_AP_noPB' : float, acoustic phonon scattering rate without Pauli blocking in s^-1.
        - 'rate_OP_PB' : float, optical phonon scattering rate with Pauli blocking in s^-1.
        - 'rate_OP_noPB' : float, optical phonon scattering rate without Pauli blocking in s^-1.
        - 'omega_ba' : float, characteristic acoustic phonon frequency in rad/s.
        - 'n_ph_AP' : float, acoustic phonon Bose-Einstein occupation number.
        - 'N_o_OP' : float, optical phonon Bose-Einstein occupation number.

    Raises
    ------
    ValueError
        If any parameter is None, negative, zero, or non-finite.
    """
    import numpy as np

    HBAR = 1.054571817e-34
    K_B = 1.380649e-23

    if Ek is None or not np.isfinite(Ek) or Ek <= 0:
        raise ValueError("Ek must be a positive finite float.")
    if mu is None or not np.isfinite(mu):
        raise ValueError("mu must be a finite float.")
    if T is None or not np.isfinite(T) or T <= 0:
        raise ValueError("T must be a positive finite float.")
    if m_star is None or not np.isfinite(m_star) or m_star <= 0:
        raise ValueError("m_star must be a positive finite float.")
    if Va is None or not np.isfinite(Va) or Va <= 0:
        raise ValueError("Va must be a positive finite float.")
    if rho_b is None or not np.isfinite(rho_b) or rho_b <= 0:
        raise ValueError("rho_b must be a positive finite float.")
    if v_lb is None or not np.isfinite(v_lb) or v_lb <= 0:
        raise ValueError("v_lb must be a positive finite float.")
    if hw_o is None or not np.isfinite(hw_o) or hw_o <= 0:
        raise ValueError("hw_o must be a positive finite float.")
    if Vo is None or not np.isfinite(Vo) or Vo <= 0:
        raise ValueError("Vo must be a positive finite float.")

    def fermi_dirac(E):
        arg = (E - mu) / (K_B * T)
        arg = np.clip(arg, -100.0, 100.0)
        return float(1.0 / (1.0 + np.exp(arg)))

    f_Ek = fermi_dirac(Ek)

    # 1. Acoustic Phonon (AP)
    omega_ba = float(4.0 * np.sqrt(2.0 * m_star * Ek) * v_lb / (np.pi * HBAR))
    hw_ba = HBAR * omega_ba

    arg_be_ap = np.clip(hw_ba / (K_B * T), 1e-12, 100.0)
    n_ph = float(1.0 / (np.exp(arg_be_ap) - 1.0))

    f_ap_abs = fermi_dirac(Ek + hw_ba)
    f_ap_em = fermi_dirac(Ek - hw_ba)
    gamma_AP = float((Va**2 * m_star * omega_ba) / (4.0 * rho_b * (v_lb**2) * (HBAR**2)))

    term_ap_abs_PB = n_ph * (1.0 - f_ap_abs)
    term_ap_em_PB = (n_ph + 1.0) * (1.0 - f_ap_em) if (Ek >= hw_ba) else 0.0
    rate_AP_PB = float(gamma_AP * (term_ap_abs_PB + term_ap_em_PB) / (1.0 - f_Ek))
    rate_AP_noPB = float(gamma_AP * (n_ph + ((n_ph + 1.0) if Ek >= hw_ba else 0.0)))

    # 2. Optical Phonon (OP)
    omega_o = float(hw_o / HBAR)
    arg_be_op = np.clip(hw_o / (K_B * T), 1e-12, 100.0)
    No = float(1.0 / (np.exp(arg_be_op) - 1.0))

    f_op_abs = fermi_dirac(Ek + hw_o)
    f_op_em = fermi_dirac(Ek - hw_o)
    gamma_OP = float((Vo**2 * m_star) / (2.0 * rho_b * omega_o * (HBAR**2)))

    term_op_abs_PB = No * (1.0 - f_op_abs)
    term_op_em_PB = (No + 1.0) * (1.0 - f_op_em) if (Ek >= hw_o) else 0.0
    rate_OP_PB = float(gamma_OP * (term_op_abs_PB + term_op_em_PB) / (1.0 - f_Ek))
    rate_OP_noPB = float(gamma_OP * (No + ((No + 1.0) if Ek >= hw_o else 0.0)))

    return {
        'rate_AP_PB': float(rate_AP_PB),
        'rate_AP_noPB': float(rate_AP_noPB),
        'rate_OP_PB': float(rate_OP_PB),
        'rate_OP_noPB': float(rate_OP_noPB),
        'omega_ba': float(omega_ba),
        'n_ph_AP': float(n_ph),
        'N_o_OP': float(No)
    }
```

---

# Main Problem: Unified Carrier Transport Relaxation and Mobility Synthesis in Bilayer Graphene

## Prompt
Accurately predicting charge carrier transport in bilayer graphene devices requires synthesizing elastic disorder scattering and inelastic electron–phonon interactions into an effective transport relaxation time and carrier mobility. In practical suspended and substrate-supported BLG field-effect transistors, the primary momentum-relaxing mechanisms include elastic neutral impurity (NI) scattering, quasi-elastic acoustic phonon (AP) scattering, and threshold-activated optical phonon (OP) scattering.

Assuming independent scattering channels, the total momentum relaxation rate is evaluated via Matthiessen’s rule:
\[
\frac{1}{\tau_{tot}} = \frac{1}{\tau_{NI}} + \frac{1}{\tau_{AP}^{PB}} + \frac{1}{\tau_{OP}^{PB}}
\]
The corresponding carrier momentum relaxation time is \( \tau_{tot} = \left(\frac{1}{\tau_{tot}}\right)^{-1} \). Under the relaxation-time approximation of the Boltzmann transport equation, the energy-dependent carrier drift mobility \( \mu(E_k) \) is given by:
\[
\mu(E_k) = \frac{e \, \tau_{tot}(E_k)}{m^*}
\]
where \( e = 1.602176634 \times 10^{-19} \text{ C} \) is the elementary charge and \( m^* \) is the effective electron mass in BLG. In experimental transport literature, mobilities are conventionally expressed in practical units of \( \text{cm}^2 / (\text{V}\cdot\text{s}) \), where \( 1 \text{ m}^2/(\text{V}\cdot\text{s}) = 10^4 \text{ cm}^2/(\text{V}\cdot\text{s}) \).

To isolate the macroscopic influence of the Pauli exclusion principle on the electron-phonon thermalization dynamics, define the Pauli-blocking suppression factor \( \chi_{PB} \) across the combined inelastic channels:
\[
\chi_{PB} = \frac{\tau_{AP}^{PB-1} + \tau_{OP}^{PB-1}}{\tau_{AP}^{noPB-1} + \tau_{OP}^{noPB-1}}
\]
and the fractional contribution of inelastic scattering to total carrier momentum relaxation:
\[
f_{inelastic} = \frac{\tau_{AP}^{PB-1} + \tau_{OP}^{PB-1}}{\tau_{tot}^{-1}}
\]

Write a Python function `blg_carrier_transport_synthesis` that embeds `blg_elastic_disorder_rates` and `blg_inelastic_phonon_rates` to calculate the total transport scattering rate, total relaxation time, carrier mobility (in \( \text{cm}^2 / (\text{V}\cdot\text{s}) \)), Pauli-blocking ratio, inelastic scattering fraction, and the detailed mechanism breakdown.

```python
def blg_carrier_transport_synthesis(Ek, mu, T, m_star=3.0060966e-32, n_NI=1.0e15, a_NI=0.5e-9, V0_NI=1.602176634e-19, Va=3.0441356e-18, rho_b=15.2e-7, v_lb=2.12e4, hw_o=3.204353268e-20, Vo=4.005441585e-9):
    """
    Synthesize elastic disorder and inelastic phonon scattering to compute total relaxation rate and carrier mobility.

    Parameters
    ----------
    Ek : float
        Carrier kinetic energy in Joules (J). Must be a positive finite float.
    mu : float
        Chemical potential in Joules (J). Must be a finite float.
    T : float
        Absolute temperature in Kelvin (K). Must be a positive finite float.
    m_star : float, optional
        Effective electron mass in kg. Default is 3.0060966e-32 kg.
    n_NI : float, optional
        Sheet density of neutral impurities in m^-2. Default is 1.0e15 m^-2.
    a_NI : float, optional
        Gaussian correlation length of the impurity potential in m. Default is 0.5e-9 m.
    V0_NI : float, optional
        Amplitude of the neutral impurity potential in Joules (J). Default is 1.602176634e-19 J.
    Va : float, optional
        Acoustic deformation potential in Joules (J). Default is 3.0441356e-18 J.
    rho_b : float, optional
        Areal mass density of BLG in kg/m^2. Default is 15.2e-7 kg/m^2.
    v_lb : float, optional
        Longitudinal acoustic phonon velocity in m/s. Default is 2.12e4 m/s.
    hw_o : float, optional
        Optical phonon energy in Joules (J). Default is 3.204353268e-20 J.
    Vo : float, optional
        Optical deformation potential in J/m. Default is 4.005441585e-9 J/m.

    Returns
    -------
    dict
        A dictionary containing:
        - 'rate_tot' : float, total transport momentum relaxation rate in s^-1.
        - 'tau_tot' : float, total transport relaxation time in seconds (s).
        - 'mobility_cm2_Vs' : float, carrier mobility in practical units of cm^2 / (V*s).
        - 'chi_PB' : float, dimensionless Pauli-blocking suppression factor for inelastic channels.
        - 'f_inelastic' : float, dimensionless fraction of total scattering from inelastic channels.
        - 'breakdown' : dict, containing individual scattering rates:
            - 'rate_NI' : float, neutral impurity rate in s^-1.
            - 'rate_AP_PB' : float, acoustic phonon rate with PB in s^-1.
            - 'rate_OP_PB' : float, optical phonon rate with PB in s^-1.

    Raises
    ------
    ValueError
        If any input parameter is None, negative, zero, or non-finite.
    """
    pass
```

## Background
Determining the operational limits of bilayer graphene field-effect transistors requires establishing a unified hierarchy of carrier momentum relaxation mechanisms across operating temperatures and carrier densities. Within the semiclassical Boltzmann transport equation, Matthiessen’s rule enables the concurrent evaluation of distinct physical processes acting on the electronic system. While neutral impurities provide an energy-insensitive, elastic baseline resistivity arising from short-range disorder, phonon scattering introduces a strongly temperature- and energy-dependent inelastic resistive channel.

By synthesizing these mechanisms into macroscopic transport observables, device engineers can assess how quantum-mechanical phenomena dictate practical electrical properties. In particular, computing the carrier mobility directly links microscopic scattering cross-sections to experimental device performance. Because optical phonon emission is kinematically locked until carrier energies exceed 200 meV, low-field transport at cryogenic and intermediate temperatures is predominantly co-governed by neutral impurities and acoustic phonons. Furthermore, tracking the Pauli blocking ratio provides direct physical insight into degenerate state-filling effects, quantifying the degree to which quantum degeneracy shelters carriers from thermal and optical phonon momentum degradation.

## Testing Template
```python
def test_case_1():
    import numpy as np
    Ek = 0.1 * 1.602176634e-19
    mu = 0.05 * 1.602176634e-19
    T = 150.0
    out = blg_carrier_transport_synthesis(Ek, mu, T)
    assert np.isclose(out['rate_tot'], 3.374140090216979e+12, rtol=1e-4)
    assert np.isclose(out['tau_tot'], 2.9637180830144297e-13, rtol=1e-4)
    assert np.isclose(out['mobility_cm2_Vs'], 15795.899004701321, rtol=1e-3)
    assert np.isclose(out['chi_PB'], 0.9966851123349783, rtol=1e-3)
    assert np.isclose(out['f_inelastic'], 0.10779359232417556, rtol=1e-3)

def test_case_2():
    import numpy as np
    # Unit check: ensure mobility is reported in cm^2/(V*s), not m^2/(V*s)
    Ek = 0.1 * 1.602176634e-19
    mu = 0.05 * 1.602176634e-19
    T = 150.0
    out = blg_carrier_transport_synthesis(Ek, mu, T)
    wrong_mobility_m2 = out['mobility_cm2_Vs'] * 1e-4
    assert not np.isclose(out['mobility_cm2_Vs'], wrong_mobility_m2)
    assert out['mobility_cm2_Vs'] > 1e3

def test_case_3():
    import numpy as np
    # High energy behavior: optical phonon activation above 200 meV
    Ek = 0.25 * 1.602176634e-19
    mu = 0.05 * 1.602176634e-19
    T = 200.0
    out = blg_carrier_transport_synthesis(Ek, mu, T)
    assert out['breakdown']['rate_OP_PB'] > 1e10
    assert out['f_inelastic'] > 0.50

def test_case_4():
    import numpy as np
    # Matthiessen's rule consistency check
    Ek = 0.1 * 1.602176634e-19
    mu = 0.05 * 1.602176634e-19
    T = 150.0
    out = blg_carrier_transport_synthesis(Ek, mu, T)
    summed_rate = out['breakdown']['rate_NI'] + out['breakdown']['rate_AP_PB'] + out['breakdown']['rate_OP_PB']
    assert np.isclose(out['rate_tot'], summed_rate, rtol=1e-9)

def test_case_5():
    import numpy as np
    # Physical bounds check
    Ek = 0.08 * 1.602176634e-19
    mu = 0.05 * 1.602176634e-19
    T = 120.0
    out = blg_carrier_transport_synthesis(Ek, mu, T)
    assert 0.0 < out['f_inelastic'] < 1.0
    assert 0.0 < out['chi_PB'] <= 1.0
    assert out['tau_tot'] > 0.0
    assert out['mobility_cm2_Vs'] > 0.0

def test_case_6():
    import numpy as np
    # Verification of breakdown sub-keys
    Ek = 0.1 * 1.602176634e-19
    mu = 0.05 * 1.602176634e-19
    T = 150.0
    out = blg_carrier_transport_synthesis(Ek, mu, T)
    assert 'breakdown' in out
    for subkey in ['rate_NI', 'rate_AP_PB', 'rate_OP_PB']:
        assert subkey in out['breakdown']
        assert isinstance(out['breakdown'][subkey], float)

def test_case_7():
    import numpy as np
    bad_inputs = [None, -0.05, 0.0, float('nan')]
    for bad in bad_inputs:
        try:
            blg_carrier_transport_synthesis(bad, 0.05 * 1.602176634e-19, 150.0)
            assert False, "Should raise on bad energy input"
        except (ValueError, TypeError):
            pass
```

## Solution
```python
def blg_elastic_disorder_rates(Ek, m_star=3.0060966e-32, n_NI=1.0e15, a_NI=0.5e-9, V0_NI=1.602176634e-19, U0_contact=None):
    """
    Compute elastic neutral impurity scattering rates and transport kernels in bilayer graphene.

    Parameters
    ----------
    Ek : float
        Carrier kinetic energy in Joules (J). Must be a positive finite float.
    m_star : float, optional
        Effective electron mass in kilograms (kg). Default is 0.033 * m_e = 3.0060966e-32 kg.
    n_NI : float, optional
        Sheet density of neutral impurities in m^-2. Default is 1.0e15 m^-2 (10^11 cm^-2).
    a_NI : float, optional
        Gaussian correlation length of the impurity potential in meters (m). Default is 0.5e-9 m (0.5 nm).
    V0_NI : float, optional
        Amplitude of the neutral impurity potential in Joules (J). Default is 1.0 eV = 1.602176634e-19 J.
    U0_contact : float or None, optional
        Integrated contact potential strength in J*m^2. If provided, evaluates contact-limit rate. Default is None.

    Returns
    -------
    dict
        A dictionary containing:
        - 'rate_NI' : float, Gaussian neutral impurity scattering rate in s^-1.
        - 'rate_NI_contact' : float, contact-limit scattering rate in s^-1 (0.0 if U0_contact is None).
        - 'kF' : float, Fermi wavevector in m^-1.
        - 'x_NI' : float, dimensionless argument 2*(kF*a_NI)^2.
        - 'R_BLG' : float, dimensionless BLG chiral transport kernel value.
        - 'R_2DEG' : float, dimensionless conventional 2DEG transport kernel value.

    Raises
    ------
    ValueError
        If any input parameter is None, negative, zero, or non-finite.
    """
    import numpy as np
    from scipy.special import iv

    HBAR = 1.054571817e-34

    if Ek is None or not np.isfinite(Ek) or Ek <= 0:
        raise ValueError("Carrier energy Ek must be a positive finite float.")
    if m_star is None or not np.isfinite(m_star) or m_star <= 0:
        raise ValueError("Effective mass m_star must be a positive finite float.")
    if n_NI is None or not np.isfinite(n_NI) or n_NI < 0:
        raise ValueError("Neutral impurity density n_NI must be a non-negative finite float.")
    if a_NI is None or not np.isfinite(a_NI) or a_NI < 0:
        raise ValueError("Correlation length a_NI must be a non-negative finite float.")
    if V0_NI is None or not np.isfinite(V0_NI):
        raise ValueError("Potential amplitude V0_NI must be a finite float.")

    kF = float(np.sqrt(2.0 * m_star * Ek) / HBAR)
    x_NI = float(2.0 * (kF * a_NI)**2)

    term_blg = 0.5 * iv(0, x_NI) - 0.75 * iv(1, x_NI) + 0.5 * iv(2, x_NI) - 0.25 * iv(3, x_NI)
    R_BLG_val = float(np.exp(-x_NI) * term_blg)

    term_2deg = iv(0, x_NI) - iv(1, x_NI)
    R_2DEG_val = float(np.exp(-x_NI) * term_2deg)

    D_EF = 2.0 * m_star / (np.pi * (HBAR**2))
    U0 = 2.0 * np.pi * (a_NI**2) * V0_NI
    rate_NI = float((2.0 * np.pi * n_NI / HBAR) * (U0**2) * D_EF * R_BLG_val)

    rate_NI_contact = 0.0
    if U0_contact is not None:
        if not np.isfinite(U0_contact) or U0_contact < 0:
            raise ValueError("U0_contact must be a non-negative finite float.")
        rate_NI_contact = float(2.0 * n_NI * m_star * (U0_contact**2) / (HBAR**3))

    return {
        'rate_NI': float(rate_NI),
        'rate_NI_contact': float(rate_NI_contact),
        'kF': float(kF),
        'x_NI': float(x_NI),
        'R_BLG': float(R_BLG_val),
        'R_2DEG': float(R_2DEG_val)
    }

def blg_inelastic_phonon_rates(Ek, mu, T, m_star=3.0060966e-32, Va=3.0441356e-18, rho_b=15.2e-7, v_lb=2.12e4, hw_o=3.204353268e-20, Vo=4.005441585e-9):
    """
    Compute inelastic acoustic and optical phonon scattering rates with and without Pauli blocking.

    Parameters
    ----------
    Ek : float
        Carrier kinetic energy in Joules (J). Must be a positive finite float.
    mu : float
        Chemical potential in Joules (J). Must be a finite float.
    T : float
        Absolute temperature in Kelvin (K). Must be a positive finite float.
    m_star : float, optional
        Effective electron mass in kg. Default is 3.0060966e-32 kg.
    Va : float, optional
        Acoustic deformation potential in Joules (J). Default is 19.0 eV = 3.0441356e-18 J.
    rho_b : float, optional
        Areal mass density of BLG in kg/m^2. Default is 15.2e-7 kg/m^2.
    v_lb : float, optional
        Longitudinal acoustic phonon velocity in m/s. Default is 2.12e4 m/s.
    hw_o : float, optional
        Optical phonon energy in Joules (J). Default is 200 meV = 3.204353268e-20 J.
    Vo : float, optional
        Optical deformation potential in J/m. Default is 2.5 eV/Angstrom = 4.005441585e-9 J/m.

    Returns
    -------
    dict
        A dictionary containing:
        - 'rate_AP_PB' : float, acoustic phonon scattering rate with Pauli blocking in s^-1.
        - 'rate_AP_noPB' : float, acoustic phonon scattering rate without Pauli blocking in s^-1.
        - 'rate_OP_PB' : float, optical phonon scattering rate with Pauli blocking in s^-1.
        - 'rate_OP_noPB' : float, optical phonon scattering rate without Pauli blocking in s^-1.
        - 'omega_ba' : float, characteristic acoustic phonon frequency in rad/s.
        - 'n_ph_AP' : float, acoustic phonon Bose-Einstein occupation number.
        - 'N_o_OP' : float, optical phonon Bose-Einstein occupation number.

    Raises
    ------
    ValueError
        If any parameter is None, negative, zero, or non-finite.
    """
    import numpy as np

    HBAR = 1.054571817e-34
    K_B = 1.380649e-23

    if Ek is None or not np.isfinite(Ek) or Ek <= 0:
        raise ValueError("Ek must be a positive finite float.")
    if mu is None or not np.isfinite(mu):
        raise ValueError("mu must be a finite float.")
    if T is None or not np.isfinite(T) or T <= 0:
        raise ValueError("T must be a positive finite float.")
    if m_star is None or not np.isfinite(m_star) or m_star <= 0:
        raise ValueError("m_star must be a positive finite float.")
    if Va is None or not np.isfinite(Va) or Va <= 0:
        raise ValueError("Va must be a positive finite float.")
    if rho_b is None or not np.isfinite(rho_b) or rho_b <= 0:
        raise ValueError("rho_b must be a positive finite float.")
    if v_lb is None or not np.isfinite(v_lb) or v_lb <= 0:
        raise ValueError("v_lb must be a positive finite float.")
    if hw_o is None or not np.isfinite(hw_o) or hw_o <= 0:
        raise ValueError("hw_o must be a positive finite float.")
    if Vo is None or not np.isfinite(Vo) or Vo <= 0:
        raise ValueError("Vo must be a positive finite float.")

    def fermi_dirac(E):
        arg = (E - mu) / (K_B * T)
        arg = np.clip(arg, -100.0, 100.0)
        return float(1.0 / (1.0 + np.exp(arg)))

    f_Ek = fermi_dirac(Ek)

    # 1. Acoustic Phonon (AP)
    omega_ba = float(4.0 * np.sqrt(2.0 * m_star * Ek) * v_lb / (np.pi * HBAR))
    hw_ba = HBAR * omega_ba

    arg_be_ap = np.clip(hw_ba / (K_B * T), 1e-12, 100.0)
    n_ph = float(1.0 / (np.exp(arg_be_ap) - 1.0))

    f_ap_abs = fermi_dirac(Ek + hw_ba)
    f_ap_em = fermi_dirac(Ek - hw_ba)
    gamma_AP = float((Va**2 * m_star * omega_ba) / (4.0 * rho_b * (v_lb**2) * (HBAR**2)))

    term_ap_abs_PB = n_ph * (1.0 - f_ap_abs)
    term_ap_em_PB = (n_ph + 1.0) * (1.0 - f_ap_em) if (Ek >= hw_ba) else 0.0
    rate_AP_PB = float(gamma_AP * (term_ap_abs_PB + term_ap_em_PB) / (1.0 - f_Ek))
    rate_AP_noPB = float(gamma_AP * (n_ph + ((n_ph + 1.0) if Ek >= hw_ba else 0.0)))

    # 2. Optical Phonon (OP)
    omega_o = float(hw_o / HBAR)
    arg_be_op = np.clip(hw_o / (K_B * T), 1e-12, 100.0)
    No = float(1.0 / (np.exp(arg_be_op) - 1.0))

    f_op_abs = fermi_dirac(Ek + hw_o)
    f_op_em = fermi_dirac(Ek - hw_o)
    gamma_OP = float((Vo**2 * m_star) / (2.0 * rho_b * omega_o * (HBAR**2)))

    term_op_abs_PB = No * (1.0 - f_op_abs)
    term_op_em_PB = (No + 1.0) * (1.0 - f_op_em) if (Ek >= hw_o) else 0.0
    rate_OP_PB = float(gamma_OP * (term_op_abs_PB + term_op_em_PB) / (1.0 - f_Ek))
    rate_OP_noPB = float(gamma_OP * (No + ((No + 1.0) if Ek >= hw_o else 0.0)))

    return {
        'rate_AP_PB': float(rate_AP_PB),
        'rate_AP_noPB': float(rate_AP_noPB),
        'rate_OP_PB': float(rate_OP_PB),
        'rate_OP_noPB': float(rate_OP_noPB),
        'omega_ba': float(omega_ba),
        'n_ph_AP': float(n_ph),
        'N_o_OP': float(No)
    }

def blg_carrier_transport_synthesis(Ek, mu, T, m_star=3.0060966e-32, n_NI=1.0e15, a_NI=0.5e-9, V0_NI=1.602176634e-19, Va=3.0441356e-18, rho_b=15.2e-7, v_lb=2.12e4, hw_o=3.204353268e-20, Vo=4.005441585e-9):
    """
    Synthesize elastic disorder and inelastic phonon scattering to compute total relaxation rate and carrier mobility.

    Parameters
    ----------
    Ek : float
        Carrier kinetic energy in Joules (J). Must be a positive finite float.
    mu : float
        Chemical potential in Joules (J). Must be a finite float.
    T : float
        Absolute temperature in Kelvin (K). Must be a positive finite float.
    m_star : float, optional
        Effective electron mass in kg. Default is 3.0060966e-32 kg.
    n_NI : float, optional
        Sheet density of neutral impurities in m^-2. Default is 1.0e15 m^-2.
    a_NI : float, optional
        Gaussian correlation length of the impurity potential in m. Default is 0.5e-9 m.
    V0_NI : float, optional
        Amplitude of the neutral impurity potential in Joules (J). Default is 1.602176634e-19 J.
    Va : float, optional
        Acoustic deformation potential in Joules (J). Default is 3.0441356e-18 J.
    rho_b : float, optional
        Areal mass density of BLG in kg/m^2. Default is 15.2e-7 kg/m^2.
    v_lb : float, optional
        Longitudinal acoustic phonon velocity in m/s. Default is 2.12e4 m/s.
    hw_o : float, optional
        Optical phonon energy in Joules (J). Default is 3.204353268e-20 J.
    Vo : float, optional
        Optical deformation potential in J/m. Default is 4.005441585e-9 J/m.

    Returns
    -------
    dict
        A dictionary containing:
        - 'rate_tot' : float, total transport momentum relaxation rate in s^-1.
        - 'tau_tot' : float, total transport relaxation time in seconds (s).
        - 'mobility_cm2_Vs' : float, carrier mobility in practical units of cm^2 / (V*s).
        - 'chi_PB' : float, dimensionless Pauli-blocking suppression factor for inelastic channels.
        - 'f_inelastic' : float, dimensionless fraction of total scattering from inelastic channels.
        - 'breakdown' : dict, containing individual scattering rates:
            - 'rate_NI' : float, neutral impurity rate in s^-1.
            - 'rate_AP_PB' : float, acoustic phonon rate with PB in s^-1.
            - 'rate_OP_PB' : float, optical phonon rate with PB in s^-1.

    Raises
    ------
    ValueError
        If any input parameter is None, negative, zero, or non-finite.
    """
    import numpy as np

    E_CHARGE = 1.602176634e-19

    if Ek is None or not np.isfinite(Ek) or Ek <= 0:
        raise ValueError("Ek must be a positive finite float.")
    if mu is None or not np.isfinite(mu):
        raise ValueError("mu must be a finite float.")
    if T is None or not np.isfinite(T) or T <= 0:
        raise ValueError("T must be a positive finite float.")
    if m_star is None or not np.isfinite(m_star) or m_star <= 0:
        raise ValueError("m_star must be a positive finite float.")

    res_el = blg_elastic_disorder_rates(Ek=Ek, m_star=m_star, n_NI=n_NI, a_NI=a_NI, V0_NI=V0_NI)
    res_ph = blg_inelastic_phonon_rates(Ek=Ek, mu=mu, T=T, m_star=m_star, Va=Va, rho_b=rho_b, v_lb=v_lb, hw_o=hw_o, Vo=Vo)

    rate_NI = res_el['rate_NI']
    rate_AP_PB = res_ph['rate_AP_PB']
    rate_AP_noPB = res_ph['rate_AP_noPB']
    rate_OP_PB = res_ph['rate_OP_PB']
    rate_OP_noPB = res_ph['rate_OP_noPB']

    rate_tot = float(rate_NI + rate_AP_PB + rate_OP_PB)
    tau_tot = float(1.0 / rate_tot)

    # Mobility in cm^2 / (V * s): 1 m^2/(V*s) = 1e4 cm^2/(V*s)
    mu_SI = E_CHARGE * tau_tot / m_star
    mobility_cm2_Vs = float(mu_SI * 1e4)

    rate_inelastic_PB = rate_AP_PB + rate_OP_PB
    rate_inelastic_noPB = rate_AP_noPB + rate_OP_noPB
    chi_PB = float(rate_inelastic_PB / rate_inelastic_noPB) if rate_inelastic_noPB > 0 else 1.0
    f_inelastic = float(rate_inelastic_PB / rate_tot)

    return {
        'rate_tot': float(rate_tot),
        'tau_tot': float(tau_tot),
        'mobility_cm2_Vs': float(mobility_cm2_Vs),
        'chi_PB': float(chi_PB),
        'f_inelastic': float(f_inelastic),
        'breakdown': {
            'rate_NI': float(rate_NI),
            'rate_AP_PB': float(rate_AP_PB),
            'rate_OP_PB': float(rate_OP_PB)
        }
    }
```
