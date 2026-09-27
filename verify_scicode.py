import numpy as np
from scipy.special import iv

HBAR = 1.054571817e-34
M_E = 9.1093837e-31
M_STAR_DEFAULT = 0.033 * M_E
E_CHARGE = 1.602176634e-19
K_B = 1.380649e-23

# ==========================================
# SP1 SOLUTION
# ==========================================
def blg_elastic_disorder_rates(Ek, m_star=M_STAR_DEFAULT, n_NI=1e15, a_NI=0.5e-9, V0_NI=1.0*E_CHARGE, U0_contact=None):
    if Ek is None or not np.isfinite(Ek) or Ek <= 0:
        raise ValueError("Carrier energy Ek must be a positive finite float.")
    if m_star is None or not np.isfinite(m_star) or m_star <= 0:
        raise ValueError("Effective mass m_star must be a positive finite float.")
    if n_NI is None or not np.isfinite(n_NI) or n_NI < 0:
        raise ValueError("Neutral impurity density n_NI must be non-negative finite float.")
    if a_NI is None or not np.isfinite(a_NI) or a_NI < 0:
        raise ValueError("Correlation length a_NI must be non-negative finite float.")
    if V0_NI is None or not np.isfinite(V0_NI):
        raise ValueError("Potential amplitude V0_NI must be finite float.")
    
    kF = float(np.sqrt(2.0 * m_star * Ek) / HBAR)
    x_NI = float(2.0 * (kF * a_NI)**2)
    
    term = 0.5 * iv(0, x_NI) - 0.75 * iv(1, x_NI) + 0.5 * iv(2, x_NI) - 0.25 * iv(3, x_NI)
    R_BLG_val = float(np.exp(-x_NI) * term)
    R_2DEG_val = float(np.exp(-x_NI) * (iv(0, x_NI) - iv(1, x_NI)))
    
    D_EF = 2.0 * m_star / (np.pi * (HBAR**2))
    U0 = 2.0 * np.pi * (a_NI**2) * V0_NI
    rate_NI = float((2.0 * np.pi * n_NI / HBAR) * (U0**2) * D_EF * R_BLG_val)
    
    rate_NI_contact = 0.0
    if U0_contact is not None:
        if not np.isfinite(U0_contact) or U0_contact < 0:
            raise ValueError("U0_contact must be non-negative finite float.")
        rate_NI_contact = float(2.0 * n_NI * m_star * (U0_contact**2) / (HBAR**3))
        
    return {
        'rate_NI': float(rate_NI),
        'rate_NI_contact': float(rate_NI_contact),
        'kF': float(kF),
        'x_NI': float(x_NI),
        'R_BLG': float(R_BLG_val),
        'R_2DEG': float(R_2DEG_val)
    }

# ==========================================
# SP2 SOLUTION
# ==========================================
def blg_inelastic_phonon_rates(Ek, mu, T, m_star=M_STAR_DEFAULT, Va=19.0*E_CHARGE, rho_b=15.2e-7, v_lb=2.12e4, hw_o=0.200*E_CHARGE, Vo=2.5e10*E_CHARGE):
    if Ek is None or not np.isfinite(Ek) or Ek <= 0:
        raise ValueError("Ek must be a positive finite float.")
    if mu is None or not np.isfinite(mu):
        raise ValueError("mu must be a finite float.")
    if T is None or not np.isfinite(T) or T <= 0:
        raise ValueError("T must be a positive finite float.")
    if m_star is None or not np.isfinite(m_star) or m_star <= 0:
        raise ValueError("m_star must be positive.")
    if Va is None or not np.isfinite(Va) or Va <= 0:
        raise ValueError("Va must be positive.")
    if rho_b is None or not np.isfinite(rho_b) or rho_b <= 0:
        raise ValueError("rho_b must be positive.")
    if v_lb is None or not np.isfinite(v_lb) or v_lb <= 0:
        raise ValueError("v_lb must be positive.")
    if hw_o is None or not np.isfinite(hw_o) or hw_o <= 0:
        raise ValueError("hw_o must be positive.")
    if Vo is None or not np.isfinite(Vo) or Vo <= 0:
        raise ValueError("Vo must be positive.")

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

# ==========================================
# MP SOLUTION
# ==========================================
def blg_carrier_transport_synthesis(Ek, mu, T, m_star=M_STAR_DEFAULT, n_NI=1e15, a_NI=0.5e-9, V0_NI=1.0*E_CHARGE, Va=19.0*E_CHARGE, rho_b=15.2e-7, v_lb=2.12e4, hw_o=0.200*E_CHARGE, Vo=2.5e10*E_CHARGE):
    if Ek is None or not np.isfinite(Ek) or Ek <= 0:
        raise ValueError("Ek must be a positive finite float.")
    if mu is None or not np.isfinite(mu):
        raise ValueError("mu must be a finite float.")
    if T is None or not np.isfinite(T) or T <= 0:
        raise ValueError("T must be a positive finite float.")
    if m_star is None or not np.isfinite(m_star) or m_star <= 0:
        raise ValueError("m_star must be positive.")
        
    res_el = blg_elastic_disorder_rates(Ek=Ek, m_star=m_star, n_NI=n_NI, a_NI=a_NI, V0_NI=V0_NI)
    res_ph = blg_inelastic_phonon_rates(Ek=Ek, mu=mu, T=T, m_star=m_star, Va=Va, rho_b=rho_b, v_lb=v_lb, hw_o=hw_o, Vo=Vo)
    
    rate_NI = res_el['rate_NI']
    rate_AP_PB = res_ph['rate_AP_PB']
    rate_AP_noPB = res_ph['rate_AP_noPB']
    rate_OP_PB = res_ph['rate_OP_PB']
    rate_OP_noPB = res_ph['rate_OP_noPB']
    
    rate_tot = float(rate_NI + rate_AP_PB + rate_OP_PB)
    tau_tot = float(1.0 / rate_tot)
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

print("Running SP1 test cases...")
def test_case_1():
    import numpy as np
    out = blg_elastic_disorder_rates(0.1 * 1.602176634e-19)
    assert np.isclose(out['rate_NI'], 3.010429408887473e+12, rtol=1e-5)
    assert np.isclose(out['kF'], 294303550.5772216, rtol=1e-5)
    assert np.isclose(out['R_BLG'], 0.46358936248662097, rtol=1e-5)
test_case_1()

def test_case_2():
    import numpy as np
    # Discriminative: Chiral BLG vs conventional 2DEG kernel
    out = blg_elastic_disorder_rates(0.1 * 1.602176634e-19)
    assert not np.isclose(out['R_BLG'], out['R_2DEG'], rtol=1e-2)
    assert out['R_BLG'] < out['R_2DEG']
test_case_2()

def test_case_3():
    import numpy as np
    # Short-range contact limit test
    out = blg_elastic_disorder_rates(0.05 * 1.602176634e-19, a_NI=1e-12, U0_contact=1.57e-37)
    expected_contact = 2.0 * 1e15 * (0.033 * 9.1093837e-31) * (1.57e-37**2) / (1.054571817e-34**3)
    assert np.isclose(out['rate_NI_contact'], expected_contact, rtol=1e-5)
    assert np.isclose(out['R_BLG'], 0.5, atol=1e-3)
test_case_3()

def test_case_4():
    import numpy as np
    # Bundled invalid inputs
    bad_inputs = [None, -1.0, 0.0, float('nan'), float('inf')]
    for bad in bad_inputs:
        try:
            blg_elastic_disorder_rates(bad)
            assert False, "Should raise for bad Ek"
        except (ValueError, TypeError):
            pass
test_case_4()

print("SP1 tests passed.")

print("Running SP2 test cases...")
def test_case_5():
    import numpy as np
    out = blg_inelastic_phonon_rates(0.1 * 1.602176634e-19, 0.05 * 1.602176634e-19, 150.0)
    assert np.isclose(out['rate_AP_PB'], 3.637106721906641e+11, rtol=1e-4)
    assert np.isclose(out['rate_AP_noPB'], 3.6492034234896735e+11, rtol=1e-4)
    assert out['rate_OP_PB'] < 1e5  # Highly suppressed below 200 meV
test_case_5()

def test_case_6():
    import numpy as np
    # Optical phonon emission threshold test (200 meV threshold)
    out_below = blg_inelastic_phonon_rates(0.15 * 1.602176634e-19, 0.05 * 1.602176634e-19, 200.0)
    out_above = blg_inelastic_phonon_rates(0.25 * 1.602176634e-19, 0.05 * 1.602176634e-19, 200.0)
    assert out_below['rate_OP_PB'] < 1e7
    assert out_above['rate_OP_PB'] > 1e10
    assert not np.isclose(out_below['rate_OP_PB'], out_above['rate_OP_PB'], rtol=1e-1)
test_case_6()

def test_case_7():
    import numpy as np
    bad_inputs = [(None, 0.05, 150.0), (0.1, None, 150.0), (0.1, 0.05, -10.0), (0.1, 0.05, 0.0)]
    for ek_b, mu_b, t_b in bad_inputs:
        try:
            blg_inelastic_phonon_rates(ek_b, mu_b, t_b)
            assert False, "Should raise on bad parameters"
        except (ValueError, TypeError):
            pass
test_case_7()

print("SP2 tests passed.")

print("Running MP test cases...")
def test_case_8():
    import numpy as np
    out = blg_carrier_transport_synthesis(0.1 * 1.602176634e-19, 0.05 * 1.602176634e-19, 150.0)
    assert np.isclose(out['rate_tot'], 3.374140090216979e+12, rtol=1e-4)
    assert np.isclose(out['mobility_cm2_Vs'], 15795.899004701321, rtol=1e-3)
    assert 0.95 < out['chi_PB'] <= 1.0
    assert 0.05 < out['f_inelastic'] < 0.20
test_case_8()

def test_case_9():
    import numpy as np
    # Regression test at higher energy 0.25 eV where optical emission turns on
    out = blg_carrier_transport_synthesis(0.25 * 1.602176634e-19, 0.05 * 1.602176634e-19, 200.0)
    assert out['breakdown']['rate_OP_PB'] > 1e10
    assert out['mobility_cm2_Vs'] > 0.0
    # Negative test ensuring mobility is in cm^2/(V*s), not m^2/(V*s)
    wrong_mobility_m2 = out['mobility_cm2_Vs'] * 1e-4
    assert not np.isclose(out['mobility_cm2_Vs'], wrong_mobility_m2)
test_case_9()

print("All tests passed completely!")
