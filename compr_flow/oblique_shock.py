import math


def fun_theta_beta_m(mach: float, beta:float, gamma=1.4) -> float:
    """

    Calculate theta-beta-mach equation

    Parameters
    ----------
    mach:float
        Mach number
    beta: float
        oblique shock angle
    gamma: float
        isentropic exponent

    Returns
    -------
    theta: float
        tangent of an corner angle
    """
    gp1 = gamma + 1.0

    mult1 = mach * mach * math.sin(beta) * math.sin(beta) - 1.0
    mult2 = 1.0 + (gp1 / 2.0 - math.sin(beta)*math.sin(beta)) * mach*mach

    return mult1 / mult2 / math.tan(beta)


def calc_mach_os(mach: float, beta:float, theta:float, gamma=1.4) -> float:
    """

    Calculate mach number after oblique shock

    Parameters
    ----------
    mach:float
        Mach number
    beta: float
        oblique shock angle
    theta: float
        tangent of an corner angle
    gamma: float
        isentropic exponent

    Returns
    -------
    mach2: float
        Mach number after an oblique shock
    """
    gm1_2 = (gamma - 1.0) / 2.0
    m_sin_2 = mach*mach * math.sin(beta)*math.sin(beta)

    mult1 = 1.0 + gm1_2 * m_sin_2
    mult2 = gamma * m_sin_2 - gm1_2

    return math.sqrt(mult1 / mult2) / math.sin(beta - theta)


def rise_pres(mach: float, beta:float, gamma=1.4) -> float:
    """

    Calculate rise in pressure after oblique shock

    Parameters
    ----------
    mach:float
        Mach number
    beta: float
        oblique shock angle
    gamma: float
        isentropic exponent

    Returns
    -------
    mach2: float
        Mach number after an oblique shock
    """
    m_sin_2 = mach * mach * math.sin(beta) * math.sin(beta)

    return (2.0 * gamma / (gamma + 1.0) * m_sin_2) - ((gamma-1.0) / (gamma + 1.0))