import math


def pi(var_lambda: float, gamma=1.4) -> float:
    """"
    Returns gas dynamic function pi

    Parameters
    ----------
    var_lambda: float
        superficial velocity of flow
    gamma: float
        isentropic exponent (default 1.4)

    Returns
    -------
    pi: float
        values of GD function
    """
    gp1 = gamma + 1.0
    gm1 = gamma - 1.0

    power = gamma / gm1

    return math.pow(1.0 - (gm1 / gp1) * var_lambda*var_lambda, power)


def tau(var_lambda: float, gamma=1.4) -> float:
    """"
    Returns gas dynamic function tau

    Parameters
    ----------
    var_lambda: float
        superficial velocity of flow
    gamma: float
        isentropic exponent (default 1.4)

    Returns
    -------
    tau: float
        values of GD function
    """
    gp1 = gamma + 1.0
    gm1 = gamma - 1.0

    return 1.0 - (gm1 / gp1) * var_lambda*var_lambda


def q(var_lambda: float, gamma=1.4) -> float:
    """"
    Returns gas dynamic function q

    Parameters
    ----------
    var_lambda: float
        superficial velocity of flow
    gamma: float
        isentropic exponent (default 1.4)

    Returns
    -------
    q: float
        values of GD function
    """
    gp1 = gamma + 1.0
    gm1 = gamma - 1.0

    power = 1.0 / gm1

    mult1 = math.pow(gp1 / 2.0, power)
    mult2 = math.pow(1.0 - (gm1 / gp1) * var_lambda*var_lambda, power)

    return mult1 * mult2 * var_lambda


def calc_mach(var_lambda: float, gamma=1.4) -> float:
    """"
    Returns Mach number instead superficial velocity

    Parameters
    ----------
    var_lambda: float
        superficial velocity of flow
    gamma: float
        isentropic exponent (default 1.4)

    Returns
    -------
    M: float
        Mach number
    """
    gp1 = gamma + 1.0
    gm1 = gamma - 1.0

    mult1 = var_lambda*var_lambda * 2.0
    mult2 = gp1 - gm1 * var_lambda*var_lambda

    return math.sqrt(mult1 / mult2)


def calc_lambda(var_mach: float, gamma=1.4) -> float:
    """"
    Returns superficial velocity instead Mach number

    Parameters
    ----------
    var_mach: float
        Mach number
    gamma: float
        isentropic exponent (default 1.4)

    Returns
    -------
    lambda: float
        superficial velocity
    """
    gp1 = gamma + 1.0
    gm1 = gamma - 1.0

    mult1 = var_lambda*var_lambda * 2.0
    mult2 = gp1 - gm1 * var_lambda*var_lambda

    return math.sqrt(mult1 / mult2)
