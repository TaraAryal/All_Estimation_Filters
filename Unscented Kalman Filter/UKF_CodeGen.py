# File name:		UKF_CodeGen.py
# Written by:		Niranjan Bhujel
# Date:			    28-December-2021
# Description:	    File to generate C code for UKF


import casadi as ca
from pynlcontrol import Estimator, BasicUtils


def Fc(x, u):
    """
    Function that returns right hand side of state equations

    Parameters
    ----------
    x : ca.SX.sym array
        State vector
    u : ca.SX.sym array
        Control input vector

    Returns
    -------
    ca.SX.sym
        Right hand side of state equation

    Parameters to be estimated are also augmented as states whose derivative is zero 
    """
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    M = x[3]
    D = x[4]

    u = u[0]

    Ki = 2.0
    Rp = 0.05
    Tg = 0.2

    return ca.vertcat(
        x2,
        x3,
        -Ki/(M*Tg)*x1-(D/(M*Tg)+1/(Rp*M*Tg))*x2-(D/M+1/Tg)*x3-1/(M*Tg)*u,
        0,
        0
    )


def Hc(x):
    """
    Function that returns measurement variable in terms of x.

    Parameters
    ----------
    x : ca.SX.sym array
        State vector

    Returns
    -------
    ca.SX.sym
        Measurement variable

    Measurement in this case is $\Delta\omega$ only which is second state variable.
    """
    return x[1]


Q11 = ca.SX.sym('Q11')
Q22 = ca.SX.sym('Q22')
Q33 = ca.SX.sym('Q33')
Q44 = ca.SX.sym('Q44')
Q55 = ca.SX.sym('Q55')

Qw = BasicUtils.directSum([Q11, Q22, Q33, Q44, Q55])

R = ca.SX.sym('R')
Rv = BasicUtils.directSum([R])

Ts = 0.02
In, Out, InName, OutName = Estimator.UKF(5, 1, 1, Fc, Hc, Qw, Rv, Ts)

UKF_Func = ca.Function(
    'UKF_Func',
    In + [Q11, Q22, Q33, Q44, Q55] + [R],
    Out,
    InName + ['Q11', 'Q22', 'Q33', 'Q44', 'Q55'] + ['R'],
    OutName
)

BasicUtils.Gen_Code(
    func=UKF_Func,
    filename='UKF_Code',
    mex=False,
    printhelp=True
)
