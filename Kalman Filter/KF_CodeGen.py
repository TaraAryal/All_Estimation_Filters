# File name:		KF_CodeGen.py
# Written by:		Niranjan Bhujel
# Date:			    28-December-2021
# Description:	    File to generate C code for KF


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
    """
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]

    u = u[0]

    Ki = 2.0
    Rp = 0.05
    Tg = 0.2
    M = 4.0
    D = 1.5

    return ca.vertcat(
        x2,
        x3,
        -Ki/(M*Tg)*x1-(D/(M*Tg)+1/(Rp*M*Tg))*x2-(D/M+1/Tg)*x3-1/(M*Tg)*u,
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

# A, B, C and D matrices
x = ca.SX.sym('x', 3)
u = ca.SX.sym('u', 1)

dotx = Fc(x, u)
y = Hc(x)

A = ca.jacobian(dotx, x)
B = ca.jacobian(dotx, u)

C = ca.jacobian(y, x)
D = ca.jacobian(y, u)


# Components of process noise covariances
Q11 = ca.SX.sym('Q11')
Q22 = ca.SX.sym('Q22')
Q33 = ca.SX.sym('Q33')

Qw = BasicUtils.directSum([Q11, Q22, Q33])

# Components of mesurement noise covariances
R = ca.SX.sym('R')
Rv = BasicUtils.directSum([R])

Ts = 0.02

In, Out, InName, OutName = Estimator.KF(A, B, C, D, Qw, Rv, Ts)

KF_Func = ca.Function(
    'KF_Func',
    In + [Q11, Q22, Q33] + [R],
    Out,
    InName + ['Q11', 'Q22', 'Q33'] + ['R'],
    OutName
)
BasicUtils.Gen_Code(
    func=KF_Func,
    filename='KF_Code',
    mex=False,
    printhelp=True
)
