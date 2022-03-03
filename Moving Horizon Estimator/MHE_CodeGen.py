# File name:		MHE_CodeGen.py
# Written by:		Niranjan Bhujel
# Date:			    27-February-2022
# Description:	    File to generate C code for MHE with arrival cost


from pynlcontrol import BasicUtils, Estimator
import casadi as ca


def Fc(x, u, p):
    """
    Function that returns right hand side of state equations

    Parameters
    ----------
    x : ca.SX.sym array
        State vector
    u : ca.SX.sym array
        Control input vector
    p : ca.SX.sym array
        Unknown parameters to be estimated

    Returns
    -------
    ca.SX.sym
        Right hand side of state equation

    Unlike KF family, parameters are passed separately.
    """
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]

    u = u[0]

    M = p[0]
    D = p[1]

    Ki = 2.0
    Rp = 0.05
    Tg = 0.2

    return ca.vertcat(
        x2,
        x3,
        -Ki/(M*Tg)*x1-(D/(M*Tg)+1/(Rp*M*Tg))*x2-(D/M+1/Tg)*x3-1/(M*Tg)*u,
    )


def Hc(xk, p):
    """
    Function that returns measurement variable in terms of x and p.

    Parameters
    ----------
    x : ca.SX.sym array
        State vector

    p : ca.SX.sym array
        Unknown parameter vector

    Returns
    -------
    ca.SX.sym
        Measurement variable

    Measurement in this case is $\Delta\omega$ only which is second state variable.
    """
    return xk[1]


Ts = 0.02
N = 10

nX = 3
nU = 1
nY = 1
nP = 2

# Arrival Cost
In1, Out1, InName1, OutName1 = Estimator.arrivalCost(
    nX=nX,
    nU=nU,
    nY=nY,
    nP=nP,
    Fc=Fc,
    Hc=Hc,
    Ts=Ts,
    Integrator='rk4')

Arrival_Func = ca.Function('Arrival_Func', In1, Out1, InName1, OutName1)

Opts = {
    'qpsol_options': {'print_header': False,
                      'print_iter': False,
                      'print_info': False,
                      # 'constr_viol_tol': 1e-5,
                      # 'dual_inf_tol': 1e-5,
                      }
}

In2, Out2, InName2, OutName2 = Estimator.simpleMHE(
    nX=nX,
    nU=nU,
    nY=nY,
    nP=nP,
    N=N,
    Fc=Fc,
    Hc=Hc,
    Ts=Ts,
    pLow=[0.2, 0],
    pUpp=[10, 10],
    arrival=True,
    GGN=True,
    Options=Opts,
)

MHE_func = ca.Function('MHE_func', In2, Out2, InName2, OutName2)

# Combine Arrival cost and MHE in one function
xL = ca.MX.sym('xL', nX)
pL = ca.MX.sym('pL', nP)
xLb = ca.MX.sym('xLb', nX)
pLb = ca.MX.sym('pLb', nP)

PL = ca.MX.sym('PL', (nX+nP), (nX+nP))

uL = ca.MX.sym('uL', nU)
yL = ca.MX.sym('yL', nY)

VL = ca.MX.sym('VL', nY, nY)
WL = ca.MX.sym('WL', nX, nX)
Wp = ca.MX.sym('Wp', nP, nP)

zGuess = ca.MX.sym('zGuess', nX*(N+1)+nP)
um = ca.MX.sym('um', nU, N)
ym = ca.MX.sym('ym', nY, N+1)

zOp = ca.MX.sym('zOp', nX*(N+1)+nP)

Out1Eval = Arrival_Func(xL, pL, xLb, pLb, PL, uL, yL, VL, WL, Wp)

Out2Eval = MHE_func(
    zGuess, um, ym, Out1Eval[0], Out1Eval[1], Out1Eval[2], VL, WL, zOp)

ArrivalMHE = ca.Function(
    'ArrivalMHE',
    [xL, pL, xLb, pLb, PL, uL, yL, VL, WL, Wp, zGuess, um, ym, zOp],
    [Out1Eval[0], Out1Eval[1], Out1Eval[2], Out2Eval[0],
        Out2Eval[1], Out2Eval[2], Out2Eval[3], Out2Eval[4]],
    ['xL', 'pL', 'xLb', 'pLb', 'PL', 'uL', 'yL', 'VL',
        'WL', 'Wp', 'zGuess', 'um', 'ym', 'zOp'],
    ['xLbn', 'pLbn', 'PLn', 'zOut', 'p_hat', 'x_hat', 'Cost', 'xLout']
)

BasicUtils.Gen_Code(
    func=ArrivalMHE,
    filename='MHEArrival',
    mex=False,
    printhelp=True
)
