#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
CM (Carr and Madan) method
==========
The method comes from Carr and Madan

References
----------
.. [1] Carr,  P.,  &  Madan,  D.  (1999).  Option  valuation  using  the  fastFourier transform. 
    The Journal of Computational Finance, 2(4), 61–73.https://doi.org/10.21314/jcf.1999.043
"""

import numpy as np
from numpy.fft import fft

__all__=['FFT', 'Filon', 'COS']

def FFT(model, S, K, alpha=3, N=2**15):
    """CM method using Fast Fourier Transform in Discrete Fourier inverse.
    Inverts characteristic function to obtain the density.
    Parameters
    ----------
    model : your characteristic function model
            It contains characteristic function dependent only on u
    S : float, 
        Initial stock price
    K : float, 
        The strike price 
    N : int, optional
        Number of discrete points for evaluation , it should be the power of 2
    eps : float,
        the lambda value in the paper
    alpha : positive float
        change the value of it will influence the accuracy of the result, shouldn't be too large
    Returns
    --------
    price of the option.
    
    Notes
    -----
    `charfun` method (risk-neutral conditional chracteristic function)
    of `model` instance should depend on
    one argument only (array_like) and should return
    array_like of the same dimension.
    """
        
    if not hasattr(model, 'charfun'):
        raise Exception('Characteristic function is not available!')
        
    # risk free interest
    r = model.riskfree
    
    # maturity date
    T = model.maturity
    
    # log of strike price
    k = np.log(K)
    
    # log of initial price
    s = np.log(S)
    
    # value of eta
    eta = 0.25
    
    # value of eps
    eps = 2 * np.pi / (N * eta)
   
    #value of b
    b = 0.5 * N * eps
    
    # value of u
    u = np.arange(1, N + 1, 1)
    
    # value of v
    vo = eta * (u - 1)
    
    # value of v after modification
    v = vo - (alpha + 1) * 1j
    
    # value for modified charfun
    modfun = np.exp(-r * T) * (model.charfun(v)*np.exp(1j*v*s)/(alpha ** 2 + alpha - vo ** 2 + 1j * (2 * alpha + 1) * vo))
    
    # kronecker delta 
    delt = np.zeros(N, dtype=np.float)
    delt[0] = 1
    
    # index
    j = np.arange(1, N + 1, 1)
    
    # value for simpsons weight
    SimpsonW = (3 + (-1) ** j - delt) / 3
    
    # value of each term
    FFTFunc = np.exp(1j * b * vo) * modfun * eta * SimpsonW

    # doing fft
    payoff = (fft(FFTFunc)).real
    
    # vector for values
    values = np.exp(-alpha * k) / np.pi * payoff
    
    # index for the value we want
    pos = int((k+b) / eps)
    
    # final result
    return values[pos]



"""
Filon's method
==========
The method comes from Filon. It is a good way to calculate the integral of oscilating functions.



References
----------
.. [1] Kuznetsov, A., Kyprianou, A. E., \& Rivero, V. (2012). 
    The Theory of Scale Functions for Spectrally Negative Lévy Processes. 
    
    [2] Filon, L. N. G. (1930). III.—On a Quadrature Formula for Trigonometric Integrals. 
    Proceedings of the Royal Society of Edinburgh, 49(1), 38–47. 
    https://doi.org/10.1017/s0370164600026262
"""




def Filon(model, S, K, N=2**18, b=100, alpha=3):
    """Filon's method.
    Parameters
    ----------
    model : instance of specific model class
        The method depends on availability of characteristic function
    S : array_like
        The initial stock price
    K : array_like
        The strike price
    N : int
        Number of points on the grid. The more the better, but slower.
    b : float
        upper limit for the integral.
    alpha : positive float
        change the value of it will influence the accuracy of the result, shouldn't be too large    
    Returns
    -------
    array_like
         asset price
        
    Notes
    -----
    `charfun` method (risk-neutral conditional chracteristic function)
    of `model` instance should depend on
    one argument only (array_like) and should return
    array_like of the same dimension.
    """
    
    if not hasattr(model, 'charfun'):
        raise Exception('Characteristic function is not available!')
        
    # risk free
    r = model.riskfree
    
    # maturity date
    T = model.maturity
    
    # log of strike price
    k = np.log(K)
    
    # factor for the whole function
    factor = np.exp(-alpha*k)/np.pi*np.exp(-r*T)
    
    # final result
    inte = (filon_cos(0,b,N,k ,S ,alpha,model)+filon_sin(0,b,N,k, S, alpha,model))
    return factor*inte  
    
# function used to calculate functions contain cos terms.
def filon_cos(a,b,N,x,S,alpha,model):
    h = (b-a)/(2*N)
    return (h*A(h*x)*(f(b, S, alpha,model)*np.sin(b*x)-f(a, S, alpha,model)*np.cos(a*x))            +h*B(h*x)*(c_cos(a,b,2,x,N, S, alpha,model)-1/2*(f(b, S, alpha,model)*np.cos(b*x)-f(a, S, alpha,model)*np.cos(a*x)))            +h*C(h*x)*c_cos(a,b,1,x,N, S, alpha,model))



# function used to calculate functions contain sin terms.
def filon_sin(a,b,N,x,S,alpha,model):
    h = (b-a)/(2*N)
    return (h*A(h*x)*(g(a, S, alpha,model)*np.cos(a*x)-g(b, S, alpha,model)*np.cos(b*x))            +h*B(h*x)*((c_sin(a,b,2,x,N, S, alpha,model)))            +h*C(h*x)*c_sin(a,b,1,x,N, S, alpha,model))


# values for A, B and C.
def A(x):
    return 1/x+np.sin(2*x)/(2*x**2)-2*np.sin(x)**2/(x**3)

def B(x):
    return 2*(((1+np.cos(x)**2)/x**2)-np.sin(2*x)/(x**3))

def C(x):
    return 4*(np.sin(x)/(x**3)-np.cos(x)/(x**2))

# cos terms 
def c_cos(a,b,j,x,N, S, alpha,model):
    uvec = np.linspace(a,b,num=2*N)
    uodd = uvec[1::2]
    ueven = uvec[0::2]
    if j==1:
        u = uodd
    else:
        u = ueven
    return np.dot(f(u, S, alpha,model)*np.ones_like(u),np.cos(x*u))

# sin terms
def c_sin(a,b,j,x,N, S, alpha,model):
    uvec = np.linspace(a,b,num=2*N)
    uodd = uvec[1::2]
    ueven = uvec[0::2]
    if j==1:
        u = uodd
    else:
        u = ueven
    value = np.dot(g(u, S, alpha,model)*np.ones_like(u),np.sin(x*u))
    return value



# function used to calculate cos integral
def f(x, S, alpha,model):
    return (model.charfun(x-1j*(alpha+1))*np.exp(np.log(S)*1j*(x-1j*(alpha+1)))/(alpha**2+alpha-x*x+1j*x*(2*alpha+1))).real

# function used to calculate sin integral
def g(x, S, alpha,model):
    return (model.charfun(x-1j*(alpha+1))*np.exp(np.log(S)*1j*(x-1j*(alpha+1)))/(alpha**2+alpha-x*x+1j*x*(2*alpha+1))).imag



"""
COS method
==========
The method comes from Fang and Oosterlee [1], and the formula I am using is from [2].

References
----------
.. [1] Fang, F., & Oosterlee, C. W. (2009).
    A Novel Pricing Method for European Options
    Based on Fourier-Cosine Series Expansions.
    *SIAM Journal on Scientific Computing*, 31(2), 826. doi:10.1137/080718061
    <http://ta.twi.tudelft.nl/mf/users/oosterle/oosterlee/COS.pdf>
    
    [2] Ruijter, M. J., & Oosterlee, C. W. (2015). Notes on the BENCHOP imple-mentations for the COS method.
"""

import numpy as np



def COS(model,S, K, N=2**8, L=10):
    """COS method.
    Parameters
    ----------
    model : instance of specific model class
        The method depends on availability of two methods:
            - charfun
            - cos_restriction
    S : array_like
        The initial stock price
    K : array_like
        The strike price
    N : int
        Number of points on the grid. The more the better, but slower.
        
    Returns
    -------
    array_like
        Option premium normalized by asset price
        
    Notes
    -----
    `charfun` method (risk-neutral conditional chracteristic function)
    of `model` instance should depend on
    one argument only (array_like) and should return
    array_like of the same dimension.
    `cos_restriction` method of `model` instance takes `maturity`
    and `riskfree` as array arguments,
    and returns two corresponding arrays (a, b).
    """
    
    # Error alert
    if not hasattr(model, 'charfun'):
        raise Exception('Characteristic function is not available!')
    if not hasattr(model, 'cos_restriction'):
        raise Exception('COS restriction is not available!')
    
    # unit of lower limit and upper limit from model
    c1, c2, c4 =  model.cos_restriction()
    
    # lower limit and upper limit
    a = c1 - L * (c2+c4**0.5)**0.5
    b = c1 + L * (c2+c4**0.5)**0.5
    
    # riskfree parameter
    r = model.riskfree
    
    # maturity time
    T = model.maturity
    
    # log moneyness
    x = np.log(S/K)
    
    # weight vector for each add term, the first term is a half
    wvec =  np.append(.5, np.ones(N-1))
    
    # vector of j
    jvec = np.arange(N)
    
    # modified j
    modj = jvec*np.pi/(b-a)
    
    # terms of phi
    phivec = model.charfun(modj)*np.exp(1j*jvec*np.pi*(x-a)/(b-a))
    
    # parameters for Chi and Psi
    arg = (modj, 0, b, a, b)
    
    # terms of u
    uvec = 2/(b-a)*K*(Chi(*arg)-Psi(*arg))
    
    # return value
    return np.exp(-r*T)*np.dot(wvec,phivec*uvec).real
    
    
    # Chi function
def Chi(j, z1, z2, a, b):
    return 1/(1+(j**2))*(np.cos(j*(z2-a))*np.exp(z2)-np.cos(j*(z1-a))*np.exp(z1)                +j*np.sin(j*(z2-a))*np.exp(z2)-j*np.sin(j*(z1-a))*np.exp(z1))
        
    # Psi function
def Psi(j ,z1, z2, a, b):
    # first value is different to others
    
    j = j[1:]
    value = (np.sin(j*(z2-a))-np.sin(j*(z1-a)))/j
    value = np.insert(z2-z1,1,value)
    return value








