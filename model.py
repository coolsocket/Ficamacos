#!/usr/bin/env python
# coding: utf-8

# In[ ]:



#Models for option pricing


"""
GBM(Geometric Brownian Motion)
=========================
"""

import numpy as np

__all__ = ['GBM', 'GBMParam','VG', 'VGParam']


class GBMParam(object):

    """Parameter storage.
    """

    def __init__(self, sigma=.2):
        """Initialize class.
        Parameters
        ----------
        sigma : positive float
        """
        self.sigma = sigma


class GBM(object):

    """Geometric Brownian Motion.
    Attributes
    ----------
    param
        Model parameters
    Methods
    -------
    charfun
        Characteristic function
    cos_restriction
        Restrictions used in COS function
    """

    def __init__(self, param, riskfree, maturity):
        """Initialize the class.
        Parameters
        ----------
        param : GBMParam instance
            Model parameters
        riskfree : float
            Risk-free rate, annualized
        maturity : float
            Fraction of a year
        """
        self.param = param
        self.riskfree = riskfree
        self.maturity = maturity

    def charfun(self, arg):
        """Characteristic function.
        Parameters
        ----------
        arg : array_like
            Grid to evaluate the function
        Returns
        -------
        array_like
            Values of characteristic function
        """
        return np.exp(arg * self.riskfree * self.maturity * 1j
                      - arg**2 * self.param.sigma**2 * self.maturity / 2)

    def cos_restriction(self):
        """Restrictions used in COS function.
        Returns
        -------
        a : float
        b : float
        """
        # Truncation rate
        L = 10
        c1 = self.riskfree * self.maturity
        c2 = self.param.sigma**2 * self.maturity

        a = c1 - L * c2**.5
        b = c1 + L * c2**.5

        return a, b
    
    

"""
VG(Variance Gamma)
=========================
"""


class VGParam(object):

    """Parameter storage.
    """

    def __init__(self,sigma=0.01,theta=0.005,v=0.1):
        """Initialize class.
        Parameters
        ----------
        a : positive float
        b : positive float
        """
        self.sigma = sigma
        self.v = v
        self.theta = theta


class VG(object):

    """Variance Gamma Process.
    Attributes
    ----------
    param
        Model parameters
    Methods
    -------
    charfun
        Characteristic function
    cos_restriction
        Restrictions used in COS function
    """

    def __init__(self, param, riskfree, maturity):
        """Initialize the class.
        Parameters
        ----------
        param : VGParam instance
            Model parameters
        riskfree : float
            Risk-free rate, annualized
        maturity : float
            Fraction of a year
        """
        self.param = param
        self.riskfree = riskfree
        self.maturity = maturity

    def charfun(self, arg):
        """Characteristic function.
        Parameters
        ----------
        arg : array_like
            Grid to evaluate the function
        Returns
        -------
        array_like
            Values of characteristic function
        """
        
        sigma=self.param.sigma
        theta=self.param.theta
        v = self.param.v
        T=self.maturity
        r=self.riskfree
        
        return (1-1j*arg*theta*v+1/2*sigma**2*v*arg**2)**(-T/v)*np.exp(1j*arg*T*(r- 1/v*np.log(1-theta*v-1/2*sigma**2*v)))

    def cos_restriction(self):
        """Restrictions used in COS function.
        Returns
        -------
        a : float
        b : float
        """
        # parameters
        sigma=self.param.sigma
        theta=self.param.theta
        v = self.param.v
        T=self.maturity
        r=self.riskfree
        
        # Truncation rate
        c1 = r * T
        c2 = sigma**2 * T + v * theta **2 *T 
        c4 = 3*(1+2*v/T-v*sigma**4*(c2**(-2))*T)

        return c1, c2, c4

