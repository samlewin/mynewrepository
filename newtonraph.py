# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 20:27:03 2018

@author: sam
"""
import math
TOL=10**-6
def newton_raph(func,deriv,x0,tol,maxiterations=100,toPrint=False):
    """Performs the Newton-Raphson method for a pre-defined function
    func with pre-defined derivative deriv. x0 is the starting value 
    and tol is the accuracy. Will run up to maxiterations number of 
    iterations unless the derivative becomes unbounded."""
    small=10**-4
    nextiter=x0
    count=0
    if toPrint:
        print('Initial value: '+str(nextiter))
    for i in range(maxiterations):
        if abs(deriv(nextiter))<small:
            raise ValueError('No convergence')
        else:
            new=max(TOL,nextiter-func(nextiter)/deriv(nextiter))
            count+=1
            if toPrint:
                print('Next iteration: ',new)
            if 2*abs(new-nextiter)/(new+nextiter)<tol:
                return (new,count)
            nextiter=new
    return 'Did not converge in '+str(maxiterations)+' iterations'
