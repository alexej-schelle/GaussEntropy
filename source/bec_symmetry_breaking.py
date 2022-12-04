############################################################################################################################################################################
#
#
# -*- coding: utf-8 -*-
#
#
############################################################################################################################################################################
#                                                                                                                                              
#																										  
#
#   Markov sampling method (generalised Metropolis-Hastings sampling algorithm) for calculating the
#   integrated (random) condensate wave field (order parameter) for a given total particle number, 
#   temperature and trap frequency of a Bose-Einstein condensate.
#
#   This software version corresponds to publication of Fluctuation and Noise Letters, Vol. 16, No. 01, 1750009 (2017).
#  
#   The Python routine bec_symmytry_breaking.py calculates :
#
#   1. Condensate field modes   	    	            	    	            	    	        
#   2. Condensate wave fields at equilibrium (symmetry aspects)    
#   3. Condensate wave field propagations (from low energy to high energy)    
#   4. Partial phase distributions of real valued field modes    
#   5. Partial phase distributions of imaginary valued field modes    
#   6. Total phases of the condensate wave field
#   7. Chemical potentials of the condensate 
#
# 
#   For proper installation, please replace the path /your-installation-path/ with your installation path in the lines 235 - 316 of this code.
#
#
#  * :
# 
#   License Copyright:  Dr. A. Schelle, Sudetenstr. 76, 87600 Kaufbeuren 
#   License Type :      MIT license (2017)
#   License Contact:    E-Mail : alexej.schelle@gmail.com
# 
#   ** : 
#
#   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files 
#   (the "Software" bec_symmetry_breaking.py), to deal in the Software without restriction, including without limitation the rights to use, 
#   copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
#   furnished to do so, subject to the following conditions:
# 
#   The above copyright notice (*) and this permission notice (**) shall be included in all copies or substantial portions of the Software.
# 
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
#   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#   Free support is provided per email at support@krealix.de.   
#   For more explicit consulting and discussions, please schedule a meeting at https://calendly.com/alexej-schelle/.  
#
############################################################################################################################################################################


import os
import sys
import math
import random
import numpy
import numpy as np
import pylab
import matplotlib.pyplot as plt
import operator


maxmode = 50 # Typical mode size for analysis : 50 - 500 modes.
ptn = 2000 # Typical particle number : 10^3 - 10^5
sample = 10000 # Typical sample size : 10^5 - 10^8

omx = 2.0*math.pi*25.0 # Trap frequency in x direction
omy = 2.0*math.pi*75.0 # Trap frequency in y direction 
omz = 2.0*math.pi*125.0 # Trap frequency in z direction

start_temp = 5.0 # in units of nK

rhbkb = 7.63822291E-3 # in units of nK
norm = 0.0
int_norm = 0
drop = 0
prob = 0.0
phase_0 = 0.0
z_start = 5.0
mu_start = 0.0
mu_k = 0.0

print ('Trap depth [nK]: '), maxmode*omz*rhbkb # around few muK
print ('Critical temperature [nK]: '), rhbkb*pow(ptn,1.0/3.0)*pow(omx*omy*omz,1.0/3.0)/pow(1.202,1.0/3.0) # around 1.15 mK

en_x = ['']*maxmode
en_y = ['']*maxmode
en_z = ['']*maxmode

pols_x = ['']*maxmode
pols_y = ['']*maxmode
pols_z = ['']*maxmode

pols = ['']*maxmode

x = ['']*ptn
y = ['']*ptn

mu_x = [] # collect the real part of the integrated wave field
mu_y = [] # collect the imaginary part of the integrated wave field

mu_x_collect = [] # collect the real part of the integrated wave field
mu_y_collect = [] # collect the imaginary part of the integrated wave field

mu_x_prop_collect = [] # collect the real part of the integrated wave field
mu_y_prop_collect = [] # collect the imaginary part of the integrated wave field

phase_collect = [] # collect the phase of the integrated wave field

mu_chemical_x = [] # collect the chemical potential of the condensate field
mu_chemical_y = [] # collect the chemical potential of the condensate field

for l in range(1, sample):
          
    drop = 1
    z = random.uniform(1, ptn) # Random sample of non-condensate particle number
    temp = start_temp
    
    mu = 0.0
    norm = 0.0
    
    print('Sample step Nr. ') + str(l)
                  
    for k in range(1, maxmode):
        
        en_x[k] = rhbkb*k*omx/temp # Energy in x direction
        en_y[k] = rhbkb*k*omy/temp # Energy in y direction
        en_z[k] = rhbkb*k*omz/temp # Energy in z direction
     
        pols_x[k] = 1.0/(1.0-math.exp(-en_x[k])) # Complex poles in x direction
        pols_y[k] = 1.0/(1.0-math.exp(-en_y[k])) # Complex poles in y direction
        pols_z[k] = 1.0/(1.0-math.exp(-en_z[k])) # Complex poles in z direction
                    
    for k in range(1, maxmode):
    
        pols[k-1] = pols_x[maxmode-k]*pols_y[maxmode-k]*pols_z[maxmode-k]-1.0 # General poles
    
    pols[maxmode-1] = -(ptn-z_start)
    
    prob = complex(0.0,0.0)
    p = ['']*maxmode
    phase_0 = 0.0
   
    x = numpy.roots(pols) # Complex roots of the number conserving equation

    for k in range(0,len(x)):

        p[k] = random.uniform(0.0,1.0)   # Random amplitudes - set p[k] = delta(k-k_random) for single BEC spectrum
        norm = norm + p[k]*p[k]*(x[k].real**2+x[k].imag**2) # Total norm
       
    norm = math.sqrt(norm)
    
    for k in range(0,len(x)): # Calculate phase of the condensate wave field
        	
        p[k] = p[k]/norm # Random amplitudes
    	
        if (operator.gt(x[k].real**2 + x[k].imag**2,0.0)):
                
            if (operator.gt(x[k].real,0.0)):
            
                phase_0 = math.atan(x[k].imag/x[k].real)
        
            if (operator.iand(operator.lt(x[k].real,0.0),operator.ge(x[k].imag,0.0))):
            
                phase_0 = math.atan(x[k].imag/x[k].real) + math.pi
    
            if (operator.iand(operator.lt(x[k].real,0.0),operator.lt(x[k].imag,0.0))):
	            
                phase_0 = math.atan(x[k].imag/x[k].real) - math.pi
    
            if (operator.iand(operator.eq(x[k].real,0.0),operator.gt(x[k].imag,0.0))):
            
                phase_0 = 0.5*math.pi

            if (operator.iand(operator.eq(x[k].real,0.0),operator.lt(x[k].imag,0.0))):
            
                phase_0 = -0.5*math.pi
        	       
            mu += x[k]*p[k] # Random amplitudes times phases
            mu_k = x[k]*p[k] # Condensate field modes

            mu_x_prop_collect.append(1.0*mu.real) # Collect condensate field propogation in ascending mode direction
            mu_y_prop_collect.append(1.0*mu.imag) # Collect condensate field propagation in ascending mode direction
    	
            mu_x_collect.append(1.0*mu_k.real) # Collect condensate field modes in x direction
            mu_y_collect.append(1.0*mu_k.imag) # Collect condensate field modes in y direction
    	
            mu_chemical_x.append(x[k].real) # Collect condensate chemical potential - real parts
            mu_chemical_y.append(x[k].imag) # Collect condensate chemical potential - imaginary parts
        
            if operator.ne(phase_0,0.0):
    	
                prob += complex(0.5*p[k]*p[k]*math.log(math.fabs(x[k].real**2 + x[k].imag**2)), p[k]*p[k]*phase_0) # Calculate transition probability
    	
            prob = math.sqrt(prob.real**2 + prob.imag**2) 
    	
    if (operator.gt(min((math.exp(prob))/(math.exp(mu_start)),1.00),random.uniform(0.00,1.00))): # Condition for transition to another state at equilibrium
 
        mu_start = prob
        z_start = z
        drop = 0
	    
    if (operator.ne(drop,1)):	
		
        mu_x.append(1.0*mu.real) # Collect condensate field modes in x direction
        mu_y.append(1.0*mu.imag) # Collect condensate field modes in y direction
            							    							    							    							
        if (operator.gt(mu.real,0.0)): # First quarter plane
            
            phase_0 = math.atan(mu.imag/mu.real) + math.pi
            phase_collect.append(phase_0/math.pi)
        			
        if (operator.iand(operator.lt(mu.real,0.0),operator.ge(mu.imag,0.0))): # Second quarter plane
            
    	    phase_0 = math.atan(mu.imag/mu.real) + math.pi + math.pi
    	    phase_collect.append(phase_0/math.pi)
    
        if (operator.iand(operator.lt(mu.real,0.0),operator.lt(mu.imag,0.0))): # Third quarter plane
	            
            phase_0 = math.atan(mu.imag/mu.real) - math.pi + math.pi
            phase_collect.append(phase_0/math.pi)
    
        if (operator.iand(operator.eq(mu.real,0.0),operator.gt(mu.imag,0.0))): # Fourth quarter plane
            
            phase_0 = 0.5*math.pi + math.pi
            phase_collect.append(phase_0/math.pi)

        if (operator.iand(operator.eq(mu.real,0.0),operator.lt(mu.imag,0.0))): # Fifth quartr plane
            
            phase_0 = -0.5*math.pi + math.pi
            phase_collect.append(phase_0/math.pi)
    	    	        

# Plot condensate field modes   	    	            	    	            	    	        

#plt.figure(51)
plt.hist2d(mu_x_collect, mu_y_collect, bins = 500, normed = True)
plt.tick_params(axis='both', which='major', labelsize = 16)
plt.xlabel('$Re(\Psi_0)$', fontsize = 18)
plt.ylabel('$Im(\Psi_0)$', fontsize = 18)
cbar = plt.colorbar()
cbar.ax.set_ylabel('$\pi_e[Re(\Psi_0),Im(\Psi_0))]$')
plt.savefig('/Users/AS_Scientific_Analytics/Desktop/Journal_of_Physics/bec_symmetry_breaking_github/fig_field_modes.png')


# Plot condensate wave field at equilibrium    

# plt.figure(1)
# plt.hist2d(mu_x, mu_y, bins = 300, normed = True)
# plt.tick_params(axis='both', which='major', labelsize = 16)
# plt.xlabel('$Re(\Psi_0)$', fontsize = 18)
# plt.ylabel('$Im(\Psi_0)$', fontsize = 18)
# cbar = plt.colorbar()
# cbar.ax.set_ylabel('$\pi_e[Re(\Psi_0),Im(\Psi_0))]$')
# plt.xlim([-1.00, 1.00])
# plt.ylim([-1.00, 1.00])
# plt.savefig('/Users/AS_Scientific_Analytics/Desktop/Journal_of_Physics/bec_symmetry_breaking/fig_field_thermal.png')

# Plot non-condensate wave field propagation    

#plt.figure(53)
#plt.hist2d(mu_x_prop_collect, mu_y_prop_collect, bins = 500, normed = True)
#plt.tick_params(axis='both', which='major', labelsize = 16)
#plt.xlabel('$Re(\Psi_0)$', fontsize = 18)
#plt.ylabel('$Im(\Psi_0)$', fontsize = 18)
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('$\pi_e[Re(\Psi_0),Im(\Psi_0))]$')
#plt.savefig('/Users/AS_Scientific_Analytics/Desktop/bec_symmetry_breaking_github/fig_field_propagation.png')


# Plot phase distribution    

#plt.figure(54)
#plt.hist(mu_x, bins = 300, normed = True)
#plt.tick_params(axis='both', which='major', labelsize = 16)
#plt.xlabel('$Re(\Psi_0)$', fontsize = 18)
#plt.ylabel('$\pi_e[Re(\Psi_0)]$', fontsize = 18)
#plt.xlim([-1.00, 1.00])
#plt.savefig('/Users/AS_Scientific_Analytics/Desktop/bec_symmetry_breaking_github/fig_distribution_real.png')


# Plot phase distribution    

#plt.figure(55)
#plt.hist(mu_y, bins = 300, normed = True)
#plt.tick_params(axis='both', which='major', labelsize = 16)
#plt.xlabel('$Im(\Psi_0)$', fontsize = 18)
#plt.ylabel('$\pi_e[Im(\Psi_0)]$', fontsize = 18)
#plt.xlim([-1.00, 1.00])
#plt.savefig('/Users/AS_Scientific_Analytics/Desktop/bec_symmetry_breaking_github/fig_distribution_imag.png')


# Plot phase of wave field

#plt.figure(56)
#plt.hist2d(phase_collect, phase_collect, bins = 500, normed = True)
#plt.tick_params(axis='both', which='major', labelsize = 16)
#plt.xlabel('$\phi_0$', fontsize = 18)
#plt.ylabel('$\pi_e(\phi_0)$', fontsize = 18)
#plt.xlim([0.00, 2.00])
#plt.ylim([0.00, 2.00])
#cbar.ax.set_ylabel('$\pi_p[\phi_0] [\pi]$')
#plt.savefig('/Users/AS_Scientific_Analytics/Desktop/bec_symmetry_breaking_github/fig_distribution_phase.png')


# Chemical potential of the condensate

#plt.figure(57)
#plt.hist2d(mu_chemical_x, mu_chemical_y, bins = 200, normed = True)
#plt.tick_params(axis='both', which='major', labelsize = 16)
#plt.xlabel('$Re(\mu_0)$', fontsize = 18)
#plt.ylabel('$Im(\mu_0)$', fontsize = 18)
#cbar.ax.set_ylabel('$\pi_e[Re(\mu_0),Im(\mu_0)]$')
#plt.colorbar()
#plt.savefig('/Users/AS_Scientific_Analytics/Desktop/bec_symmetry_breaking_github/fig_distribution_chemical.png')
