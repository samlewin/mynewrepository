from godunov_with_comparison import MLgodunov
from HLLCgodunov import HLLCgodunov
import initial_data as i
import numpy as np
import pylab as plt
import solve, time

#Time the process
result = MLgodunov(np.array(i.initrhovals),np.array(i.inituvals),np.array(i.initpvals),i.gamma,i.T,i.dx)

result1 = HLLCgodunov(np.array(i.initrhovals),np.array(i.inituvals),np.array(i.initpvals),i.gamma,i.T,i.dx)

#Pressure
fig1 = plt.figure()
plt.plot(i.xvals,result1[2],'--k',
         label='HLLC')
plt.plot(i.xvals,result[2],'-g',
         label='ML')
yvals=[solve.solve(xi-i.x0,i.T,i.uL,i.uR,i.pL,i.pR,i.rhoL,i.rhoR,i.gamma,TOL=10**-6)[2] for xi in i.xvals]
plt.plot(i.xvals,yvals,'r',label='Exact Solution')
plt.legend()
plt.xlim((0,1))
plt.xlabel('x')
plt.ylabel('p')
plt.xticks((0,0.5,1))
plt.title('Pressure Plot ' +'\n' + 'Time='+str(i.T)+', Steps='+
          str(result[3]))

fig1.savefig("pressure.pdf", bbox_inches='tight')

##Velocity
fig2 = plt.figure()
plt.plot(i.xvals,result[0],'-g',
         label='ML')
plt.plot(i.xvals,result1[0],'--k',
         label='HLLC')
yvals=[solve.solve(xi-i.x0,i.T,i.uL,i.uR,i.pL,i.pR,i.rhoL,i.rhoR,i.gamma,TOL=10**-6)[0] for xi in i.xvals]
plt.plot(i.xvals,yvals,'r',label='Exact Solution')
plt.legend()
plt.xlim((0,1))
plt.xticks((0,0.5,1))
plt.xlabel('x')
plt.ylabel('u')
plt.title('Velocity Plot ' +'\n' + 'Time='+str(i.T)+', Steps='+
          str(result[3]))

fig2.savefig("velocity.pdf", bbox_inches='tight')

##Density
fig3 = plt.figure()
plt.plot(i.xvals,result[1],'-g',
         label='ML')
plt.plot(i.xvals,result1[1],'--k',
         label='HLLC')
yvals=[solve.solve(xi-i.x0,i.T,i.uL,i.uR,i.pL,i.pR,i.rhoL,i.rhoR,i.gamma,TOL=10**-6)[1] for xi in i.xvals]
plt.plot(i.xvals,yvals,'r',label='Exact Solution')
plt.legend()
plt.xlim((0,1))
plt.xticks((0,0.5,1))
plt.xlabel('x')
greek_letterz=[chr(code) for code in range(945,970)]
plt.ylabel(greek_letterz[16])
plt.title('Density Plot ' +'\n' + 'Time='+str(i.T)+', Steps='+
          str(result[3]))

fig3.savefig("density.pdf", bbox_inches='tight')

##Errors
fig4 = plt.figure()
plt.plot(i.xvals,result[4][:-1],'g', label='HLLC')
plt.plot(i.xvals,result[5][:-1],'r', label='ML')
plt.legend()
plt.title('Mean HLLC error vs ML error')
plt.xlabel('x')
plt.ylabel('error')

fig4.savefig("errors.pdf", bbox_inches='tight')
plt.show()
