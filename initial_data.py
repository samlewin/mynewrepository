#Import the initial data
test0 = input("Enter name of test file:")
test1 = input("Enter x0:")
test2 = input("Enter mesh divisions:")
file = open(test0,'r')
constants = []
for line in file:
    constants.append(line.strip('\n'))
file.close()
uL, uR, pL, pR, rhoL, rhoR, gamma, T =(float(constants[1][3:]), float(constants[4][3:]),
                             float(constants[2][3:]), float(constants[5][3:]),
                             float(constants[0][5:]), float(constants[3][5:]),
                             float(constants[-2][6:]), float(constants[-1][2:]))
def initdata(x):
    if x < x0:
        return (uL, rhoL, pL)
    else:
        return (uR, rhoR, pR)
mesh, x0 = int(test2), float(test1)
dx = 1/mesh
xvals = [dx*x for x in range(1,mesh + 1)] #these are the x_(1/2),...,x_(100+1/2)
inituvals = [initdata(xvals[0])[0]] + [initdata(xi)[0] for xi in xvals] + [initdata(xvals[-1])[0]]
initrhovals = [initdata(xvals[0])[1]] + [initdata(xi)[1] for xi in xvals] + [initdata(xvals[-1])[1]]
initpvals = [initdata(xvals[0])[2]] + [initdata(xi)[2] for xi in xvals] + [initdata(xvals[-1])[2]]
