#!/usr/bin/python

# Exact solution for the Riemann problem - Sod shock tube

from numpy import *
#from pylab import *
from matplotlib import rc, rcParams
from matplotlib.pyplot import *

x=arange(-10,10,0.1)

x = linspace(-10,10,200)
#print x
print size(x)

#print x

# Initial conditions for the Riemann problem
#-----------------------------------
#	d4	|  d1
#	u4	|  u1
#	p4	|  p1	
#-----------------------------------
d4 = 1.0 #kg/m^3
u4 = 0.0 #m/s
p4 = 100000 #N/m^2
d1 = 0.125
u1 = 0.0
p1 = 10000

#Constants
gamma = 1.4

a4 = sqrt(gamma*p4/d4)
#print a4
a1 = sqrt(gamma*p1/d1)
#print a1

p4byp1 = p4/p1

#Find p2/p1

for p2byp1 in arange(1,4,0.01):
	fac1 = (p2byp1-1)/sqrt((gamma+1)/(2.0*gamma)*(p2byp1-1)+1);
	fac2 = -2*gamma/(gamma-1);
	p4byp1check = p2byp1*( 1.0 + (gamma-1)/(2*a4)*(u4-u1-a1/gamma*fac1) )**fac2;
	#print p2byp1, p4byp1check, p4byp1
        #print p2byp1
	if(abs(p4byp1check-p4byp1) <= 0.05):
		print p2byp1,p4byp1check
  		break

p2 = p2byp1*p1	
p3 = p2
u3 = u4+2*a4/(gamma-1)*(1-(p3/p4)**((gamma-1)/(2*gamma)))
u2 = u3
print u3	   

a3 = (gamma-1)/2*(u4-u3+2*a4/(gamma-1))
print a3 

#Find the shock speed S

S=u1 + a1*sqrt(( gamma+1)/(2*gamma)*(p2byp1-1)+1 )
print S

#Determine the different regions and write the exact solution
t = 0.01
p = zeros(size(x))
u = zeros(size(x))
d = zeros(size(x))
for i in arange(0,size(x),1):
	if(x[i]/t <= (u4-a4)):
		p[i] = p4
   		u[i] = u4
		d[i] = d4
		#plot(x,p,'o')
		#print "here"
	if(x[i]/t >= (u4-a4) and x[i]/t <= (u3-a3)):
		u[i] = 2.0/(gamma+1)*(x[i]/t+(gamma-1)/2*u4+a4)
		a = u[i] - x[i]/t
		p[i] = p4*(a/a4)**(2*gamma/(gamma-1))
		d[i] = gamma*p[i]/a**2
		#plot(x,p,'o')
		#print "here too"
		dval = d[i]
	if(x[i]/t>=(u3-a3) and x[i]/t<=u3):
		p[i] = p3
		u[i] = u3
		d[i] = dval
		#plot(x,p,'o')
	if(x[i]/t>=u3 and x[i]/t<=S):
		p[i] = p2
		u[i] = u2
		d[i] = d1*(u1-S)/(u2-S)
		#plot(x,p,'o')
	if(x[i]/t>=S):
		p[i] = p1
		u[i] = u1
		d[i] = d1
		#plot(x,p,'o')

plot(x,d)	
#axis([-10,10,p1-10000,p4+10000])
#print p[199]
xlabel(r'$x$(m)')
ylabel(r'$p(N/m^2)$')
show()
savefig('PressureExactRiemann.eps')	
array = column_stack((x,p))
print array
savetxt('SodPressureExact.dat',array)
array = column_stack((x,u))
savetxt('SodVelocityExact.dat',array)
array = column_stack((x,d))
savetxt('SodDensityExact.dat',array)
#f = open("values.dat", "wt")
#f.write("x y\n");
#for i in linspace(0,size(x)-1,size(x)):
#	f.write("%g %g\n" % (x[i],p[i]))

