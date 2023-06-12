import matplotlib.pyplot as plt
import numpy as np
import glob


for itern in np.arange(0,80,1):
	if(itern<=9):
		solution_file = "../Solution_000" + str(itern) + ".txt";
	elif(itern<=99):
		solution_file = "../Solution_00" + str(itern) + ".txt";
	elif(itern<=999):
		solution_file = "../Solution_0" + str(itern) + ".txt";
	elif(itern<=9999):
		solution_file = "../Solution_" + str(itern) + ".txt";

	solution = np.loadtxt(solution_file);
	print(solution_file)
	plt.figure(1)
	
	pressure = (solution[:,3] - 0.5*solution[:,2]**2/solution[:,1])*0.4;
	uvel = solution[:,2]/solution[:,1];

	#plt.plot(solution[:,0], solution[:,1])
	plt.plot(solution[:,0], pressure)
	plt.axis("tight")
	plt.pause(0.2)
	plt.clf()

plt.show()	

plt.figure(2)
plt.plot(solution[:,0], solution[:,1])
exact = np.loadtxt("SodDensityExact.dat")
plt.plot(exact[:,0], exact[:,1],'r')
plt.savefig("Comparison_Density.png")

plt.figure(3)
plt.plot(solution[:,0], pressure)
exact = np.loadtxt("SodPressureExact.dat")
plt.plot(exact[:,0], exact[:,1],'r')
plt.savefig("Comparison_Pressure.png")

plt.figure(4)
plt.plot(solution[:,0], uvel)
exact = np.loadtxt("SodVelocityExact.dat")
plt.plot(exact[:,0], exact[:,1],'r')
plt.savefig("Comparison_Velocity.png")
plt.show()



#print(solution_file)
	
#plt.figure(1)
#plt.plot(solution[:,0], solution[:,1])
#plt.show()
#plt.pause(0.2)
#plt.clf()

