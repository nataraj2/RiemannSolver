import matplotlib.pyplot as plt
import numpy as np
import glob

import glob

def count_files_starting_with(directory, prefix):
    file_pattern = directory + '/' + prefix + '*'
    matching_files = glob.glob(file_pattern)
    return len(matching_files)

# Specify the directory and prefix
directory = '../'
prefix = 'Solution'

# Count the files starting with the specified prefix
count = count_files_starting_with(directory, prefix)

print("Number of files starting with '{}': {}".format(prefix, count))

fig, axs = plt.subplots(2, 2)
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.3, hspace=0.3)
for itern in np.arange(0,count,1):
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
	
	pressure = (solution[:,3] - 0.5*solution[:,2]**2/solution[:,1])*0.4;
	uvel = solution[:,2]/solution[:,1];

	
	axs[0, 0].cla()
	axs[0, 0].plot(solution[:,0], solution[:,1])
	axs[0, 0].axis("tight")
	axs[0, 0].set_xlabel(r'$x$')
	axs[0, 0].set_ylabel(r'$Density$')	
	plt.pause(0.2)

	axs[0, 1].cla()
	axs[0, 1].plot(solution[:,0], uvel)
	axs[0, 1].axis("tight")
	axs[0, 1].set_xlabel(r'$x$')
	axs[0, 1].set_ylabel(r'$Velocity$')
	plt.pause(0.2)

	axs[1, 0].cla()
	axs[1, 0].plot(solution[:,0], pressure)
	axs[1, 0].axis("tight")
	axs[1, 0].set_xlabel(r'$x$')
	axs[1, 0].set_ylabel(r'$Pressure$')
	plt.pause(0.2)

	axs[1, 1].cla()
	axs[1, 1].plot(solution[:,0], 1.0/0.4*pressure/solution[:,1])
	axs[1, 1].axis("tight")
	axs[1, 1].set_xlabel(r'$x$')
	axs[1, 1].set_ylabel(r'$C_vT$')
	plt.pause(0.2)

plt.show()	


