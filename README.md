This project consists of plotting the B-field around a dipole (the earth), and modeling the trajectory of a charged particle approaching this dipole.
I firstly implemented both these tasks in python. Due to the inefficiency of the python-program i also implemented the trajectory-model in c++. 

To plot the B-field, run the 'field_plot.py'-file using Python (3.8)
To plot a particle-trajectory using the python-implementation, one needs to specify the appropriate particle-variables in 'path_plot.py'. The earth is described by the variables in the earth.py. You should then execute the 'path_plot'-python-script (3.8).

In order to plot a particle-trajectory using my c++-implementation, one needs to specify the appropriate variables in the 'in.txt'-file. Do not change the descriptive words, only their values. Furthermore, one should execute the 'run_script.sh'-shellscript. This will initiate an interactive python-session, and display a path-plot.
