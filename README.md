# Improving-Lagrangian-mesh-grid
In this project we want to improve a Lagrangian mesh grid. It can be useful when we want to simulate the motion of a fluid for example.  
At the beginning we have Cartesian grid with markers(yellow circles) at the center of cells (in blue). Then the markers of the grid are going to move according to the following velocity field :  
v(x,y) = (0 , 0) if sqrt(x²+y²) >= 1  
v(x,y) = ((-1+sqrt(x²+y²))*y , -(-1+sqrt(x²+y²))*x) if sqrt(x²+y²) < 1  
The current position of each marker is indicated by a a black +.  
The goal is to have at each time step at least one marker inside each cell. If it is not the case we have to add new ones in a Cartesian grid more accurate than the one we have at the beginning. Their initial positions (at time t=0) are indicated by green circles and their current position by red +.  
To reach this goal, we use a method using the closest markers of the empty cell to know where we can add a marker at t=0 in order to have a marker in each empty cell.  
You can run the program by executing the exec.py file. At the beginning of the file you will see a 2-boolean list named execute. The first boolean is for executing the programm on a given situation : we move the initial markers nb_iter times and then we use the algorithm. You can set detail to True if you want the algorithm to run step by step. If it is the case, each time you close the window an other one will be opened and you will see a pink cross which indicate the exact position we compute for the initial position of the marker we add.  
The second boolean is for real time simulation. Set it to True and you will see the computer adding new markers when it's necessary.  
The algo.py file contains all the stuff we need and use in exec.py.  
This project uses numpy and matplotlib.
