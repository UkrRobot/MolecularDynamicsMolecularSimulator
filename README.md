# MolecularDynamicsMolecularSimulator

The method of the molecular dynamics it consists of numerical solution of Newton’s equations for all atoms—at the given instant, the incremental values of coordinates and velocities of all particles are calculated using discretizated Newton’s differential equations from their values at the previous time moment. 
First of all, it is necessary to define the system model, which is the subject of modeling. For the study of qualitative properties of the system of many particles, the assumptions that dynamics is classical and molecules are chemically inert are right enough. 
Model corresponds to “simple” fluids, which consist of neutral atoms or molecules, the liquid argon in this case. It is known that argon atoms feebly attract each other at comparatively large r; it is caused mainly by the cross-polarization of atoms; the resulting force of attraction is termed the van der Waals force. The Lennard-Jones potential describes well such interaction; 
The interaction potential between two particles is calculated in the program by means of a function, which returns the value of the interaction energy between two particles depending on the distance between them, dr5squart(dx21dy2) (dx and dy are the differences of the particle coordinates).
double potential_direct(double dx,double dy) { double d2,d6; d2 = double(1.0/(dx*dx+dy*dy)); // divided by square of distance d6 = double(d2*d2*d2); // obtaining of the 6th degree of expression double potential = (potential_a*d6 - potential_b)*d6; return   potential; }
The disposition of particles in the restricted volume may be random or in the knots of some implicit grid. It is also necessary to choose the boundary conditions. There is a possibility to choose one of two types of boundary conditions in the program: elastic walls or periodic conditions. The elastic potential for boundaries is chosen in the form of Hooke’s potential.
For evaluating the particles’ motion and calculating the pressure produced by them, calculations of forces are required. This is fulfilled by the procedures similar in code to the functions evaluating the potentials given above. Using procedures instead of functions as in the case of the potential energy calculations is caused by the fact that the force is the vectorial magnitude, and it is necessary to count its.
The velocity form of Verlet’s algorithm is applied most often in scientific works. The more exact Bimann’s algorithm. The additional array of accelerations (accelerations are equal to forces because the reduced mass is equal to 1) calculated on the previous timestep is used in the procedure (three arrays: accel_1 [i], accel_2 [i], and accel_3 [i] are used in it): The procedure (void) Step is preceded by procedure Act_atoms, which is called for calculating of acceleration arrays: accel_1, accel_2, accel_3 (accel[i].vx i accel[i].vy)
As the program provides application of dimensionless quantities, it is necessary to provide relations between dimensional and dimensionless quantities or, in other words, define the unit values for transforming dimensionless quantities into real ones.

Recommended Experiments
1. Analyze the basic algorithm and its implementation. Justify connections of real physical characteristics used for building the model and the dimensionless parameters, which are applied in this program. 
2. Find out how the time step influences the accuracy of calculations. 
3. Define the basic characteristics of the studied system by analyzing the modeling results. Convert the values of the pressure and diffusivity to SI units. 
4. Study the temperature influence on the structure, that is, on the placement of atoms and RPDF (consider the temperature from 10 to 100 K). Study and explain how boundary conditions influence the structure. 
5. Construct the void for the graph of the velocity distribution of particle. 
6. Elaborate the void for definition of the velocity auto-correlation function and provide a calculation of the self-diffusion coefficient using the autocorrelation function. 
7. Develop the void for the system thermostatic control. Define the “melting” temperature for the system that contacts with the thermostat.

![alt text](https://github.com/UkrRobot/MolecularDynamicsMolecularSimulator/blob/master/MolecularDynamicsScrn.png)

This program was developed under the guidance of Professor A. Ovrutsky.
