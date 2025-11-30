# SpringerDartBlasterSim
Springer Dart Blaster Simulator

All models are wrong. Some models are useful.

This project simulates a typical springer dart blaster, where a compressed spring pushes on a plunger, forcing air into a barrel and propelling a foam dart forwards.
There are a lot of variables that affect the exit velocity of the dart. This simulator attempts to account for the most common variables, while simplifyhing the complexity of others.

This model is in a prototype state. It functions and results in realistic exit velocities, but has not been verified against real world experimental setups.
In it's current state, it is recommended to only use the model to compare relative performance between different blaster configurations. Do not trust that the absolute FPS is correct.

This model simulates two cylinders connected by a volume of "dead space" where the cylinders do not travel. This would be sufficient is the event took place slowly, however a significant limiting factor of dart blasters comes from flow through a small orifice.
The model used here is the isentropic model for compressible subsonic flow through an orifice: https://web.cecs.pdx.edu/~gerry/class/ME322/notes/pdf/compressibleMdot.pdf
 
The orifice is simply the smallest cross sectional area in the flow. Oftentimes this may be at the back of the barrel, in whcih case the dead volume of the barrel is 0. 

<img width="1802" height="611" alt="image" src="https://github.com/user-attachments/assets/7aef20bc-0439-41b9-a289-321d5e659c80" />

The model currently assumes that all air seals are perfect, and mass is conserved within the system.


The simulation results are shown using matplotlib, and the exit fps is printed to the console. The remaining plunger travel when the dart exits the barrel is also shown.
 
<img width="800" height="921" alt="image" src="https://github.com/user-attachments/assets/b468b317-f38e-4938-9238-24627486b7cb" />
