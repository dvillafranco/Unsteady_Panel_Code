%% Notes on folder "Airfoil Case 1"

Main code is called run_comp_af_old.m
	a. casc_justaf
		- sets up discretization for the airfoil and calculates panel lengths etc 
		
	b. find_gams_justaf_only.m
		- sets up matrices to solve the steady problem first and obtain the circulation
		around the airfoil for the steady case 
		
	c. psionsufrace_original.m
		- gets psi on surface in order to get velocity at the trailing edge. This is what
		we use to define the first position of the wake point. From the obtained velocity
		on the panels, the pressure coefficient and the lift coefficient are also 
		calculated.
		
	d. find_gams_justaf_only_withvortex_shed.m
		- sets up the number of time points to loop through for shedding wake.
		The imposed vortex's strength, and position are defined in this code. This code
		sets up the matrices necessary to solve the "unsteady" problem
		
		i. psionsurface
			- recalculates the velocity at the trailing edge of the airfoil.
			Uses this velocity to 'convect' the wake downstream. The imposed
			vortex is convected with the free-stream velocity. Since the velocities
			have been recalculated each loop. The coefficient of pressure and also the 
			coefficient of lift for each 'time-step' is found in this code. 
			
			x_fixed(nx+1) is the form of the stored wake points. If we are recalculating 
			the trailing edge velocity at each time step then only the 'newest' shed 
			vortex should be convected off the airfoil with that velocity. How should the 
			remaining vortices moved in the fixed wake case? With free-stream?
			
		
General Notes:
	- This code solves the unsteady problem for an airfoil with a vortex imposed 
	upstream. That vortex is then convected downstream with the free-stream velocity. 
	In the case of the wake convection, the trailing edge velocity is recalculated each 
	time and the wake is convected with that velocity. 
	
	
	
	
