%%Free_Wake_Single_Airfoil_README

%% Notes on folder "Airfoil_Unsteady_FreeWake

Main code is called run_comp_af_old.m
	a. casc_justaf
		- sets up discretization for the airfoil and calculates panel lengths etc 
		
	b. find_gams_justaf_only.m
		- sets up matrices to solve the steady problem first and obtain the circulation
		around the airfoil for the steady case 
		
	c. psionsurface_original.m
		- gets psi on surface in order to get velocity at the trailing edge. This is what
		we use to define the first position of the wake point. From the obtained velocity
		on the panels, the pressure coefficient and the lift coefficient are also 
		calculated.
		
	
	d. find_gams_justaf_withvortex;
		- sets up the number of time points to loop through for shedding wake.
		The imposed vortex's strength, and position are defined in this code. This code
		sets up the matrices necessary to solve the "unsteady" problem
		
		Note that gam_dim is the vector which holds the strength of all the vortices 
		in the wake. It is ordered such that the 'newest' vortex is always at the front 
		of this vector. 
		
			i. psionsurface.m
			- recalculates the velocity at the trailing edge of the airfoil. (V_end)
			
			- Note that it is at the end of this code that the lift is calculated. 
			
				(i) wake_convect.m
				- calculates the u and v velocity components for each shed
				vortex in the wake length
				- Take into consideration some shed vortex in the wake. In order to... 
				figure
				out what velocity this vortex will convect with, we sum the influence of 
				all OTHER vortices in the wake and also add the influence of the free-
				stream
				and the imposed vortex along with the influence of the panel strengths.
				
				(ii) psiinfield_vortex.m
				This calculates the influence of all the vortices in the wake and then 
				figures out how the imposed vortex should convect downstream. 
	

Note: It seems that the formulation for the corrected wake and gamma things work out
properly when the signs and 2*pi's are introduced into the equations. 

Need to run this for a thick airfoil and then figure out what the lift response looks like.
If we do get a dip in the lift at the trailing edge 


It seems the like gam_dim values near the trailing edge are very large. As a results of
this, they are affecting the recovery of the lift in the wake. The lift coefficient 
curve is following the trend that we would expect for the unsteady lift problem
however it does not recover to the value as seen in previous works. 
