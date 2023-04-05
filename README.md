# PRE
 
D_G_main.m : set the dimensional trait space, calculates states variables (investments, affinities), calls an ODE and an ODE function to solve the differential equations 

D_G_ode_func.m : computes for each time step growth, mortality rates and the differential equations governing biomass and nutrients concentration evolution

D_G_physical_settings.m : computes the environmental conditions, according to the type of simulation (steady state/seasonal cycle) and the time range

daily_insolation.m : function used to have the light for each day (can vary according to the chosen latitude)

D_G_diversity.m : computes the diversity index over a seasonal cycle

D_G_extra_output.m : computes the trophic transfer efficiencies and the productivity 

D_G_plot.m : plots all the figures 
