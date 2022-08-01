# Mobility-Modelling
Mobility Modeling is a smart-cities project built with C that aims to model mobility and traffic flow using probabilistic stochastic approaches
on multiple road types( Highways, intersections..etc ) to measure different measures such as average speed and risk of collision rate..etc


Here underneath is the list of files of the project with the description of each of it :

    Sparse.c: this file will generate the transition matrixs for each configuration based on the modeled markov chain
             and transform each transition matrix to a sparse matrix for perfomance purposes.
    GTH.c: this file will resolve numerically each sparse matrix generated in the previous step and calculate the stationnarry distributions of each configuration
    Average_Speed.c: This program will calculate the average speed from the stationary distributions calculated previously.
    Average_Speed.p: it is a 'GNUPLOT' file which will create graphs from the data generated in the previous step (The files will be created in the 'OUTPUTS' folder)

The program as sent varies one parameter at a time, and on this code, it varies Vc from 0 to Vmax.
