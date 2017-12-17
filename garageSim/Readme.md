## Simulation of parking garage using Monte Carlo method
The file parkingGarageParallel.c is a simulation of cars arriving and leaving at a parking garage, parallelised using OpenMP and MPI.

### Installation and running:
1. Make sure to have MPI and OpenMP installed.
2. for building the executable, run the following command:
gcc -fopenmp -o parkingGarageParallel parkingGarageParallel.c -lm
3. For running the simulation run:
./parkingGarageParallel <simulation size> <no of stalls in the garage> <mean time between car arrivals> <mean time of stay in the garage>
