windDrivenRainFoam

An open-source solver for wind-driven rain based on OpenFOAM

More information: https://carmeliet.ethz.ch/research/downloads

Solver modifications:
1. LTS with localEulerDdtScheme
2. fvOption support in alpha and Urain equations
3. Update calculateFhd to create a DATABASE for discrete DFh and sum DFh in each phase diameter interval
previously each interval was calculated with interval's DF (best) * diameter's Vt <- Less accurate for large intervals?
4. Allow custom fhd input in transportProperties in case raindrop size distribution doesnt follow Best (1950)
