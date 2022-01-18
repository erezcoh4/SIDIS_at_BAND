piCharge = #1
echo run GEMC for ${piCharge} acceptance map

Nevents = #2
echo ${Nevents} events

/u/home/cohen/clas12Tags/4.4.1/source/gemc /u/home/cohen/clas12Tags/gcards/band-clas12.gcard -USE_GUI=0 -N=${Nevents} -INPUT_GEN_
FILE="LUND, /volatile/clas12/users/ecohen/GEMC/LUND/10.2/AcceptanceCorrection/${piCharge}/ee${piCharge}_p_uniform_distribution.dat" -OUTPUT="evio, /volatile/clas12/users/ecohen/GEMC/evio/10.2/AcceptanceCorrection/${piCharge}/ee${piCharge}_p_uniform_distribution.ev"

