#source /site/12gev_phys/softenv.csh 2.4
#source /group/clas12/packages/setup.csh
#module load sqlite/4.4.1
#setenv CCDB_CONNECTION sqlite:////volatile/clas12/users/tkutz/GEMC/sqlite_databases/ccdb_latest.sqlite
#setenv RCDB_CONNECTION sqlite:////volatile/clas12/users/tkutz/GEMC/sqlite_databases/rcdb_latest.sqlite
#setenv GEMC /u/home/cohen/clas12Tags/4.4.1/source
#setenv GEMC_DATA_DIR /u/home/cohen/clas12Tags/4.4.1
#setenv GEMC_DATA_DIR /work/clas12/users/tkutz/gemc/clas12Tags/4.4.1
#setenv GEMC_VERSION 4.4.1
#setenv FIELD_DIR /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/soft/noarch/data/magfield/ascii/

source /group/clas12/packages/setup.csh
module load clas12/pro

cd SIDIS_at_BAND
pwd
git pull --no-edit

echo "done pulling from repo, ready for SIDIS at BAND"
ll
