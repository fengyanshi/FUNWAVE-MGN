#
# Template:  OpenMPI, High-Bandwidth (Infinipath PSM) Variant
# Revision:  $Id: openmpi-psm.qs 483 2014-01-21 21:07:12Z frey $
#
# Usage:
# 1. Modify "NPROC" in the -pe line to reflect the number
#    of processors desired
# 2. Modify the value of "MY_EXE" to be your MPI program and any
#    arguments to be passed to it.
# 3. Uncomment the WANT_CPU_AFFINITY line if you want Open MPI to
#    bind workers to processor cores.  Can increase your program's
#    efficiency.
# 4. Uncomment the SHOW_MPI_DEBUGGING line if you want very verbose
#    output written to the Grid Engine output file by OpenMPI.
# 
#$ -pe openmpi 100
#$ -l  psm_endpoints=1
#
#$ -l standby=1
#
# If you want an email message to be sent to you when your job ultimately
# finishes, edit the -M line to have your email address and change the
# next two lines to start with #$ instead of just #
# -m eas
# -M fyshi@udel.edu
#

#
# Setup the environment; choose the OpenMPI version that's
# right for you:
#
source /opt/shared/valet/docs/valet.sh
vpkg_require openmpi/1.4.4-gcc

#
# The MPI program to execute and any arguments to it:
#
MY_EXE="/lustre/work/kirby/btrad/CHOI/05_tslevelsMS4/src/mynest"

#
# Ask Open MPI to do processor binding?
#
WANT_CPU_AFFINITY=YES

#
# Should the job use just half of the cores allocated to it?
# (Numerically-heavy jobs on Mills may perform better using
# this option.)
#
#WANT_HALF_CORES_ONLY=YES

#
# Let's try to intelligently control how many PSM contexts
# get used by Open MPI jobs.  Uncomment the following line
# if you DO NOT want to try to control context usage.
#
#DISABLE_PSM_ADJUSTMENTS=YES

#
# Uncomment to enable lots of additional information as OpenMPI executes
# your job:
#
#SHOW_MPI_DEBUGGING=YES

##
## You should NOT need to change anything after this comment.
##
OPENMPI_FLAGS="--display-map --mca btl ^tcp"
if [ "x$WANT_CPU_AFFINITY" = "xYES" ]; then
  OPENMPI_FLAGS="${OPENMPI_FLAGS} --bind-to-core"
fi
if [ "x$WANT_HALF_CORES_ONLY" = "xYES" ]; then
  OPENMPI_FLAGS="${OPENMPI_FLAGS} --cpus-per-proc 2 --np $((NSLOTS/2)) --loadbalance"
fi
if [ "x$SHOW_MPI_DEBUGGING" = "xYES" ]; then
  OPENMPI_FLAGS="${OPENMPI_FLAGS} --debug-devel --debug-daemons --display-devel-map --display-devel-allocation --mca mca_verbose 1 --mca coll_base_verbose 1 --mca ras_base_verbose 1 --mca ras_gridengine_debug 1 --mca ras_gridengine_verbose 1 --mca btl_base_verbose 1 --mca mtl_base_verbose 1 --mca plm_base_verbose 1 --mca pls_rsh_debug 1"
  if [ "x$WANT_CPU_AFFINITY" = "xYES" -o "x$WANT_HALF_CORES_ONLY" = "xYES" ]; then
    OPENMPI_FLAGS="${OPENMPI_FLAGS} --report-bindings"
  fi
fi
## For PSM, let's smartly adjust the PSM runtime environment:
if [ ! "x$DISABLE_PSM_ADJUSTMENTS" = "xYES" ]; then
  OUT_PSM_SHAREDCONTEXTS_MAX=`/opt/shared/openmpi/extras/bin/psm-max-sharedcontexts`
  if [ $? = 0 ]; then
    echo "Adjusted maximum PSM contexts per node = $OUT_PSM_SHAREDCONTEXTS_MAX"
    PSM_SHAREDCONTEXTS_MAX=$OUT_PSM_SHAREDCONTEXTS_MAX
    export PSM_SHAREDCONTEXTS_MAX
  fi
fi

echo "GridEngine parameters:"
echo "  mpirun        = "`which mpirun`
echo "  nhosts        = $NHOSTS"
echo "  nproc         = $NSLOTS"
echo "  executable    = $MY_EXE"
echo "  MPI flags     = $OPENMPI_FLAGS"
echo "-- begin OPENMPI run --"
mpirun ${OPENMPI_FLAGS} $MY_EXE
echo "-- end OPENMPI run --"
