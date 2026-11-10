#!/bin/bash
#PBS -A YOUR_PROJECT
#PBS -N wmles_channel
#PBS -l walltime=01:00:00
#PBS -l select=1
#PBS -l place=scatter
#PBS -l filesystems=home:flare
#PBS -k doe
#PBS -j oe
#PBS -q debug
# ============================================================================
# Aurora (ALCF) job script: build the UDF, then run, on one node with 12 ranks
# (one per GPU tile). Submit from a run directory that holds the case files
# (channel.par/.udf/.oudf/.usr/.re2, logLawWallModel.hpp, udf.cmake):
#     qsub -A <project> /path/to/scripts/run_aurora.sh
# Set NEKRS_HOME to your nekRS v26.0 install and the cache paths to a
# location on /flare that you own.
# ============================================================================

export NEKRS_HOME=${NEKRS_HOME:-$HOME/.local/nekrs-v26.0}
CACHE_ROOT=${CACHE_ROOT:-$PBS_O_WORKDIR/cache}
export XDG_CACHE_HOME=$CACHE_ROOT
export SYCL_CACHE_DIR=$CACHE_ROOT/libsycl_cache
export OCCA_CACHE_DIR=$CACHE_ROOT/occa
export NEKRS_CACHE_DIR=$CACHE_ROOT/nekrs

cd "$PBS_O_WORKDIR"
echo "Jobid: $PBS_JOBID"
echo "Running on nodes $(cat "$PBS_NODEFILE")"

module load cmake

export NEKRS_GPU_MPI=0
export MPICH_GPU_SUPPORT_ENABLED=0
export NEKRS_CACHE_BCAST=1
export NEKRS_LOCAL_TMP_DIR=/tmp/
export OCCA_DPCPP_COMPILER_FLAGS="-O3 -fsycl -fsycl-targets=intel_gpu_pvc -ftarget-register-alloc-mode=pvc:auto -fma"
export FI_CXI_RX_MATCH_MODE=hybrid   # required by parRSB
export UR_L0_USE_COPY_ENGINE=0
export ZES_ENABLE_SYSMAN=1
export ROMIO_HINTS=.romio_hint
cat > .romio_hint <<HINT
romio_cb_write enable
romio_cb_read enable
cb_buffer_size 16777216
cb_config_list *:8
striping_unit 2097152
striping_factor 16
HINT

CASE=channel
RANKS_PER_NODE=12
CPU_BIND_LIST="1-8:9-16:17-24:25-32:33-40:41-48:53-60:61-68:69-76:77-84:85-92:93-100"
GPUS_PER_NODE=6
TILES_PER_GPU=2

# bind each rank to one GPU tile
cat > .lhelper <<HELPER
#!/bin/bash
gpu_id=\$(((PALS_LOCAL_RANKID / ${TILES_PER_GPU}) % ${GPUS_PER_NODE}))
tile_id=\$((PALS_LOCAL_RANKID % ${TILES_PER_GPU}))
export ZE_AFFINITY_MASK=\$gpu_id.\$tile_id
"\$@"
HELPER
chmod u+x .lhelper

echo "==================== BUILD ===================="
mpiexec -n ${RANKS_PER_NODE} -ppn ${RANKS_PER_NODE} --cpu-bind=list:${CPU_BIND_LIST} -- \
  ./.lhelper "$NEKRS_HOME/bin/nekrs" --setup "$CASE" --backend dpcpp --device-id 0 --build-only ${RANKS_PER_NODE}
BUILD_RC=$?
echo "BUILD_EXIT_CODE=$BUILD_RC"
[ $BUILD_RC -ne 0 ] && exit $BUILD_RC

echo "==================== RUN ===================="
mpiexec -n ${RANKS_PER_NODE} -ppn ${RANKS_PER_NODE} --cpu-bind=list:${CPU_BIND_LIST} -- \
  ./.lhelper "$NEKRS_HOME/bin/nekrs" --setup "$CASE" --backend dpcpp --device-id 0
RC=$?
echo "RUN_EXIT_CODE=$RC"
exit $RC
