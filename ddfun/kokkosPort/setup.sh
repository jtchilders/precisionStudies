# source /vast/projects/datascience/parton/conda/2025-01/bin/activate
module load cuda/12.3.0 gcc/12.2.0 cmake/3.28.3

# setup Kokkos
# BASEPATH=/vast/projects/datascience/parton/kokkos/kokkos-4.5.00/Kokkos_ARCH_VOLTA70/Release
BASEPATH=/vast/projects/datascience/parton/kokkos/kokkos-4.5.01/Kokkos_ARCH_AMPERE80/Release
export INSTPATH=install
export KOKKOS_HOME=$BASEPATH/kokkos/$INSTPATH
export KOKKOSKERNELS_HOME=$BASEPATH/kokkos-kernels/$INSTPATH
export CMAKE_PREFIX_PATH=$KOKKOS_HOME/lib64/cmake/Kokkos:$KOKKOSKERNELS_HOME/lib64/cmake/KokkosKernels
export LD_LIBRARY_PATH=$KOKKOS_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$KOKKOS_HOME/bin:$PATH