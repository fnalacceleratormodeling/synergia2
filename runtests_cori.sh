#!/bin/bash -l
 
#SBATCH -q debug    
#SBATCH -C knl,quad,cache
#SBATCH -N 1
#SBATCH -t 00:20:00 
#SBATCH -J runtests
#SBATCH -o runtests.out
#SBATCH -e runtests.err
 
echo "allocated nodes: " $SLURMD_NODENAME

#. ~/cori/load_gcc_synergia_modules.sh
#. ../../setup.sh

. ~/cori/load_synergia_modules_intel_mic_2019.sh
module unload darshan

cd src/synergia/foundation/tests
srun -n 1 ./test_four_momentum
srun -n 1 ./test_physical_constants
srun -n 1 ./test_reference_particle
srun -n 1 ./test_distribution

cd ../../utils/tests
srun -n 1 ./test_command_line_arg
srun -n 1 ./test_logger
srun -n 1 ./test_complex_error_function
srun -n 1 ./test_floating_point
srun -n 1 ./test_fast_int_floor

srun -n 1 ./test_commxx_mpi
srun -n 2 ./test_commxx_mpi
srun -n 4 ./test_commxx_mpi
srun -n 7 ./test_commxx_mpi
srun -n 12 ./test_commxx_mpi

srun -n 1 ./test_commxx_divider_mpi
srun -n 4 ./test_commxx_divider_mpi

srun -n 1 ./test_logger_mpi
srun -n 2 ./test_logger_mpi
srun -n 3 ./test_logger_mpi

srun -n 1 ./test_hdf5_file
srun -n 1 ./test_hdf5_writer
srun -n 1 ./test_hdf5_serial_writer
srun -n 1 ./test_hdf5_chunked_array2d_writer

srun -n 1 ./test_distributed_fft3d
srun -n 1 ./test_distributed_fft3d_mpi
srun -n 2 ./test_distributed_fft3d_mpi
srun -n 3 ./test_distributed_fft3d_mpi
srun -n 4 ./test_distributed_fft3d_mpi

srun -n 1 ./test_distributed_fft2d
srun -n 1 ./test_distributed_fft2d_mpi
srun -n 2 ./test_distributed_fft2d_mpi
srun -n 3 ./test_distributed_fft2d_mpi
srun -n 4 ./test_distributed_fft2d_mpi

srun -n 1 ./test_xml_serialization

srun -n 1 ./test_multi_array_offsets

srun -n 1 ./test_multi_array_serialization

srun -n 1 ./test_lsexpr

cd ../../lattice/tests

srun -n 1 ./test_chef_utils
srun -n 1 ./test_lattice_element
srun -n 1 ./test_lattice_element_slice
srun -n 1 ./test_lattice
srun -n 1 ./test_lattice_diagnostics
srun -n 1 ./test_lattice_diagnostics_mpi
srun -n 2 ./test_lattice_diagnostics_mpi
srun -n 1 ./test_element_adaptor
srun -n 1 ./test_mad8_adaptor_map
srun -n 1 ./test_madx_adaptor_map
srun -n 1 ./test_chef_lattice
srun -n 1 ./test_chef_lattice_section
srun -n 1 ./test_madx_parser
srun -n 1 ./test_madx_reader

cd ../../collective/tests

srun -n 1 ./test_rectangular_grid_domain
srun -n 1 ./test_deposit
srun -n 1 ./test_deposit_xyz
srun -n 1 ./test_rectangular_grid
srun -n 1 ./test_distributed_rectangular_grid
srun -n 1 ./test_distributed_rectangular_grid_mpi
srun -n 2 ./test_distributed_rectangular_grid_mpi
srun -n 3 ./test_distributed_rectangular_grid_mpi
srun -n 4 ./test_distributed_rectangular_grid_mpi

srun -n 1 ./test_interpolate_rectangular_zyx
srun -n 1 ./test_interpolate_rectangular_xyz
srun -n 1 ./test_space_charge_3d_open_hockney
srun -n 1 ./test_space_charge_3d_open_hockney_mpi
srun -n 2 ./test_space_charge_3d_open_hockney_mpi
srun -n 3 ./test_space_charge_3d_open_hockney_mpi
srun -n 4 ./test_space_charge_3d_open_hockney_mpi
srun -n 1 ./test_space_charge_3d_open_hockney2
srun -n 1 ./test_space_charge_3d_open_hockney3
srun -n 1 ./test_space_charge_3d_open_hockney4

srun -n 1 ./test_space_charge_2d_open_hockney
srun -n 1 ./test_space_charge_2d_open_hockney_mpi
srun -n 2 ./test_space_charge_2d_open_hockney_mpi
srun -n 3 ./test_space_charge_2d_open_hockney_mpi
srun -n 4 ./test_space_charge_2d_open_hockney_mpi
srun -n 1 ./test_space_charge_2d_open_hockney2
srun -n 1 ./test_space_charge_2d_open_hockney2_mpi
srun -n 2 ./test_space_charge_2d_open_hockney2_mpi
srun -n 3 ./test_space_charge_2d_open_hockney2_mpi
srun -n 4 ./test_space_charge_2d_open_hockney2_mpi
srun -n 1 ./test_space_charge_2d_open_hockney3
srun -n 1 ./test_space_charge_2d_open_hockney4
srun -n 1 ./test_space_charge_2d_bassetti_erskine
srun -n 1 ./test_impedance
srun -n 2 ./test_impedance
srun -n 3 ./test_impedance
srun -n 4 ./test_impedance
srun -n 1 ./test_wake_field

