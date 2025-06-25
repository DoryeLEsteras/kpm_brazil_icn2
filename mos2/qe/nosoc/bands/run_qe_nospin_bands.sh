
       export OMP_NUM_THREADS=1;mpirun --allow-run-as-root -np 20  pw.x -i            mos2.bands.in >             mos2.bands.out
       export OMP_NUM_THREADS=1;mpirun --allow-run-as-root -np 20  bands.x -i            mos2.bs.in >            mos2.bs.out
       