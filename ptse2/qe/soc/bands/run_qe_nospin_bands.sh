
       export OMP_NUM_THREADS=1;mpirun --allow-run-as-root -np 10  pw.x -i            ptse2.bands.in >             ptse2.bands.out
       export OMP_NUM_THREADS=1;mpirun --allow-run-as-root -np 10  bands.x -i            ptse2.bs.in >            ptse2.bs.out
       