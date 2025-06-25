
           export OMP_NUM_THREADS=1;mpirun --allow-run-as-root -np 10  pw.x -i            ptse2.nscf.in >            ptse2.nscf.out
           export OMP_NUM_THREADS=1;mpirun --allow-run-as-root -np 10  projwfc.x -i            ptse2.proj.in >            ptse2.proj.out
       