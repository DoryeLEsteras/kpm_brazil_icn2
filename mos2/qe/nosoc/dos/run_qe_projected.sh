
           export OMP_NUM_THREADS=1;mpirun --allow-run-as-root -np 20  pw.x -i            mos2.nscf.in >            mos2.nscf.out
           export OMP_NUM_THREADS=1;mpirun --allow-run-as-root -np 20  projwfc.x -i            mos2.proj.in >            mos2.proj.out
       