nvcc -DGRB --compiler-options "-fopenmp" *.c *.cpp *.cu -I/opt/gurobi752/linux64/include/  -L/opt/gurobi752/linux64/lib/  -lgurobi75 -o chavExecGpu

