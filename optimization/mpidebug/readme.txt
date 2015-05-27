This set of functions (mpidebug/*) was developed with the aim to debug the robustMPI algorithm.
While it must work in this simple setting, they are not optimally efficient nor general.

For larger applications it is encouraged to compile wrappers as mex files to the MPI libraries as done in:
https://www.hpc.ntnu.no/display/hpc/Matlab+for+HPC

Another alternative (that I did not test) is to interface PMATLAB:
http://www.ll.mit.edu/mission/cybersec/softwaretools/pmatlab/pmatlab.html

This small set was inspired in the alternatives above.
