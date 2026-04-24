# multispring-spmv-kernels

This repository explains the sparse matrix-vector product kernels in the Artifact Description (AD)/ Artifact Evaluation (AE) appendix for paper "Heterogeneous Computing for High-Throughput 3D Ensemble Simulations: Efficient Resource Utilization via CPU-GPU Co-execution" submitted for SC26 Technical Paper.

## Description

This paper addresses the pressing HPC need for scalable and energy-efficient ensemble workloads in the "AI for Science" era. While modern supercomputers feature heterogeneous architectures, the prevailing GPU-centric paradigm severely underutilizes CPU cores and CPU memory resources. We tackle this by proposing a framework that leverages modern high-bandwidth CPU-GPU interconnects to maximize full-node resource utilization, achieving dramatic improvements in both time-to-solution and energy-to-solution for power-constrained next-generation supercomputing. Here, we introduce a novel class of adaptive CPU-GPU co-execution algorithms for non-linear time-evolution problems. Moving beyond the traditional view of CPU-GPU data transfer as a mere bottleneck, our algorithm integrally manages memory hierarchies to distribute computational loads. To ensure performance portability across rapidly evolving hardware, we propose multiple load distribution schemes together with a performance model for choosing the suitable scheme. By demonstrating a 2.79x speedup and 2.1x higher energy efficiency over standard GPU-only implementations in an end-to-end neural network inverse analysis workflow, this paper provides a robust algorithmic foundation for complex scientific applications on modern heterogeneous systems.

The findings of this paper were attained by a software implementation and evaluation on computer environments. The matrix-vector product kernels in this repository is used for the evaluation of Table IV in the paper, which are the most computationally costly kernels in the baseline and proposed methods. Kernels not included in Table IV are memory bandwidth bound and take less time than the matrix-vector product kernel in both the baseline and proposed methods; thus, the performance trends of the entire application programs can be estimated by the kernels in this repository.

Below we explain the details of the sparse matrix-vector product kernels.

## Requirement
* GH200 node with nvhpc complier

## Source codes
* `main.F90`: main code
* `kernelCRS.F`: CRS-based kernel code in OpenMP/OpenACC
* `kernelEBE.F`: EBE-based kernel code in OpenACC
* `Makefile`: Makefile used for building program

## Compile
Code can be compiled by typing below in the source directory.
```
make clobber
make all
```
Three executable files will be generated as below
```
./kernelCRS_CPU_OPENMP.exe
./kernelCRS_GPU_OPENACC.exe
./kernelEBE_GPU_OPENACC.exe
```

## Test data
Test data can be downloaded from https://www.eri.u-tokyo.ac.jp/cshpc/share/SC26/data.tar.gz (6.0 GB) and extracted as `tar xzf data.tar.gz`. This dataset is based on a finite-element mesh with 7,759,800 second-order tetrahedral elements and 32,425,518 degrees of freedom (DOF). Since we cannot post the proprietary dataset of actual site used in the paper, we used a synthetic dataset with similar characteristics. Contents of `./data` is as below.
* `./data/setting.dat`: setting data of finite-element model (number of nodes, number of elements, number of materials, number of nonzero CRS components).
* `./data/material.dat`: material properties of finite-element model.
* `./data/conn.bin`: connectivity data of finite-element model.
* `./data/coor.bin`: coordinate data of finite-element model.
* `./data/num.bin`: material number of each element of finite-element model.
* `./data/crsptr.bin`: pointer array for matrix in CRS format.
* `./data/crsind.bin`: index array for matrix in CRS format.
* `./data/crsval.bin`: value array for matrix in CRS format.

## Run
Kernels programs can be run as below.
```
export OMP_NUM_THREADS=70
./kernelCRS_CPU_OPENMP.exe
./kernelCRS_GPU_OPENACC.exe
./kernelEBE_GPU_OPENACC.exe
```

Results measured on a Single-GH200 node with one 72-core ARMv9-a Grace CPU@3.1 GHz, 480 GB LPDDR5X memory, and one H100 96 GB GPU, with nvhpc/24.1, NVIDIA Driver Version 535.104.05, CUDA Version 12.2 on Ubuntu 22.04.4 LTS is shown below:

==== ./kernelCRS_CPU_OPENMP.exe ====
```
 using CRS with NVEC            1
 using OpenMP
 matvec took   8.0765962600708008E-002
 matvec took   7.8464984893798828E-002
 matvec took   7.8108072280883789E-002
 matvec took   7.7730178833007813E-002
 matvec took   7.7975034713745117E-002
 matvec took   7.7946186065673828E-002
 matvec took   7.8166961669921875E-002
 matvec took   7.7924013137817383E-002
 matvec took   7.7954053878784180E-002
 matvec took   7.7886104583740234E-002
 norm            0    5661.157519040894
 norm            1    1211885855437854.
```

==== ./kernelCRS_GPU_OPENACC.exe ====
```
 using CRS with NVEC            1
 using OpenACC
 matvec took   1.0486125946044922E-002
 matvec took   1.0435819625854492E-002
 matvec took   1.0427951812744141E-002
 matvec took   1.0430812835693359E-002
 matvec took   1.0428905487060547E-002
 matvec took   1.0430812835693359E-002
 matvec took   1.0431051254272461E-002
 matvec took   1.0428905487060547E-002
 matvec took   1.0434150695800781E-002
 matvec took   1.0432958602905273E-002
 norm            0    5661.157519041320
 norm            1    1211885855437870.
```

==== ./kernelEBE_GPU_OPENACC.exe ====
```
 using EBE with NVEC            1
 using OpenACC
 matvec took   8.4710121154785156E-003
 matvec took   8.3971023559570313E-003
 matvec took   8.3909034729003906E-003
 matvec took   8.3870887756347656E-003
 matvec took   8.3990097045898438E-003
 matvec took   8.3999633789062500E-003
 matvec took   8.3971023559570313E-003
 matvec took   8.3949565887451172E-003
 matvec took   8.3930492401123047E-003
 matvec took   8.3999633789062500E-003
 norm            0    5661.157519041320
 norm            1    1211885874177068.
```

## Validation of results

The correctness of results can be checked by comparing the `norm` values.

## Interpretation of results

Elapsed time for the matrix vector product is given by the number `matvec took` in seconds. 10 iterations are conducted to check the fluctuation in elapsed time. As the dataset provided for the AD/AE Appendix is not exactly the same dataset used in Table IV in the paper, the measurements are slightly different, however, we can grasp the overall performance characteristics in which the SpMV elapsed is reduced by use of GPU and with EBE-based methods.  The FLOP count and memory transfer size can be measured using nvprof, and scaled with the measured values without nvprof to obtain values in Table IV.

## License

multispring-spmv-kernels, version 1.0.0 (c) 2026 Tsuyoshi Ichimura et al. multispring-spmv-kernels is freely distributable under the terms of an MIT-style license.
