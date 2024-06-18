# parallel_final
并行期末大作业
## 相关指令
py生成矩阵指令

```shell
python3 script.py --n 16 --output output.txt
```

> pyhon3 generator.py --n 256 --output ./test-case/256.txt

mpi 编译：

```shell
mpicc -g -o Gauss_mpi Gauss_mpi.c
```

mpi运行

```shell
mpirun -n [proccess] ./Gauss_mpi < [test file] > [res file]
```

> mpirun -n 8 ./Gauss_mpi < ./test-case/256.txt > ./res/256_mpi_8.txt

