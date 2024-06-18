#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

const int debug = 0;

struct Result {
    double start_time, end_time;
    double *result_matrix;
    int dimension;
};

// 打印函数，根据 Debug 常量的值来决定是否打印
void DPrintf(const char *format, ...) {
    if (debug == 1) {
        va_list args;
        va_start(args, format);
        vprintf(format, args);
        va_end(args);
        fflush(stdout);
    }
    return;
}

// 截取数组片段
void sliceArray(double *source, int start, int length, double *destination) {
    for (int i = 0; i < length; i++) {
        destination[i] = source[start + i];
    }
}

// 打印数组工具
void printDoubleArray(const double *array, size_t size) {
    if (debug == 0) {
        return;
    }
    DPrintf("[");
    for (size_t i = 0; i < size; ++i) {
        DPrintf("%f", array[i]);
        if (i < size - 1)
            DPrintf(", ");
    }
    DPrintf("]\n");
}

// 初始化矩阵,并初始化右侧为单位矩阵
void initMatrix(int n, double *destination) {
    int row_len = 2 * n;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < row_len; ++j) {
            if (j < n) {
                scanf("%lf", &destination[i * row_len + j]);
            } else {
                if (j == (i + n)) {
                    destination[i * row_len + j] = 1;
                } else if (j >= n) {
                    destination[i * row_len + j] = 0;
                }
            }
        }
    }
}

// 输出结果
void printRes(struct Result res) {
    printf("MPI process took %.6f seconds\n", res.end_time - res.start_time);
    printf("%d\n", res.dimension);
    int row_len = 2 * res.dimension;
    for (int i = 0; i < res.dimension; ++i) {
        for (int j = res.dimension; j < res.dimension * 2; ++j) {
            printf("%lf ", res.result_matrix[i * row_len + j]);
        }
        printf("\n");
    }
}

int main(void) {
    int n;
    double *mat = NULL;

    MPI_Init(NULL, NULL);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // 获取起始时间
    double start_time = MPI_Wtime();

    if (rank == 0) {
        scanf("%d", &n);
        mat = (double*)malloc(2 * n * n * sizeof(double));
        initMatrix(n, mat);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int row_len = 2 * n;
    int local_rows = n / size;

    double *pivot_row = (double*)malloc(row_len * sizeof(double));
    // 选取主行并进行广播
    if (rank == 0) {
        sliceArray(mat, 0, row_len, pivot_row);
        DPrintf("init pivot by rank[0]\n");
    }
    MPI_Bcast(pivot_row, 2 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    DPrintf("after bcast rank[%d] n[%d] local_rows[%d]\n", rank, n, local_rows);

    // scatter分发矩阵
    double *local_mat = (double*)malloc((row_len) * (local_rows) * sizeof(double));
    MPI_Scatter(mat, (row_len) * local_rows, MPI_DOUBLE,
                local_mat, (row_len) * local_rows, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    DPrintf("local_mat rank[%d] n[%d]  local_rows[%d] size[%d]\n", rank, n, local_rows, size);
    printDoubleArray(local_mat, row_len * local_rows);

    // 格式化左侧为对角矩阵
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < local_rows; ++j) {
            int real_j = rank * local_rows + j;
            if (real_j != i) {
                double d = local_mat[j * row_len + i] / pivot_row[i];
                DPrintf("i[%d] real_j[%d] d[%lf]\n", i, real_j, d);
                for (int k = 0; k < row_len; ++k) {
                    local_mat[j * row_len + k] -= pivot_row[k] * d;
                }
            }
        }

        // 如果当前进程占有下一个主行则更新
        if ((i + 1) / local_rows == rank) {
            int local_row = (i + 1) % local_rows;
            sliceArray(local_mat, local_row * row_len, row_len, pivot_row);
        }

        // 同步点
        MPI_Barrier(MPI_COMM_WORLD);
        // 更新主行（由持有的线程进行广播）
        if (i != n - 1) {
            MPI_Bcast(pivot_row, row_len, MPI_DOUBLE, (i + 1) / local_rows, MPI_COMM_WORLD);
        }
    }
    DPrintf("diagonal matrix rank[%d]\n", rank);

    // 左侧格式化为单位矩阵
    for (int i = 0; i < local_rows; ++i) {
        double d = local_mat[i * row_len + (rank * local_rows + i)];
        for (int j = 0; j < row_len; ++j) {
            local_mat[i * row_len + j] = local_mat[i * row_len + j] / d;
        }
    }

    // gather
    double *result_mat;
    if (rank == 0) {
        result_mat = (double*)malloc(row_len * n * sizeof(double));
    }

    MPI_Barrier(MPI_COMM_WORLD);
    DPrintf("rank[%d]:\n", rank);
    MPI_Gather(local_mat, row_len * local_rows, MPI_DOUBLE,
               result_mat, row_len * local_rows, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    // 获取结束时间
    double end_time = MPI_Wtime();

    MPI_Finalize();

    struct Result result = {start_time, end_time, result_mat, n};

    if (rank == 0) {
        printRes(result);
        free(result_mat);
    }

    free(local_mat);
    if (rank == 0) {
        free(mat);
    }

    return 0;
}
