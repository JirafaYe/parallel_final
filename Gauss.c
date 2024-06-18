#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>

const int debug=0;


// 打印函数，根据 Debug 常量的值来决定是否打印
void DPrintf(const char *format, ...) {
    if (debug==1) {
        va_list args;
        va_start(args, format);
        vprintf(format, args);
        va_end(args);
        fflush(stdout);
    }
    return;
}


int main(void) {
    int n;

    double start_time = clock(); // 获取开始时间

    // 输入矩阵大小
    scanf("%d", &n);

    // 分配矩阵内存
    double **mat = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i <n; ++i) {
        mat[i] = (double*)malloc(2 * n * sizeof(double));
    }
    
    // 输入矩阵
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            scanf("%lf", &mat[i][j]);
        }
    }
    
    // 初始化右侧为单位矩阵
    for (int i = 0; i < n; ++i) {
        for (int j = n; j < 2 * n; ++j) {
            if (j == (i + n)) {
                mat[i][j] = 1;
            } else if (j >= n) {
                mat[i][j] = 0;
            }
        }
    }

    // 化简为对角阵
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                double d = mat[j][i] / mat[i][i];
                DPrintf("i[%d] j[%d] d[%lf]\n",i,j,d);
                for (int k = 0; k < 2 * n; ++k) {
                    mat[j][k] -= mat[i][k] * d;
                }
            }
        }    
    }

    // 化简为单位矩阵
    for (int i = 0; i < n; ++i) {
        double d = mat[i][i];
        for (int j = 0; j < 2 * n; ++j) {
            mat[i][j] /= d;
        }
    }

    double end_time = clock(); // 获取结束时间

    // 输出结果
    //clock返回为时钟周期数，所以除`CLOCKS_PER_SEC`求出实际秒数
    printf("process took %.6f seconds\n", ((double) (end_time - start_time)) / CLOCKS_PER_SEC);

    printf("%d\n", n);
    for (int i = 0; i < n; ++i) {
        for (int j = n; j < 2 * n; ++j) {
            printf("%lf ", mat[i][j]);
        }
        printf("\n");
    }
    
    free(mat);

    return 0;
}
