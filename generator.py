import numpy as np
import argparse

def generate_invertible_matrix(n, output_file):
    """
    生成一个n x n的可逆矩阵,并将结果输出到指定文件中。
    
    参数:
    n (int): 矩阵的大小
    output_file (str): 输出文件的路径
    """
    # 生成随机矩阵
    m = np.random.rand(n, n)
    
    # 在对角线上放置行和的绝对值
    mx = np.sum(np.abs(m), axis=1)
    np.fill_diagonal(m, mx)
    
    # 将结果输出到文件
    with open(output_file, 'w') as f:
        f.write(f"{n}\n")
        np.savetxt(f, m, fmt='%.4f')
        # f.write(f"\nn = {n}\n")
    
    print(f"矩阵已成功输出到 {output_file}")
    return m

if __name__ == "__main__":
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='生成可逆矩阵并输出到文件')
    parser.add_argument('--n', type=int, required=True, help='矩阵的大小')
    parser.add_argument('--output', type=str, required=True, help='输出文件路径')
    args = parser.parse_args()
    
    # 调用函数生成矩阵
    generate_invertible_matrix(args.n, args.output)