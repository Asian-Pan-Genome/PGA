import argparse
import sys
from Bio import SeqIO

def main():
    # 1. 设置命令行参数解析
    parser = argparse.ArgumentParser(
        description="过滤掉蛋白质序列内部（中间）含有 '*' 的序列，允许末尾存在 '*'"
    )
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="输入的 FASTA 文件路径 (Input FASTA file)"
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="输出的干净 FASTA 文件路径 (Output FASTA file)"
    )
    
    args = parser.parse_args()

    clean_records = []
    total_count = 0
    removed_count = 0

    print(f"正在读取并处理文件: {args.input} ...")

    # 2. 解析并过滤序列
    try:
        for record in SeqIO.parse(args.input, "fasta"):
            total_count += 1
            seq_str = str(record.seq)
            
            # 核心逻辑：切片去掉最后一个字符，检查中间是否含有 '*'
            if "*" not in seq_str[:-1]:
                clean_records.append(record)
            else:
                removed_count += 1

        # 3. 写入新文件
        SeqIO.write(clean_records, args.output, "fasta")

        # 4. 打印统计报告
        print("\n" + "="*30)
        print("处理完成！统计报告如下：")
        print(f"• 原始序列总数: {total_count} 条")
        print(f"• 过滤异常序列: {removed_count} 条 (内部含 *)")
        print(f"• 保留干净序列: {len(clean_records)} 条")
        print(f"• 结果已保存至: {args.output}")
        print("="*30)

    except FileNotFoundError:
        print(f"错误：找不到输入文件 '{args.input}'，请检查路径是否正确。", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"运行过程中发生错误: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
