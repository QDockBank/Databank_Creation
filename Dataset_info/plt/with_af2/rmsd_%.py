# --*-- conding:utf-8 --*--
# @Time : 4/3/25 3:59 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : rmsd_%.py

# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

QUANTUM_FILE = "result_summary/q_rmsd_summary.txt"  # quantum方法的结果文件
AF3_FILE = "result_summary/af2_rmsd_summary.txt"  # af3方法的结果文件
OUTPUT_FILE = "all/rmsd_compare.txt"  # 对比输出文件


# =====================================================

def load_rmsd_data(file_path):
    """
    从形如:
        1e1x   1.163
        1e2k   2.835
    的文件中解析出一个 dict: { pdb_id: rmsd_val, ... }。
    跳过空行或注释行(#开头)。
    """
    data = {}
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line_s = line.strip()
            if not line_s or line_s.startswith("#"):
                continue
            cols = line_s.split()
            if len(cols) < 2:
                continue
            pdb_id = cols[0]
            try:
                rmsd_val = float(cols[1])
                data[pdb_id] = rmsd_val
            except ValueError:
                pass
    return data


def main():
    # 1) 读取 quantum 与 af3 的 RMSD 数据
    quantum_data = load_rmsd_data(QUANTUM_FILE)
    af3_data = load_rmsd_data(AF3_FILE)

    # 统计计数
    quantum_better = 0
    af3_better = 0
    tie_count = 0

    # 准备输出
    with open(OUTPUT_FILE, "w", encoding="utf-8") as out_f:
        # 写个标题(可选)
        out_f.write("# Compare RMSD fig between quantum and af3\n")
        out_f.write("# Format: pdb_id   quantum=...   af2=...   better=...\n\n")

        # 找到共同的 pdb_id 集合
        common_ids = set(quantum_data.keys()) & set(af3_data.keys())
        # 如果你想对所有出现的 pdb_id 做对比(哪家没有就记录N/A),
        # 也可用 set(quantum_data.keys()) | set(af3_data.keys())，这里示例只对共有ID

        for pdb_id in sorted(common_ids):
            q_rmsd = quantum_data[pdb_id]
            a_rmsd = af3_data[pdb_id]

            if q_rmsd < a_rmsd:
                better_method = "quantum"
                quantum_better += 1
            elif q_rmsd > a_rmsd:
                better_method = "af2"
                af3_better += 1
            else:
                better_method = "tie"
                tie_count += 1

            line_out = (f"{pdb_id}\tquantum={q_rmsd:.3f}\taf2={a_rmsd:.3f}\t"
                        f"better={better_method}\n")
            out_f.write(line_out)

        # 再写一下统计信息
        total_compared = len(common_ids)
        # 避免除0
        if total_compared == 0:
            out_f.write("\n# No common pdb_id found between two files.\n")
        else:
            q_percent = quantum_better / total_compared * 100
            af3_percent = af3_better / total_compared * 100
            tie_percent = tie_count / total_compared * 100

            out_f.write(f"\n# Total proteins compared: {total_compared}\n")
            out_f.write(f"# quantum better: {quantum_better} ({q_percent:.1f}%)\n")
            out_f.write(f"# af2 better: {af3_better} ({af3_percent:.1f}%)\n")
            if tie_count > 0:
                out_f.write(f"# tie: {tie_count} ({tie_percent:.1f}%)\n")

    print(f"对比完成！结果写入 {OUTPUT_FILE}")

if __name__ == "__main__":
    main()

