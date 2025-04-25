# --*-- conding:utf-8 --*--
# @Time : 4/3/25 4:15 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : group_rmsd.py

import os
import re


INDEX_FILE = "group_index.txt"
COMPARE_FILE = "all/rmsd_compare.txt"

OUTPUT_S_FILE = "S/rmsd_compare.txt"
OUTPUT_M_FILE = "M/rmsd_compare.txt"
OUTPUT_L_FILE = "L/rmsd_compare.txt"



def parse_group_index(index_file):
    """
    解析贴出的分组 index 文件(含 [Group S], [Group M], [Group L])，
    返回 { "S": set([...]), "M": set([...]), "L": set([...]) }
    """
    groups_map = {"S": set(), "M": set(), "L": set()}
    current_group = None

    with open(index_file, "r", encoding="utf-8") as f:
        for line in f:
            line_s = line.strip()
            if not line_s:
                continue

            # 匹配 [Group X]
            grp_match = re.match(r'^\[Group\s+([SML])\]', line_s)
            if grp_match:
                current_group = grp_match.group(1)  # 'S' or 'M' or 'L'
                continue

            # 如果行以 "length_" 开头，则是 "length_10:" 之类，跳过
            if line_s.startswith("length_"):
                continue

            # 其余我们视为 pdb_id 行
            if current_group in groups_map:
                pdb_id = line_s.split()[0]  # 取首列作为 pdb_id
                groups_map[current_group].add(pdb_id)

    return groups_map


def parse_compare_lines(compare_file):
    """
    读取对比结果文件, 行形如:
      1e1x  quantum=1.163  af3=2.043  better=quantum

    返回一个列表 compare_data:
    [ (pdb_id, q_rmsd, a_rmsd, better_str, full_line_str), ... ]
    其中 q_rmsd, a_rmsd 是 float, better_str 是 "quantum"/"af3"/"tie"
    """
    results = []

    with open(compare_file, "r", encoding="utf-8") as f:
        for line in f:
            line_s = line.strip()
            if not line_s or line_s.startswith("#"):
                continue

            # 先分拆出 pdb_id
            parts = line_s.split(None, 1)
            if len(parts) < 2:
                continue
            pdb_id = parts[0]
            rest = parts[1]

            # 从余下字符串里提取 "quantum=xxx" "af3=yyy" "better=zzz"
            # 简单用正则
            # quantum=([-\d\.]+)
            # af3=([-\d\.]+)
            # better=(\S+)
            q_pat = re.search(r'quantum=([-\d\.]+)', rest)
            a_pat = re.search(r'af3=([-\d\.]+)', rest)
            b_pat = re.search(r'better=(\S+)', rest)

            if not q_pat or not a_pat or not b_pat:
                # 如果有行缺失信息，就跳过或继续
                continue

            try:
                q_val = float(q_pat.group(1))
                a_val = float(a_pat.group(1))
            except ValueError:
                continue

            better_str = b_pat.group(1)  # "quantum" / "af3" / "tie" ...

            results.append((pdb_id, q_val, a_val, better_str, line_s))

    return results


def write_group_file(group_lines, output_path):
    """
    对该组 group_lines (列表 of (pdb_id, q_val, a_val, better_str, line_str)),
    将 line_str 写入到 output_path 文件中，
    最后统计 quantum_better, af3_better, tie 的百分比并写到文件尾。
    """
    quantum_better = 0
    af2_better = 0
    tie_count = 0

    with open(output_path, "w", encoding="utf-8") as out_f:
        # 可写个表头
        # out_f.write("# Compare fig for this group\n\n")

        for (pdb_id, q_val, a_val, better_str, full_line) in group_lines:
            out_f.write(full_line + "\n")
            # 统计
            if better_str == "quantum":
                quantum_better += 1
            elif better_str == "af3":
                af2_better += 1
            else:
                tie_count += 1

        # 统计结束
        total = len(group_lines)
        if total == 0:
            out_f.write("\n# No data in this group.\n")
            return

        qb_pct = quantum_better / total * 100
        ab_pct = af2_better / total * 100
        tie_pct = tie_count / total * 100
        out_f.write(f"\n# total={total}\n")
        out_f.write(f"# quantum better: {quantum_better} ({qb_pct:.1f}%)\n")
        out_f.write(f"# af3 better: {af2_better} ({ab_pct:.1f}%)\n")
        if tie_count > 0:
            out_f.write(f"# tie: {tie_count} ({tie_pct:.1f}%)\n")


def main():
    # 1) 解析 index => { "S": set([...]), "M": set([...]), "L": set([...]) }
    group_map = parse_group_index(INDEX_FILE)

    # 2) 解析 compare 结果 => 返回列表 [ (pdb_id, q_val, a_val, better, line_str), ... ]
    compare_data = parse_compare_lines(COMPARE_FILE)

    # 3) 分桶: S_lines, M_lines, L_lines
    S_lines = []
    M_lines = []
    L_lines = []

    for (pdb_id, q_val, a_val, better_str, line_str) in compare_data:
        if pdb_id in group_map["S"]:
            S_lines.append((pdb_id, q_val, a_val, better_str, line_str))
        elif pdb_id in group_map["M"]:
            M_lines.append((pdb_id, q_val, a_val, better_str, line_str))
        elif pdb_id in group_map["L"]:
            L_lines.append((pdb_id, q_val, a_val, better_str, line_str))
        else:
            # 不在 S/M/L => 忽略 or 另建一个 others list
            pass

    # 4) 分别写文件
    write_group_file(S_lines, OUTPUT_S_FILE)
    write_group_file(M_lines, OUTPUT_M_FILE)
    write_group_file(L_lines, OUTPUT_L_FILE)

    print("分组完成, 并在各文件统计 quantum/af2/tie 百分比.")
    print(f"S => {OUTPUT_S_FILE}")
    print(f"M => {OUTPUT_M_FILE}")
    print(f"L => {OUTPUT_L_FILE}")


if __name__ == "__main__":
    main()
