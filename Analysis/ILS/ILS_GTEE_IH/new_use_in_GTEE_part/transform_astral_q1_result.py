import re

with open('GTEE_scored.tre', 'r') as f:
    tree = f.read()

# 把 '[...]' 注释替换为 q1*100 的整数（与 RAxML 格式一致）
def replace_annotation(match):
    q1_match = re.search(r'q1=([0-9.e+-]+)', match.group(0))
    if q1_match:
        q1 = float(q1_match.group(1))
        return str(round(q1 * 100))
    return ''

cleaned = re.sub(r"'\[[^\]]*\]'", replace_annotation, tree)

with open('GTEE_raxml_format.tre', 'w') as f:
    f.write(cleaned)
