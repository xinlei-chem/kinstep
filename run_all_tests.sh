#!/usr/bin/env bash
# KinSTEP 全量回归：8 套数值测试 + CLI 五个子命令 + 错误路径。
# 任何改动之后跑这个；它**不止看 exit code**，还抓关键数字是否退化。
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${KINSTEP_PY:-$HOME/.venvs/kinstep/bin/python}"
KS="${KINSTEP_BIN:-$HOME/.venvs/kinstep/bin/kinstep}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
# 失败日志写到 notes/ 下。★ notes/ 不随存档发布（见 RELEASE_NOTES.md 第四节），
#   所以这里必须自己建，否则干净下载后第一次失败时 cp 会再报一个无关的错。
mkdir -p "$ROOT/notes"
export OMP_NUM_THREADS=1
pass=0; fail=0; skipped=0
ok(){ printf "  ✓ %-42s %s\n" "$1" "${2:-}"; pass=$((pass+1)); }
no(){ printf "  ✗ %-42s %s\n" "$1" "${2:-}"; fail=$((fail+1)); }
# ⊘ = 跳过，不计入通过也不计入失败。退出码 77 沿用 automake 的"跳过"约定：
#   依赖第三方参考数据（CATKINAS，闭源、未授再分发，故不随本包分发）的几套
#   回归在数据缺失时跳过。详见 RELEASE_NOTES.md 与 tests/refdata.py。
sk(){ printf "  ⊘ %-42s %s\n" "$1" "${2:-}"; skipped=$((skipped+1)); }

echo "══ ① 环境 ══"
"$PY" - <<'PYEOF'
import sys, numpy, scipy
print("  python %s | numpy %s | scipy %s" % ('.'.join(map(str,sys.version_info[:3])), numpy.__version__, scipy.__version__))
for name, mod in (("mpmath（高精度，可选）","mpmath"), ("matplotlib（作图，可选）","matplotlib")):
    try:
        m=__import__(mod); print("  %-24s %s" % (name, m.__version__))
    except ImportError:
        print("  %-24s 未装（会降级，不报错）" % name)
PYEOF

echo
echo "══ ①b 全包语法扫描：中文串里的引号/反斜杠不能把代码写坏 ══"
# 两次踩过：① 中文串里混进 ASCII 双引号；② raw-string 补丁写进非 raw 串，
#   `\times` 里的 \t 被 Python 解析成制表符。两类都能被「警告当错误」抓住。
sweep=$("$PY" -I - "$ROOT" << 'PYEOF'
import sys, pathlib, warnings, ast
bad = []
root = pathlib.Path(sys.argv[1])
# 只扫自己的代码：refdata_local/ 下是第三方（CatMAP 上游副本）的参考代码，不归我们管
mine = ["kinstep", "测试", "案例", "论文"]
files = [f for d in mine for f in sorted((root / d).rglob("*.py"))]
files += sorted(root.glob("*.py"))
for f in files:
    if "__pycache__" in str(f):
        continue
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ast.parse(f.read_text(encoding="utf-8"), filename=str(f))
    except Exception as e:
        bad.append("%s: %s" % (f.name, str(e).split("\n")[0][:70]))
print("\n".join(bad) if bad else "OK")
PYEOF
)
if [ "$sweep" = "OK" ]; then
    ok "  全部 .py 语法干净（含转义警告）"
else
    no "  语法/转义问题：$(echo "$sweep" | head -3 | tr '\n' ' ')"
fi
echo
echo "══ ② 数值回归（15 套，对真值逐数字）══"
for t in t3b_prefix t4_test t5_koper_cv t6_steady t7_catkinas_401 t8_drc t9_scan t10_thermo t11_highprec t12_photo t13_sono t9b_flow t14_templates t15_nl t16_i18n t17_field; do
    log="$TMP/$t.log"
    # ★ 只看 exit code 不够：测试脚本会打印 ❌ 但照样 exit 0（静默失败）。
    #   实测教训：T5 的右端相对差从 2.4e-15 退化到 2.6e-06，exit code 全程是 0。
    #   ⇒ 输出里出现 ❌ / ✗ 一律判失败。这条覆盖所有测试，不必逐个改脚本。
    timeout 300 "$PY" "$ROOT/tests/$t.py" > "$log" 2>&1; rc=$?
    # ★ 原来写的是 `if ! cmd; then ... exit=$?`，而 `!` 会把 $? 翻成 0，
    #   所以失败信息里永远印 exit=0。先取 rc 再判断。
    if [ "$rc" = 77 ]; then
        why="$(grep -m1 -oE "⊘ 跳过本套回归：.*" "$log" | sed "s/⊘ 跳过本套回归：//")"
        sk "$t" "${why:-见 $log}"
    elif [ "$rc" != 0 ]; then
        no "$t" "exit=$rc 见 $log"; cp "$log" "$ROOT/notes/FAILED_$t.log"
    elif grep -qE '❌|✗|Traceback' "$log"; then
        no "$t" "输出里有 ❌/✗：$(grep -oE '[^ ]*❌[^，]*|[^ ]*✗[^，]*' "$log" | head -1)"
        cp "$log" "$ROOT/notes/FAILED_$t.log"
    else
        ok "$t" "$(grep -oE '最大[^，；]*[0-9]\.[0-9]+e[-+][0-9]+' "$log" | tail -1)"
    fi
done

echo
echo "══ ③ 关键数字有没有退化（只看 exit code 不够）══"
chk(){ # 名称 日志 正则 阈值 [hp]
    # 所属回归因缺第三方参考数据而跳过时，这里没有数字可抓——那是跳过，不是退化。
    # 不加这条就会把"跳过"误报成"没抓到数字"的失败。
    if grep -q "⊘ 跳过本套回归" "$TMP/$2.log" 2>/dev/null; then
        # 原因照抄日志里的那一行：现在既可能是缺参考数据，也可能是缺可选件，
        # 写死成"缺第三方参考数据"会误导读报告的人。
        w="$(grep -m1 -oE "⊘ 跳过本套回归：.*" "$TMP/$2.log" | sed "s/⊘ 跳过本套回归：//")"
        sk "$1" "所属回归已跳过（${w:-原因见日志}）"; return; fi
    # ★ 第 5 个参数 hp：这条判据的阈值**以高精度升级为前提**。强吸附支上 1−z 低于
    #   double 地板，没有 mpmath 就只能停在 1e-2 量级 —— 那是可选件没装，不是退化。
    #   （参考值随包分发之前，这两条总是先因缺 .mat 而跳过，这个分支走不到。）
    if [ "${5:-}" = hp ] && ! "$PY" -c "import mpmath" >/dev/null 2>&1; then
        sk "$1" "判据以高精度为前提，mpmath 未装（可选件）"; return; fi
    v=$(grep -oE "$3" "$TMP/$2.log" 2>/dev/null | tail -1 | grep -oE '[0-9]\.[0-9]+e[-+][0-9]+' | tail -1)
    if [ -z "$v" ]; then no "$1" "没抓到数字"; return; fi
    if "$PY" -c "import sys;sys.exit(0 if float('$v')<=$4 else 1)"; then ok "$1" "$v ≤ $4"
    else no "$1" "$v > $4 ← 退化了"; fi
}
chk "T5 Koper 右端逐点"      t5_koper_cv      '右端最大相对差 = [0-9.e+-]+'      1e-13
chk "T6 稳态 vs 闭式解"       t6_steady        '相对差 [0-9.e+-]+'               1e-14
chk "T7 覆盖度 vs CATKINAS"   t7_catkinas_401  '覆盖度 / 可逆度 最大相对差 = [0-9.e+-]+' 1e-9
chk "T8 DRC 求和规则偏差"     t8_drc           '求和规则偏差 [0-9.e+-]+'          1e-9
chk "T9 全链路 Rnet"          t9_scan          '净速率 Rnet  最大 [0-9.e+-]+'     1e-8   hp
chk "T10 熵校正 5 个常数"      t10_thermo       'ΔG⁰ 最大差 [0-9.e+-]+'            1e-8
chk "T10 全链路 Rnet"         t10_thermo       '净速率最大相对差 [0-9.e+-]+'      1e-7   hp
chk "T11 高精度 401 点 Rnet"  t11_highprec     '净速率 Rnet 最大 [0-9.e+-]+'      1e-8   hp
chk "T12 光催化电荷守恒"      t12_photo        '最大电荷失衡 = [0-9.e+-]+'        1e-10
chk "T13 超声份额反推"        t13_sono         '相对差 [0-9.e+-]+'                1e-9
chk "T9b 流动出口分压中位"    t9b_flow         '出口分压 y.5:8.[ ]+最大 [0-9.e+-]+[ ]+中位 [0-9.e+-]+' 1e-12
chk "T9b 覆盖度中位"          t9b_flow         '其它量：覆盖度中位 [0-9.e+-]+'    1e-12

echo
echo "══ ④ CLI 五个子命令 ══"
cd "$TMP"
"$KS" --help >/dev/null 2>&1 && ok "kinstep --help" || no "kinstep --help"
"$KS" drivers >/dev/null 2>&1 && ok "kinstep drivers" || no "kinstep drivers"
"$KS" templates >/dev/null 2>&1 && ok "kinstep templates" || no "kinstep templates"
"$KS" templates oer >/dev/null 2>&1 && ok "kinstep templates oer" || no "kinstep templates oer"
"$KS" new --template her t.kinstep.toml >/dev/null 2>&1 && ok "kinstep new --template" || no "kinstep new --template"
"$KS" ask "析氢反应 pH=0 电位 -0.1 V vs SHE" >/dev/null 2>&1 && ok "kinstep ask" || no "kinstep ask"
"$KS" ask "析氧反应 1.5 V vs SHE" --write ask_out.kinstep.toml >/dev/null 2>&1 && ok "kinstep ask --write" || no "kinstep ask --write"
"$KS" check ask_out.kinstep.toml >/dev/null 2>&1; [ $? = 3 ] && ok "  ask 产物 check 正确报缺能量(exit=3)" || no "  ask 产物 check"
"$KS" ask "析氢反应 pH=0 电位 -0.1 V vs RHE" 2>&1 | grep -q "重复计算" && ok "  ask 拦住 pH 双重计数" || no "  ask pH 双重计数"
"$KS" new a.kinstep.toml >/dev/null 2>&1 && ok "kinstep new" || no "kinstep new"
"$KS" new b.kinstep.toml --scan >/dev/null 2>&1 && ok "kinstep new --scan" || no "kinstep new --scan"
"$KS" check a.kinstep.toml >/dev/null 2>&1 && ok "kinstep check（单点）" || no "kinstep check（单点）"
"$KS" check b.kinstep.toml >/dev/null 2>&1 && ok "kinstep check（扫描）" || no "kinstep check（扫描）"
"$KS" explain a.kinstep.toml >/dev/null 2>&1 && ok "kinstep explain" || no "kinstep explain"
timeout 180 "$KS" run a.kinstep.toml --drc >"$TMP/run1.log" 2>&1 && ok "kinstep run --drc" || no "kinstep run --drc"
grep -q "求和规则 ΣX_B = 1.0000" "$TMP/run1.log" && ok "  run 的 DRC 求和规则 = 1" || no "  run 的 DRC 求和规则"
timeout 180 "$KS" run b.kinstep.toml >"$TMP/run2.log" 2>&1 && ok "kinstep run（401 点扫描）" || no "kinstep run（扫描）"
grep -q "收敛 401/401" "$TMP/run2.log" && ok "  扫描 401/401 全收敛" || no "  扫描收敛率"
[ -s b.kinstep_out/scan.dat ] && ok "  导出原始数据 scan.dat" || no "  scan.dat"
# matplotlib 是**可选**件（pyproject 的 plot extra）；没装就不出图，这是设计如此，
# 不是缺陷。所以这里按"有没有装"分流——缺可选件而报红是误报。
# 实测：干净 sdist 装进新 venv、只有 numpy+scipy 时，原来这一项会红。
if "$PY" -c "import matplotlib" >/dev/null 2>&1; then
    [ -s b.kinstep_out/volcano.png ] && ok "  出图 volcano.png" || no "  volcano.png"
else
    sk "  出图 volcano.png" "未装可选件 matplotlib"
fi

echo
echo "══ ④b 六个案例必须都能跑（README/EXAMPLES.md 里的数字都来自它们）══"
mkdir -p "$TMP/cases" && cd "$TMP/cases"
for c in "$ROOT"/examples/ex*.py; do
    nm="$(basename "$c" .py)"
    timeout 200 "$PY" "$c" > "$TMP/$nm.log" 2>&1; rc=$?
    if [ "$rc" = 0 ]; then
        if grep -qE "❌|Traceback" "$TMP/$nm.log"; then
            no "$nm" "输出里有 ❌/Traceback"
        else
            ok "$nm" "$(grep -oE '[0-9]\.[0-9]+e[-+][0-9]+' "$TMP/$nm.log" | tail -1)"
        fi
    else
        no "$nm" "exit=$rc 见 $TMP/$nm.log"; cp "$TMP/$nm.log" "$ROOT/notes/FAILED_$nm.log"
    fi
done
cd "$TMP"

echo
echo "══ ④c 文档里引用的关键数字不能漂（README / EXAMPLES.md）══"
dchk(){ grep -qF "$2" "$ROOT/$1" && ok "  $1 含「$2」" || no "  $1 缺「$2」"; }
dchk README.md "2.4e-15"
dchk README.md "1.4e-10"
dchk README.md "2.2e-16"
dchk EXAMPLES.md "118.3 mV/dec"
dchk EXAMPLES.md "1.394e+23"
dchk EXAMPLES.md "1.793e-27"
# ★ 论文里的数字同样不能漂：代码一改、论文就该跟着改，不能只改一头
dchk paper/manuscript.md "2.44\\times10^{-15}"
dchk paper/manuscript.md "1.41\\times10^{-10}"
dchk paper/manuscript.md "2.71\\times10^{-10}"
dchk paper/manuscript.md "3.59\\times10^{-11}"
dchk paper/manuscript.md "1.39\\times10^{23}"
dchk paper/manuscript.md "118.3"
dchk paper/SI.md "4.801095"
# ★ Fig.4 的四个口径数：今天就因为"三列均值"与"全局极值"混用，在报告里写错过两个数。
#   把它们钉在正文里 —— 重绘 Fig.4 若得出别的数，这里先红。
dchk paper/manuscript.md "median 46.5 mV dec"
dchk paper/manuscript.md "49.2% of those points"
dchk paper/manuscript.md "spans 0.5000 to 0.9999"
dchk paper/manuscript.md "4992 solutions in total"
# 断言总数也不能漂：加了测试就得回去改 SI 和正文里的这几句话
dchk paper/SI.md "82 assertions in total"
dchk paper/manuscript.md "82 assertions"
# ★ 三档配置的数字必须跟着代码走。
#   （2026-10-23 为 Fig.4 的四个口径数加了 4 项 ⇒ 76/62/61 → 80/66/65。
#     2026-10-24 参考值随包分发 + 新增 T17 场效应 ⇒ 不再有"第三方档"，
#     三档变成 81 / 81 / 76（最后一档的 5 项跳过全部是可选件缺失）。）
dchk paper/manuscript.md "77 of the 82 pass"
dchk paper/SI.md  "82 passed / 0 failed / 0 skipped"
dchk paper/SI.md  "77 passed / 0 failed / 5 skipped"

echo
echo "══ ⑤ 错误路径：报错要对，退出码要对 ══"
e(){ printf '%s' "$2" > bad.kinstep.toml; "$KS" check bad.kinstep.toml >"$TMP/e.log" 2>&1; rc=$?
     if [ "$rc" = "$3" ]; then ok "$1" "exit=$rc"; else no "$1" "exit=$rc 期望 $3"; fi; }
e "机理不配平" '[mechanism]
text = """
CO + * <-> CO2*
"""
[conditions]
T_K = 500.0
[[steps]]
Ga_fwd = 1.0
dG_rxn = -0.5
' 2
e "少给一步能量" '[mechanism]
text = """
CO + * <-> CO*
O2 + 2* <-> 2O*
"""
[conditions]
T_K = 500.0
[[steps]]
Ga_fwd = 1.0
dG_rxn = -0.5
' 2
e "缺温度 T_K" '[mechanism]
text = """
CO + * <-> CO*
"""
[conditions]
p_gas = { CO = 1.0 }
[[steps]]
Ga_fwd = 1.0
dG_rxn = -0.5
' 2
e "能量来源不可信" '[mechanism]
text = """
CO + * <-> CO*
"""
[conditions]
T_K = 500.0
[[steps]]
Ga_fwd = 0.8
dG_rxn = -0.5
source = "llm_suggested"
' 3
"$KS" check 不存在.toml >/dev/null 2>&1; [ $? = 2 ] && ok "文件不存在" "exit=2" || no "文件不存在"

echo
echo "════════════════════════════════════════"
printf "  通过 %d 项，失败 %d 项" "$pass" "$fail"
# 跳过有两种原因：缺第三方参考数据（CATKINAS），或缺可选依赖（matplotlib）。
# 标签不能只写前一种——只有 numpy+scipy 时那 15 项里有 1 项是后者。
[ "$skipped" = 0 ] && printf "\n" || printf "，跳过 %d 项（缺第三方参考数据或可选依赖）\n" "$skipped"
[ "$fail" = 0 ] && echo "  ✓ 全部通过" || echo "  ✗ 有失败，失败日志已复制到 notes/FAILED_*.log"
exit $([ "$fail" = 0 ] && echo 0 || echo 1)
