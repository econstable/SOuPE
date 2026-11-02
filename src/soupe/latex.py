"""Export the implemented Stevens expressions by inspecting Python syntax.

This renders the formula kernel for review; it is not an independent
derivation or a general symbolic algebra engine.
"""

from __future__ import annotations

import ast
import inspect
import re
from dataclasses import dataclass
from pathlib import Path

from . import _stevens as stevens_module

__all__ = ["render_latex", "write_latex"]


_OP_NAME_RE = re.compile(r"^O_(\d+)_(\d+)(?:_([cs]))?$")
_OP_TOKEN_RE = re.compile(r"^O_?(\d+)_?(\d+)(?:_?([csCS]))?$")
_OPS_ATTR_MAP = {
    "Jx": "J_x",
    "Jy": "J_y",
    "Jz": "J_z",
    "Jp": "J_+",
    "Jm": "J_-",
    "I": "I",
    "JJ": "J(J+1)",
}
_COMP_ORDER = {None: 0, "c": 1, "s": 2}


@dataclass
class OperatorDef:
    name: str
    k: int
    q: int
    component: str | None
    assignments: list[tuple[str, ast.AST]]
    ret: ast.AST


def _number_to_latex(value) -> str:
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:g}"
    if isinstance(value, complex):
        re_part = float(value.real)
        im_part = float(value.imag)
        if abs(re_part) < 1e-15:
            return f"{im_part:g} i"
        if abs(im_part) < 1e-15:
            return f"{re_part:g}"
        sign = "+" if im_part >= 0 else "-"
        return f"({re_part:g} {sign} {abs(im_part):g} i)"
    return str(value)


def _name_to_latex(name: str) -> str:
    if len(name) == 1 and name.isalpha():
        return name
    return f"\\mathrm{{{name}}}"


def _expr_to_latex(
    node: ast.AST,
    *,
    env: dict[str, ast.AST] | None = None,
    _seen: set[str] | None = None,
    preserve_names: set[str] | None = None,
) -> str:
    if isinstance(node, ast.Constant):
        return _number_to_latex(node.value)

    if isinstance(node, ast.Name):
        if preserve_names is not None and node.id in preserve_names:
            return _name_to_latex(node.id)
        if env is not None and node.id in env:
            seen = set() if _seen is None else set(_seen)
            if node.id in seen:
                return _name_to_latex(node.id)
            seen.add(node.id)
            return _expr_to_latex(
                env[node.id],
                env=env,
                _seen=seen,
                preserve_names=preserve_names,
            )
        return _name_to_latex(node.id)

    if isinstance(node, ast.Attribute):
        if isinstance(node.value, ast.Name) and node.value.id == "ops":
            mapped = _OPS_ATTR_MAP.get(node.attr)
            if mapped is not None:
                return mapped
            return f"\\mathrm{{ops.{node.attr}}}"
        return f"\\mathrm{{{ast.unparse(node)}}}"

    if isinstance(node, ast.UnaryOp):
        operand = _expr_to_latex(
            node.operand,
            env=env,
            _seen=_seen,
            preserve_names=preserve_names,
        )
        if isinstance(node.op, ast.USub):
            return f"-({operand})"
        return f"({operand})"

    if isinstance(node, ast.BinOp):
        lhs = _expr_to_latex(
            node.left, env=env, _seen=_seen, preserve_names=preserve_names
        )
        rhs = _expr_to_latex(
            node.right, env=env, _seen=_seen, preserve_names=preserve_names
        )
        if isinstance(node.op, ast.Add):
            return f"({lhs} + {rhs})"
        if isinstance(node.op, ast.Sub):
            return f"({lhs} - {rhs})"
        if isinstance(node.op, ast.Mult):
            return f"({lhs} {rhs})"
        if isinstance(node.op, ast.Div):
            return f"\\frac{{{lhs}}}{{{rhs}}}"
        if isinstance(node.op, ast.Pow):
            return f"({lhs})^{{{rhs}}}"
        return f"\\mathrm{{{ast.unparse(node)}}}"

    if isinstance(node, ast.Call):
        func_name = (
            node.func.id
            if isinstance(node.func, ast.Name)
            else ast.unparse(node.func)
        )
        args = [
            _expr_to_latex(
                a, env=env, _seen=_seen, preserve_names=preserve_names
            )
            for a in node.args
        ]
        if func_name == "_sym" and len(args) == 2:
            a, b = args
            return f"\\frac{{1}}{{2}}({a} {b} + {b} {a})"
        if func_name == "float" and len(args) == 1:
            return args[0]
        joined = ", ".join(args)
        return f"\\mathrm{{{func_name}}}({joined})"

    return f"\\mathrm{{{ast.unparse(node)}}}"


def _op_label_latex(op: OperatorDef) -> str:
    base = f"O_{{{op.k}}}^{{{op.q}}}"
    if op.component is None:
        return base
    return f"{base}(\\mathrm{{{op.component}}})"


def _split_rhs(rhs: str, max_len: int = 95) -> list[str]:
    s = " ".join(rhs.split())
    if len(s) <= max_len:
        return [s]

    breaks: list[tuple[int, int]] = []
    brace_depth = 0
    paren_depth = 0
    i = 0
    while i < len(s):
        ch = s[i]
        if ch == "{":
            brace_depth += 1
        elif ch == "}":
            brace_depth = max(0, brace_depth - 1)
        elif ch == "(":
            paren_depth += 1
        elif ch == ")":
            paren_depth = max(0, paren_depth - 1)
        elif (
            brace_depth == 0
            and i >= 1
            and i + 1 < len(s)
            and s[i] in "+-"
            and s[i - 1] == " "
            and s[i + 1] == " "
        ):
            breaks.append((i, paren_depth))
        i += 1

    chunks: list[str] = []
    start = 0
    min_chunk = 0.55
    search_right = 1.35
    while len(s) - start > max_len:
        limit = start + max_len
        min_pos = start + int(max_len * min_chunk)
        candidates = [
            (pos, depth) for pos, depth in breaks if min_pos <= pos < limit
        ]
        if candidates:
            min_depth = min(depth for _, depth in candidates)
            shallow = [pos for pos, depth in candidates if depth == min_depth]
            split_at = shallow[-1]
        else:
            right_limit = min(len(s), start + int(max_len * search_right))
            right_candidates = [
                (pos, depth)
                for pos, depth in breaks
                if limit <= pos <= right_limit
            ]
            if right_candidates:
                min_depth = min(depth for _, depth in right_candidates)
                shallow = [
                    pos
                    for pos, depth in right_candidates
                    if depth == min_depth
                ]
                split_at = shallow[0]
            else:
                space = s.rfind(
                    " ",
                    start + int(max_len * 0.6),
                    min(len(s), start + int(max_len * 1.2)),
                )
                split_at = (
                    space if space > start else min(len(s), start + max_len)
                )
        chunks.append(s[start:split_at].strip())
        start = split_at
    tail = s[start:].strip()
    if tail:
        chunks.append(tail)
    return [c for c in chunks if c]


def _collect_operator_defs() -> dict[str, OperatorDef]:
    src = inspect.getsource(stevens_module)
    tree = ast.parse(src)

    out: dict[str, OperatorDef] = {}
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        m = _OP_NAME_RE.match(node.name)
        if m is None:
            continue
        k = int(m.group(1))
        q = int(m.group(2))
        comp = m.group(3)
        assignments: list[tuple[str, ast.AST]] = []
        ret = None
        for stmt in node.body:
            if (
                isinstance(stmt, ast.Expr)
                and isinstance(stmt.value, ast.Constant)
                and isinstance(stmt.value.value, str)
            ):
                continue
            if isinstance(stmt, ast.Assign):
                if len(stmt.targets) != 1 or not isinstance(
                    stmt.targets[0], ast.Name
                ):
                    raise ValueError(
                        f"Unsupported assignment form in {node.name}"
                    )
                assignments.append((stmt.targets[0].id, stmt.value))
                continue
            if isinstance(stmt, ast.Return):
                ret = stmt.value
                continue
            raise ValueError(
                f"Unsupported statement in {node.name}: {ast.unparse(stmt)}"
            )

        if ret is None:
            raise ValueError(f"No return expression found in {node.name}")
        out[node.name] = OperatorDef(
            name=node.name,
            k=k,
            q=q,
            component=comp,
            assignments=assignments,
            ret=ret,
        )
    return out


def _sorted_ops(all_ops: dict[str, OperatorDef]) -> list[str]:
    return sorted(
        all_ops.keys(),
        key=lambda name: (
            all_ops[name].k,
            all_ops[name].q,
            _COMP_ORDER[all_ops[name].component],
        ),
    )


def _split_ops_argument(ops_arg: str) -> list[str]:
    raw = ops_arg.strip()
    if raw.startswith("[") and raw.endswith("]"):
        raw = raw[1:-1]
    if "," in raw:
        tokens = [t.strip() for t in raw.split(",")]
    else:
        tokens = [t.strip() for t in raw.split()]
    return [t for t in tokens if t]


def _expand_token(token: str, all_ops: dict[str, OperatorDef]) -> list[str]:
    if token.upper() == "ALL":
        return _sorted_ops(all_ops)

    m = _OP_TOKEN_RE.match(token)
    if m is None:
        raise ValueError(f"Invalid operator token '{token}'")

    k = int(m.group(1))
    q = int(m.group(2))
    comp = m.group(3).lower() if m.group(3) else None

    if comp is not None:
        name = f"O_{k}_{q}_{comp}"
        if name not in all_ops:
            raise ValueError(
                f"Operator '{name}' is not defined in soupe._stevens"
            )
        return [name]

    no_comp = f"O_{k}_{q}"
    if no_comp in all_ops:
        return [no_comp]

    with_comp = [f"O_{k}_{q}_c", f"O_{k}_{q}_s"]
    present = [name for name in with_comp if name in all_ops]
    if present:
        return present

    raise ValueError(f"Operator 'O_{k}_{q}' is not defined in soupe._stevens")


def _resolve_requested_ops(
    ops_arg: str, all_ops: dict[str, OperatorDef]
) -> list[str]:
    if not ops_arg or not ops_arg.strip():
        raise ValueError("Empty --ops value")
    requested: list[str] = []
    seen = set()
    for tok in _split_ops_argument(ops_arg):
        for name in _expand_token(tok, all_ops):
            if name in seen:
                continue
            seen.add(name)
            requested.append(name)
    if not requested:
        raise ValueError("No operators resolved from --ops")
    return sorted(
        requested,
        key=lambda name: (
            all_ops[name].k,
            all_ops[name].q,
            _COMP_ORDER[all_ops[name].component],
        ),
    )


def _build_tex_doc(ops: list[OperatorDef]) -> str:
    lines = [
        r"\documentclass[11pt,fleqn]{article}",
        r"\usepackage{amsmath}",
        r"\usepackage[a4paper,margin=2cm]{geometry}",
        r"\begin{document}",
        r"\section*{Stevens Operator Expressions}",
        r"Generated directly from \texttt{soupe._stevens}.",
        r"\setlength{\mathindent}{0pt}",
        r"\setlength{\parindent}{0pt}",
        "",
        r"\begin{align*}",
        r"X &\equiv J(J+1)\,I",
        r"\end{align*}",
        "",
    ]
    for op in ops:
        label = _op_label_latex(op)
        env_map = {k: v for k, v in op.assignments}
        rhs = _expr_to_latex(op.ret, env=env_map, preserve_names={"X"})
        rhs_lines = _split_rhs(rhs, max_len=95)
        lines.append(r"\begin{align*}")
        if len(rhs_lines) == 1:
            lines.append(rf"{label} &= {rhs_lines[0]}")
        else:
            lines.append(rf"{label} &= {rhs_lines[0]} \\")
            for idx, frag in enumerate(rhs_lines[1:], start=1):
                end = r" \\" if idx < len(rhs_lines) - 1 else ""
                lines.append(rf"&\quad {frag}{end}")
        lines.append(r"\end{align*}")
        lines.append("")
    lines.append(r"\end{document}")
    return "\n".join(lines) + "\n"


def render_latex(ops_arg: str = "ALL") -> str:
    """Render the coded formulas as a complete LaTeX document.

    Parameters
    ----------
    ops_arg : str, optional
        'ALL' or comma/space-separated names, e.g. '[O_2_0,O_4_3_c]'.
        Omitting c/s for a nonzero order selects both components.

    Returns
    -------
    str
        Document derived from the implemented expressions. This is an
        inspection aid, not an independent reference for the formulas.
    """
    all_ops = _collect_operator_defs()
    selected_names = _resolve_requested_ops(ops_arg, all_ops)
    selected = [all_ops[name] for name in selected_names]
    return _build_tex_doc(selected)


def write_latex(path: Path | str, ops_arg: str = "ALL") -> Path:
    """Write selected Stevens operator definitions to a .tex file.

    Parameters
    ----------
    path : pathlib.Path or str
        Output path. The parent directory must exist.
    ops_arg : str, optional
        Selection accepted by `render_latex`.

    Returns
    -------
    pathlib.Path
        Path to the written UTF-8 document.
    """
    out = Path(path)
    out.write_text(render_latex(ops_arg), encoding="utf-8")
    return out
