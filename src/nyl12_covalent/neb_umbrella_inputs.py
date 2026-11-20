"""Build old-band DFTB3/MM MD inputs with an exact frozen outer boundary.

MM waters are deliberately flexible (ntc=1, ntf=1), matching the retained
umbrella convention. qmshake=0 explicitly leaves QM protons unconstrained.
"""

from __future__ import annotations

from pathlib import Path
import re

DT_PS = 0.0005
STAGE_STEPS = {"heat": 2000, "equil": 10000, "production": 10000}


def _block(text: str, name: str) -> str:
    matches = re.findall(
        rf"&{name}\b.*?(?:\n\s*/|&end\b)", text, flags=re.IGNORECASE | re.DOTALL
    )
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one complete &{name} block")
    return matches[0]


def _value(block: str, key: str) -> str:
    found = re.findall(
        rf"\b{key}\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s/]+)", block, flags=re.IGNORECASE
    )
    if len(found) != 1:
        raise ValueError(f"Expected one explicit {key}")
    return found[0].strip("'\"")


def _number(value: str) -> float:
    return float(value.lower().replace("d", "e"))


def audited_qmmm_block(old_mdin: str) -> str:
    """Retain the audited old Hamiltonian, making qmshake=0 explicit.

    Reject another charge, method or electronic temperature. The 98-real-atom
    mask must be numeric and explicit; this function verifies its cardinality.
    """
    block = _block(old_mdin, "qmmm")
    if (
        _value(block, "qm_theory").upper() != "DFTB3"
        or _number(_value(block, "qmcharge")) != -2
        or _number(_value(block, "dftb_telec")) != 0
    ):
        raise ValueError("Expected audited DFTB3/charge -2/electronic 0 K")
    mask = _value(block, "qmmask")
    if not re.fullmatch(r"@[0-9,\-]+", mask):
        raise ValueError("Need explicit numeric old QM atom mask")
    atoms = []
    for token in mask[1:].split(","):
        ends = [int(x) for x in token.split("-")]
        if len(ends) == 1:
            atoms.extend(ends)
        elif len(ends) == 2 and ends[0] <= ends[1]:
            atoms.extend(range(ends[0], ends[1] + 1))
        else:
            raise ValueError("Invalid QM mask range")
    if len(atoms) != 98 or len(set(atoms)) != 98 or min(atoms) <= 0:
        raise ValueError("Expected 98 distinct one-based real QM atoms")
    if re.search(r"\bqmshake\s*=", block, flags=re.IGNORECASE):
        if _number(_value(block, "qmshake")) != 0:
            raise ValueError("QM proton constraints are prohibited")
        return block
    # The old minimization omitted qmshake; this is the declared MD override.
    end = block.rfind("/")
    if end < 0:
        raise ValueError("Slash-terminated audited qmmm block required")
    return block[:end] + " qmshake=0,\n" + block[end:]


def active_group_cards(active_mask: str) -> str:
    """Encode the exact numeric old selection without 256-character truncation.

    Native mdread2 calls rgroup when bellymask is empty. Each RES/ATOM card
    contains one inclusive range; END closes the group and its group list.
    Cards immediately follow qmmm, before wt namelists: the QM reader leaves
    the input stream there before mdread2 calls rgroup.
    """
    if not re.fullmatch(r"[:@][0-9,\-]+", active_mask):
        raise ValueError("Group boundary requires numeric residue or atom mask")
    kind = "RES" if active_mask.startswith(":") else "ATOM"
    cards = ["Exact old-band movable selection"]
    for token in active_mask[1:].split(","):
        pair = [int(v) for v in token.split("-")]
        if len(pair) == 1:
            pair *= 2
        if len(pair) != 2 or not 0 < pair[0] <= pair[1]:
            raise ValueError("Invalid active selection range")
        cards.append(f"{kind} {pair[0]} {pair[1]}")
    return "\n".join([*cards, "END", "END"]) + "\n"


def build_mdin(
    old_mdin: str, stage: str, seed: int, active_mask: str, restraint_path: Path
) -> str:
    """Build a native NVT umbrella stage without copying DL-FIND controls.

    Args:
        old_mdin: Audited old-band input including qmmm and dlfind blocks.
        stage: heat (1 ps), equil (5 ps), or production (5 ps/block).
        seed: Explicit independent positive signed-32-bit Langevin seed.
        active_mask: Exact old DL-FIND active selection, encoded as legacy belly group cards.
        restraint_path: Restraint file path without whitespace or newlines.

    Returns:
        Native Sander input using the old Hamiltonian and frozen exterior.
    """
    if stage not in STAGE_STEPS:
        raise ValueError("Unsupported stage")
    if type(seed) is not int or not 0 < seed <= 2147483647:
        raise ValueError("Seed must be a positive signed-32-bit integer")
    active = _value(_block(old_mdin, "dlfind"), "active")
    if (
        active_mask != active
        or not active_mask
        or any(c in active_mask for c in "'\"\n")
    ):
        raise ValueError("Active mask must exactly match audited old-band active")
    path = str(restraint_path)
    if not path or any(c.isspace() for c in path) or any(c in path for c in "'\""):
        raise ValueError("Unsafe restraint path")
    groups = active_group_cards(active_mask)
    qmmm = audited_qmmm_block(old_mdin)
    steps = STAGE_STEPS[stage]
    start = "irest=0, ntx=1, tempi=25.0" if stage == "heat" else "irest=1, ntx=5"
    ramp = (
        f"&wt type='TEMP0', istep1=0, istep2={steps}, value1=25.0, value2=300.0 /\n"
        if stage == "heat"
        else ""
    )
    return f"""Old-band two-CV umbrella {stage}; flexible MM waters, frozen exterior
&cntrl
 imin=0, {start}, nstlim={steps}, dt={DT_PS:.7f},
 ntb=1, ntp=0, igb=0, cut=10.0, iwrap=0,
 ntt=3, gamma_ln=2.0, temp0=300.0, ig={seed},
 ntc=1, ntf=1, ifqnt=1, jfastw=4, ntr=0, nmropt=1,
 ibelly=1, bellymask='',
 ntpr=100, ntwx=100, ntwr=500, ioutfm=1, ntxo=2,
/
{qmmm}
{groups}{ramp}&wt type='DUMPFREQ', istep1=10 /
&wt type='END' /
DISANG={path}
DUMPAVE=umbrella_coordinate.dat
"""
