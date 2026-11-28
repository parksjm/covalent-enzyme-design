"""Honor explicit native TemplateInfo ranges before native template feature math.

The installed Boltz featurizer uses record offsets but omits range clipping.
This isolated input adapter changes row visibility/mapping only. It delegates
all coordinate, residue, frame and mask feature calculations to native code.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Sequence
from typing import Any


def bounded_template_rows(
    data: Any, template_name: str, records: Sequence[Any], template_id: int
) -> list[dict[str, Any]]:
    """Map only explicitly aligned source tokens into corresponding query tokens.

    Args:
        data: Native Boltz Tokenized input, including unchanged template tokens.
        template_name: Native assembly/template name.
        records: Native TemplateInfo records with explicit half-open ranges.
        template_id: Native template grouping index used for visibility.

    Returns:
        Native compute_template_features row dictionaries without new coordinates.

    Raises:
        ValueError: Records are forced, malformed, overlapping or incomplete.
        KeyError: An explicitly named query/template chain is missing.
    """
    query_chains = {str(c["name"]): int(c["asym_id"]) for c in data.structure.chains}
    source = data.templates[template_name]
    source_chains = {str(c["name"]): int(c["asym_id"]) for c in source.chains}
    rows: list[dict[str, Any]] = []
    seen: set[int] = set()
    for record in records:
        if record.name != template_name or record.force:
            raise ValueError("Expected matching, nonforced template records")
        if (
            record.query_chain not in "ABCDEFGH"
            or record.template_chain not in "ABCDEFGH"
        ):
            raise ValueError("This protein-only adapter excludes ligand chains")
        qstart, qend = record.query_st, record.query_en
        tstart, tend = record.template_st, record.template_en
        if min(qstart, tstart) < 0 or qend <= qstart or tend <= tstart:
            raise ValueError("Invalid half-open template alignment range")
        if qend - qstart != tend - tstart:
            raise ValueError("Explicit segments must be ungapped and equal-length")
        qchain = query_chains[record.query_chain]
        tchain = source_chains[record.template_chain]
        query_index = {
            int(t["res_idx"]): int(t["token_idx"])
            for t in data.tokens
            if int(t["asym_id"]) == qchain
        }
        offset = tstart - qstart
        mapped = 0
        for token in data.template_tokens[template_name]:
            residue = int(token["res_idx"])
            if int(token["asym_id"]) != tchain or not tstart <= residue < tend:
                continue
            query_residue = residue - offset
            if not qstart <= query_residue < qend or query_residue not in query_index:
                raise ValueError("Explicit source/query correspondence is absent")
            query_token = query_index[query_residue]
            if query_token in seen:
                raise ValueError("Overlapping template records map a query token twice")
            seen.add(query_token)
            rows.append({"token": token, "pdb_id": template_id, "q_idx": query_token})
            mapped += 1
        if mapped != tend - tstart:
            raise ValueError(
                "Explicit alignment segment lacks observed template tokens"
            )
    if not rows:
        raise ValueError("No explicit protein template rows")
    return rows


def make_bounded_template_featurizer(
    native_compute: Callable[..., dict[str, Any]], torch_module: Any
) -> Callable[[Any, int], dict[str, Any]]:
    """Build a process-local adapter using the installed native feature calculator.

    Args:
        native_compute: Installed compute_template_features implementation.
        torch_module: Already installed engine Torch module; no new dependency.

    Returns:
        Drop-in process_template_features callable for an isolated Boltz entry.
    """

    def process(data: Any, max_tokens: int) -> dict[str, Any]:
        grouped: dict[str, list[Any]] = defaultdict(list)
        for record in data.record.templates:
            grouped[record.name].append(record)
        if not grouped:
            raise ValueError("Explicit template records are required")
        features = []
        for template_id, (name, records) in enumerate(grouped.items()):
            rows = bounded_template_rows(data, name, records, template_id)
            result = native_compute(data, rows, max_tokens)
            result["template_force"] = torch_module.tensor(False)
            result["template_force_threshold"] = torch_module.tensor(
                float("inf"), dtype=torch_module.float32
            )
            features.append(result)
        return {
            key: torch_module.stack([f[key] for f in features]) for key in features[0]
        }

    return process
