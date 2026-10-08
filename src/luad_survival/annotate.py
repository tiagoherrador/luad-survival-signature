"""Optional gene-symbol annotation (Ensembl ID -> HGNC symbol) through the public MyGene.info service."""
import re
import warnings


def strip_version(ensembl_id: str) -> str:
    """ENSG00000123456.12 -> ENSG00000123456"""
    return re.sub(r"\.\d+$", "", str(ensembl_id))


def _fetch_mygene(ids: list[str]) -> list[dict]:
    import requests
    r = requests.post("https://mygene.info/v3/query", timeout=30,
                      data={"q": ",".join(ids), "scopes": "ensembl.gene", "fields": "symbol,name", "species": "human"})
    r.raise_for_status()
    return r.json()


def ensembl_to_symbol(ids, fetch=_fetch_mygene, batch: int = 100) -> dict[str, dict]:
    """Returns {ensembl id as given: {"symbol": ..., "name": ...}}. If the lookup fails, returns {} with a warning
    (the analysis never depends on it). Only gene identifiers are sent to the service, never patient data."""
    ids = list(ids)
    base = {i: strip_version(i) for i in ids}
    out: dict[str, dict] = {}
    try:
        found: dict[str, dict] = {}
        uniq = sorted(set(base.values()))
        for k in range(0, len(uniq), batch):
            for hit in fetch(uniq[k:k + batch]):
                if "symbol" in hit and hit.get("query") not in found:
                    found[hit["query"]] = {"symbol": hit["symbol"], "name": hit.get("name", "")}
        for i, b in base.items():
            if b in found:
                out[i] = found[b]
    except Exception as e:  # network down, service change, ...
        warnings.warn(f"gene-symbol lookup skipped ({e}); Ensembl IDs are reported instead")
        return {}
    return out
