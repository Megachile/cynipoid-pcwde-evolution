"""Source-qualified parent identity shared by Table S3 and prose counts."""
def origin_parent(header: str) -> tuple[str, str]:
    """Source-qualified protein ID from the three accepted origin-header fields."""
    parts = header.split("|")
    if parts[0] == "cynipoid" and len(parts) == 2:
        return (parts[0], parts[1])
    if parts[0].endswith("_ENA") and len(parts) == 3:
        return (parts[0], parts[2])  # field 2 is a contig, field 3 the protein
    if len(parts) == 3 and parts[0] in {"bacterium", "bee", "midge", "fungus", "beetle", "plant", "metazoa", "rotifer", "nematode", "protist", "oomycete"}:
        return (parts[0], parts[1])
    raise SystemExit(f"Unresolved origin parent protein: {header}")
