AGONIST_ALIASES = {
    "compound 48/80",
    "compound 48-80",
    "compound48/80",
    "compound48-80",
    "substance p",
    "ciprofloxacin",
}

ANTAGONIST_ALIASES = {
    "qwf",
}

AGONIST_CIDS = {
    2855,  # Compound 48/80
    36511,  # Substance P
    2764,  # Ciprofloxacin
}

ANTAGONIST_CIDS = set[int]()

# PubChem's "Title" is often a long IUPAC string (e.g. CID 2855 resolves to
# "1-(2-Methylaminoethyl)-4-methoxy-3,5-bis(...)benzene" rather than the
# pharmacology name researchers actually use). Override with the familiar
# name wherever we know it, so the UI shows what the structure means in
# practice rather than its systematic chemical name.
CURATED_DISPLAY_NAMES: dict[int, str] = {
    2855: "Compound 48/80",
    36511: "Substance P",
    2764: "Ciprofloxacin",
}


def display_name_for_cid(cid: int | None, fallback: str) -> str:
    if cid is not None and cid in CURATED_DISPLAY_NAMES:
        return CURATED_DISPLAY_NAMES[cid]
    return fallback


def normalize_name(value: str) -> str:
    return " ".join(value.strip().lower().replace("_", " ").split())


def classify_curated_name(*names: str, cid: int | None = None) -> str | None:
    if cid in AGONIST_CIDS:
        return "agonist"
    if cid in ANTAGONIST_CIDS:
        return "antagonist"

    normalized_names = {normalize_name(name) for name in names if name.strip()}

    for name in normalized_names:
        for alias in AGONIST_ALIASES:
            if alias in name or name in alias:
                return "agonist"
        for alias in ANTAGONIST_ALIASES:
            if alias in name or name in alias:
                return "antagonist"

    return None


def resolve_curated_cid(name: str) -> int | None:
    normalized = normalize_name(name)

    for alias in AGONIST_ALIASES:
        if alias in normalized or normalized in alias:
            if "compound 48" in alias:
                return 2855
            if alias == "substance p":
                return 36511
            if alias == "ciprofloxacin":
                return 2764

    for alias in ANTAGONIST_ALIASES:
        if alias in normalized or normalized in alias:
            return None

    return None
