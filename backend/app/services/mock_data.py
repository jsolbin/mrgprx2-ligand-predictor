MOCK_COMPOUNDS = {
    "ciprofloxacin": {
        "name": "Ciprofloxacin",
        "cid": 2764,
        "smiles": "CC1COCCN1c2nc(N3CCNCC3)c(F)cc2C(=O)O",
        "inchi_key": "MYSWGUAQZAJSOK-UHFFFAOYSA-N",
        "molecular_weight": 331.34,
        "logp": 1.54,
        "tpsa": 74.57,
    },
    "morphine": {
        "name": "Morphine",
        "cid": 5288826,
        "smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
        "inchi_key": "BQJCRHHNABKAKU-KBQPJGBKSA-N",
        "molecular_weight": 285.34,
        "logp": 1.05,
        "tpsa": 52.93,
    },
    "vancomycin": {
        "name": "Vancomycin",
        "cid": 14969,
        "smiles": "C[C@H](NC(=O)[C@H](N)CO)C(=O)N[C@@H]1C(O)Oc2cc(O)cc(O)c2[C@H](O)[C@H](NC(=O)C(N)C(O)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]2C(O)Oc3cc(O)c(O)c([C@H](NC(=O)[C@H](N)CO)C(=O)N[C@@H](CO)C(N)=O)c3[C@H](O)[C@H](NC(=O)C(N)C(O)=O)C(=O)N[C@@H](CO)C(=O)N[C@@H](C(O)=O)c3cc(O)c(O)c(O)c3-c3c(O)cc(O)cc3O[C@H]1C2",
        "inchi_key": "XQKFTBPVXVTOEG-NQFVTFAPSA-N",
        "molecular_weight": 1449.25,
        "logp": -3.10,
        "tpsa": 576.76,
    },
    "substance p": {
        "name": "Substance P",
        "cid": 36511,
        "smiles": "NCC[C@H](N)C(=O)N[C@@H](Cc1ccccc1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N[C@@H](CCCNC(N)=N)C(=O)N[C@@H](CO)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCSC)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(N)=O",
        "inchi_key": "BDSIPYCVIZGJAV-QXLPPPPWSA-N",
        "molecular_weight": 1347.63,
        "logp": -1.80,
        "tpsa": 487.20,
    },
}


def get_demo_svg(label: str) -> str:
    safe_label = label[:32]
    return f"""
<svg width="200" height="200" viewBox="0 0 200 200" fill="none" xmlns="http://www.w3.org/2000/svg">
  <rect width="200" height="200" rx="24" fill="#FFF8EA"/>
  <circle cx="70" cy="100" r="18" fill="#1D6B52"/>
  <circle cx="130" cy="70" r="14" fill="#C9703D"/>
  <circle cx="130" cy="130" r="14" fill="#C9703D"/>
  <path d="M84 93L116 77" stroke="#10231C" stroke-width="6" stroke-linecap="round"/>
  <path d="M84 107L116 123" stroke="#10231C" stroke-width="6" stroke-linecap="round"/>
  <path d="M130 84V116" stroke="#10231C" stroke-width="6" stroke-linecap="round"/>
  <text x="100" y="180" text-anchor="middle" fill="#10231C" font-size="12" font-family="monospace">{safe_label}</text>
</svg>
""".strip()
