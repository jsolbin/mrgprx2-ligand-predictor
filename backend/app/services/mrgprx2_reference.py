"""Curated database of known MRGPRX2 agonists and antagonists.

This extends the small seed list in `predict_service.KNOWN_LIGANDS` with a
broader set of literature-reported ligands so that structure-based similarity
matching (Tanimoto over Morgan fingerprints) has more reference points to
learn from. Each entry needs a resolvable small-molecule SMILES to participate
in fingerprint comparison - large peptides/biologics whose PubChem SMILES are
unwieldy for substructure-style reasoning are still listed for name lookups,
but are flagged with `structural_reference: False` so the similarity search
can skip them gracefully.

To teach the classifier about a newly discovered ligand, append an entry here
(or submit it through the `/training/label` endpoint, which is merged in at
prediction time via `training_store.reference_ligands_from_labels`).
"""

MRGPRX2_REFERENCE_LIGANDS: list[dict] = [
    # --- Agonists -----------------------------------------------------
    {
        "name": "Substance P",
        "label": "agonist",
        "smiles": (
            "CC(C)C[C@@H](C(=O)N[C@@H](CCSC)C(=O)N)NC(=O)CNC(=O)[C@H](CC1=CC=CC=C1)"
            "NC(=O)[C@H](CC2=CC=CC=C2)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CCC(=O)N)"
            "NC(=O)[C@@H]3CCCN3C(=O)[C@H](CCCCN)NC(=O)[C@@H]4CCCN4C(=O)"
            "[C@H](CCCN=C(N)N)N"
        ),
        "rationale": "Endogenous neuropeptide and prototypical MRGPRX2 agonist; C-terminal basic/aromatic residues drive direct Gq activation and mast-cell degranulation.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Cortistatin-14",
        "label": "agonist",
        "smiles": (
            "C[C@H]([C@H]1C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CSSC[C@@H]"
            "(C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H]"
            "(C(=O)N1)CCCCN)CC2=CNC3=CC=CC=C32)CC4=CC=CC=C4)CC5=CC=CC=C5)CC(=O)N)"
            "CCCCN)NC(=O)[C@@H]6CCCN6)C(=O)N[C@@H](CCCCN)C(=O)O)CO)CO)CC7=CC=CC=C7)O"
        ),
        "rationale": "Neuropeptide agonist sharing a cyclic, cationic, aromatic-rich pharmacophore with somatostatin-family peptides reported to engage MRGPRX2.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "PACAP 1-27",
        "label": "agonist",
        "smiles": (
            "CC[C@H](C)[C@@H](C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H]([C@@H](C)O)"
            "C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CO)C(=O)N[C@@H](CC2=CC=C(C=C2)O)"
            "C(=O)N[C@@H](CO)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CC3=CC=C(C=C3)O)"
            "C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCC(=O)N)"
            "C(=O)N[C@@H](CCSC)C(=O)N[C@@H](C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCCCN)"
            "C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC4=CC=C(C=C4)O)C(=O)N[C@@H](CC(C)C)"
            "C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(C)C)"
            "C(=O)N)NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@H]"
            "(CC5=CN=CN5)N"
        ),
        "rationale": "PACAP/VIP-family neuropeptide agonist; basic and aromatic side-chain density mirrors other MRGPRX2-activating neuropeptides (PAMP, Substance P).",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "LL-37",
        "label": "agonist",
        "smiles": (
            "CC[C@H](C)[C@@H](C(=O)NCC(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCC(=O)O)"
            "C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCNC(=N)N)"
            "C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(=O)N)"
            "C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](CCCCN)"
            "C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CC2=CC=CC=C2)C(=O)N[C@@H](CC(C)C)"
            "C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(C)C)"
            "C(=O)N[C@@H](C(C)C)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CCCNC(=N)N)C(=O)N"
            "[C@@H]([C@@H](C)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CO)C(=O)O)"
            "NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCCN)NC(=O)"
            "[C@H](CO)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H]"
            "(CC4=CC=CC=C4)NC(=O)[C@H](CC5=CC=CC=C5)NC(=O)[C@H](CC(=O)O)NC(=O)"
            "CNC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)N"
        ),
        "rationale": "Cathelicidin host-defense peptide; cationic amphipathic helix is a well-documented MRGPRX2 agonist motif (alongside human beta-defensins).",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "PAMP-12",
        "label": "agonist",
        "smiles": (
            "CC(C)CC(C(=O)NC(CO)C(=O)NC(CCCNC(=N)N)C(=O)N)NC(=O)C(C)NC(=O)C"
            "(CC1=CNC2=CC=CC=C21)NC(=O)C(CCCCN)NC(=O)C(CC(=O)N)NC(=O)C(N)Cc3ccccc3"
        ),
        "rationale": "C-terminal fragment of proadrenomedullin N-terminal peptide (PAMP); short, highly cationic/aromatic peptide consistent with the MRGPRX2 agonist pharmacophore.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Human beta-defensin-2",
        "label": "agonist",
        "smiles": (
            "CC[C@H](C)[C@H]1C(=O)N[C@H]2CSSC[C@H]3C(=O)NCC(=O)N[C@H](C(=O)N4CCC"
            "[C@H]4C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CSSC[C@@H](NC1=O)"
            "C(=O)N2)C(=O)O)CC(=O)O)Cc5ccc(O)cc5)Cc6c[nH]c7ccccc67"
        ),
        "rationale": "Cationic, disulfide-rich antimicrobial peptide; reported direct MRGPRX2 agonist alongside LL-37, contributing to neurogenic-inflammation signalling.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Atracurium",
        "label": "agonist",
        "smiles": (
            "C[N+]1(CCC2=CC(=C(C=C2C1CC3=CC(=C(C=C3)OC)OC)OC)OC)CCC(=O)OCCCCCOC(=O)"
            "CC[N+]4(CCC5=CC(=C(C=C5C4CC6=CC(=C(C=C6)OC)OC)OC)OC)C"
        ),
        "rationale": "Bis-quaternary benzylisoquinolinium neuromuscular blocker; canonical example of small-molecule drugs that trigger MRGPRX2-mediated pseudo-allergic reactions.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Rocuronium",
        "label": "agonist",
        "smiles": (
            "CC(=O)O[C@H]1[C@H](C[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(C"
            "[C@@H]([C@H](C4)O)N5CCOCC5)C)C)[N+]6(CCCC6)CC=C"
        ),
        "rationale": "Aminosteroid neuromuscular blocker with a quaternary ammonium head group; clinically associated with MRGPRX2-driven anaphylactoid reactions.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Cisatracurium",
        "label": "agonist",
        "smiles": (
            "C[N@@+]1(CCC2=CC(=C(C=C2[C@H]1CC3=CC(=C(C=C3)OC)OC)OC)OC)CCC(=O)OCCCCCOC(=O)"
            "CC[N@+]4(CCC5=CC(=C(C=C5[C@H]4CC6=CC(=C(C=C6)OC)OC)OC)OC)C"
        ),
        "rationale": "Stereoisomer of atracurium; shares the bis-cationic benzylisoquinolinium scaffold linked to direct MRGPRX2 activation and histamine release.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Morphine",
        "label": "agonist",
        "smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=C[C@H]3[C@H]1C5)O",
        "rationale": "Opioid agonist; morphine and related opioids directly activate MRGPRX2 on mast cells, a documented mechanism behind opioid-induced pseudo-allergic (pruritus/flushing) reactions independent of mu-opioid signalling.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Ciprofloxacin",
        "label": "agonist",
        "smiles": "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
        "rationale": "Fluoroquinolone antibiotic; the basic piperazine ring plus fused bicyclic aromatic core is a recurring MRGPRX2 agonist signature among 'pseudo-allergy' drugs.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Levofloxacin",
        "label": "agonist",
        "smiles": "C[C@H]1COC2=C3N1C=C(C(=O)C3=CC(=C2N4CCN(CC4)C)F)C(=O)O",
        "rationale": "Fluoroquinolone agonist; methylated piperazine plus tricyclic aromatic system mirrors ciprofloxacin's MRGPRX2-activating scaffold.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Moxifloxacin",
        "label": "agonist",
        "smiles": "COC1=C2C(=CC(=C1N3C[C@@H]4CCCN[C@@H]4C3)F)C(=O)C(=CN2C5CC5)C(=O)O",
        "rationale": "Fourth-generation fluoroquinolone agonist; bulkier diazabicyclic side chain still presents the basic-amine + aromatic-core pattern that engages MRGPRX2.",
        "source": "curated",
        "structural_reference": True,
    },
    # --- Antagonists ---------------------------------------------------
    {
        "name": "QWF",
        "label": "antagonist",
        "smiles": "CC(C)(C)OC(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc1ccccc1)C(N)=O",
        "rationale": "Tripeptide-derived MRGPRX2 antagonist reference; bulky diaromatic, low-cationic scaffold blocks the orthosteric pocket without triggering activation.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Quercetin",
        "label": "antagonist",
        "smiles": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
        "rationale": "Flavonol reported to inhibit MRGPRX2-mediated mast-cell degranulation; polyphenolic flavone core lacks the basic-amine agonist trigger.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Luteolin",
        "label": "antagonist",
        "smiles": "C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O",
        "rationale": "Flavone antagonist/inhibitor of MRGPRX2-driven degranulation; shares quercetin's polyphenolic flavone scaffold without a protonatable amine.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Kaempferol",
        "label": "antagonist",
        "smiles": "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O",
        "rationale": "Flavonol structurally close to quercetin/myricetin; part of the polyphenol family reported to dampen MRGPRX2 responses.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Myricetin",
        "label": "antagonist",
        "smiles": "C1=C(C=C(C(=C1O)O)O)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O",
        "rationale": "Highly hydroxylated flavonol; extends the quercetin/luteolin/kaempferol polyphenol antagonist series associated with reduced MRGPRX2 activation.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Osthole",
        "label": "antagonist",
        "smiles": "CC(=CCC1=C(C=CC2=C1OC(=O)C=C2)OC)C",
        "rationale": "Coumarin natural product reported to suppress MRGPRX2-mediated mast-cell activation; neutral lipophilic scaffold without a basic-amine agonist trigger.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Cromoglycate (DSCG)",
        "label": "antagonist",
        "smiles": "O=C(O)c1cc(=O)c2c(OCC(O)COc3c4OC(=O)c(c4cc(c3)C(=O)O)O)ccc2o1",
        "rationale": "Bis-chromone mast-cell stabiliser (cromolyn/disodium cromoglycate); clinically used to dampen mast-cell mediator release, consistent with an MRGPRX2-blocking/antagonist role rather than a direct degranulation trigger.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Boc-QWF-OBn",
        "label": "antagonist",
        "smiles": (
            "CC(C)(C)OC(=O)N[C@@H](CCC(N)=O)C(=O)N[C@H](CC1=CN(C=O)C2=CC=CC=C12)"
            "C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)OCC1=CC=CC=C1"
        ),
        "rationale": "Boc/benzyl-protected Gln-Trp-Phe tripeptide (CAS 126088-82-2); reported MRGPRX2 antagonist tool compound whose blocked N- and C-termini distinguish it from the free-terminus aromatic-basic peptide motifs that drive agonism.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Ketotifen",
        "label": "antagonist",
        "smiles": "CN1CCC(=C2C3=C(C(=O)CC4=CC=CC=C42)SC=C3)CC1",
        "rationale": "Antihistamine/mast-cell stabiliser; tricyclic benzocycloheptathiophene scaffold reported to block MRGPRX2-linked TM6 conformational movement rather than trigger it.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Nedocromil",
        "label": "antagonist",
        "smiles": "CCCC1=C2C(=CC3=C1OC(=CC3=O)C(=O)O)C(=O)C=C(N2CC)C(=O)O",
        "rationale": "Mast-cell stabiliser structurally related to cromolyn; anionic bis-carboxylate pyranoquinolinedione scaffold consistent with an MRGPRX2-blocking rather than activating role.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Silibinin",
        "label": "antagonist",
        "smiles": "COC1=C(C=CC(=C1)[C@@H]2[C@H](OC3=C(O2)C=C(C=C3)[C@@H]4[C@H](C(=O)C5=C(C=C(C=C5O4)O)O)O)CO)O",
        "rationale": "Flavonolignan natural product; polyphenolic scaffold with no protonatable cationic group, consistent with the quercetin/luteolin polyphenol-antagonist family that dampens MRGPRX2 activation.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Dexamethasone",
        "label": "antagonist",
        "smiles": "C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@@]4([C@]3([C@H](C[C@@]2([C@]1(C(=O)CO)O)C)O)F)C",
        "rationale": "Corticosteroid repositioned against mast-cell-driven disease; suppresses downstream MRGPRX2 signalling rather than engaging the orthosteric pocket as a trigger - neutral steroid scaffold lacks an agonist-type basic amine.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Azelastine",
        "label": "antagonist",
        "smiles": "CN1CCCC(CC1)N2C(=O)C3=CC=CC=C3C(=N2)CC4=CC=C(C=C4)Cl",
        "rationale": "Antihistamine repositioned as an MRGPRX2 blocker; bulky phthalazinone scaffold reported to cause steric clash in the orthosteric pocket that prevents the activating conformational shift.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Tubocurarine",
        "label": "agonist",
        "smiles": "CN1CCC2=CC(=C3C=C2[C@@H]1CC4=CC=C(C=C4)OC5=C6[C@@H](CC7=CC(=C(C=C7)O)O3)[N+](CCC6=CC(=C5O)OC)(C)C)OC",
        "rationale": "Bis-quaternary-ammonium neuromuscular blocker; classic cause of perioperative anaphylactoid reactions through direct MRGPRX2-mediated mast-cell degranulation.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Atracurium",
        "label": "agonist",
        "smiles": "C[N+]1(CCC2=CC(=C(C=C2C1CC3=CC(=C(C=C3)OC)OC)OC)OC)CCC(=O)OCCCCCOC(=O)CC[N+]4(CCC5=CC(=C(C=C5C4CC6=CC(=C(C=C6)OC)OC)OC)OC)C",
        "rationale": "Benzylisoquinolinium neuromuscular blocker with well-documented direct histamine-releasing/MRGPRX2-agonist activity, mirroring Compound 48/80's polyamine-methoxyaromatic pharmacophore.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Rocuronium",
        "label": "agonist",
        "smiles": "CC(=O)O[C@H]1[C@H](C[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(C[C@@H]([C@H](C4)O)N5CCOCC5)C)C)[N+]6(CCCC6)CC=C",
        "rationale": "Aminosteroid neuromuscular blocker bearing a quaternary nitrogen; reported MRGPRX2 agonist linked to perioperative pseudo-allergic reactions.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Moxifloxacin",
        "label": "agonist",
        "smiles": "COC1=C2C(=CC(=C1N3C[C@@H]4CCCN[C@@H]4C3)F)C(=O)C(=CN2C5CC5)C(=O)O",
        "rationale": "Fluoroquinolone antibiotic reported to activate MRGPRX2 more potently than ciprofloxacin; bicyclic basic-amine substituent drives direct mast-cell degranulation.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Levofloxacin",
        "label": "agonist",
        "smiles": "C[C@H]1COC2=C3N1C=C(C(=O)C3=CC(=C2N4CCN(CC4)C)F)C(=O)O",
        "rationale": "Fluoroquinolone bearing an N-methylpiperazine basic amine; documented MRGPRX2 agonist linked to fluoroquinolone-associated pseudo-allergic reactions.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Vancomycin",
        "label": "agonist",
        "smiles": "C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O",
        "rationale": "Glycopeptide antibiotic; the well-known cause of \"Red Man Syndrome\", mediated by direct MRGPRX2 activation and mast-cell histamine release rather than a true IgE allergy.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Icatibant",
        "label": "agonist",
        "smiles": "C1CC[C@H]2[C@@H](C1)C[C@H](N2C(=O)[C@H]3CC4=CC=CC=C4CN3C(=O)[C@H](CO)NC(=O)[C@H](CC5=CC=CS5)NC(=O)CNC(=O)[C@@H]6C[C@H](CN6C(=O)[C@@H]7CCCN7C(=O)[C@H](CCCN=C(N)N)NC(=O)[C@@H](CCCN=C(N)N)N)O)C(=O)N[C@@H](CCCN=C(N)N)C(=O)O",
        "rationale": "Peptidic bradykinin B2-receptor antagonist bearing multiple arginine residues (net charge ~+3); the dense cationic/aromatic-basic motif is reported to directly activate MRGPRX2, in line with other cationic neuropeptide agonists.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "PAMP-12",
        "label": "agonist",
        "smiles": "CC(C)CC(C(=O)NC(CO)C(=O)NC(CCCNC(=N)N)C(=O)N)NC(=O)C(C)NC(=O)C(CC1=CNC2=CC=CC=C21)NC(=O)C(CCCCN)NC(=O)C(CC(=O)N)NC(=O)C(CC3=CNC4=CC=CC=C43)NC(=O)C(CCCCN)NC(=O)C(CCCCN)NC(=O)C(CCCNC(=N)N)NC(=O)C(CC5=CC=CC=C5)N",
        "rationale": "C-terminal fragment of proadrenomedullin (PAMP); endogenous neuropeptide agonist with the cationic, aromatic-rich pharmacophore characteristic of direct MRGPRX2 activators such as Substance P.",
        "source": "curated",
        "structural_reference": True,
    },
    # --- Nonbinders (negative controls) -------------------------------
    # Structures reported to NOT engage MRGPRX2 - the classifier needs
    # these as a third class so it can learn what "doesn't bind" looks
    # like, rather than only ever choosing between agonist/antagonist.
    {
        "name": "D-Glucose",
        "label": "nonbinder",
        "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "rationale": "Neutral carbohydrate with no charge and very low LogP; lacks any cationic or aromatic features MRGPRX2 ligands rely on.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "D-Galactose",
        "label": "nonbinder",
        "smiles": "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O",
        "rationale": "Epimer of glucose used as a structurally-close negative control; same neutral, highly polar sugar scaffold with no MRGPRX2 activity.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Sucrose",
        "label": "nonbinder",
        "smiles": "OC[C@H]1O[C@H](O[C@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@@H](O)[C@@H]1O",
        "rationale": "Larger disaccharide - more rings and mass than many ligands, but uncharged and hydrophilic, so it doesn't engage the cationic MRGPRX2 pocket.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "L-Alanine",
        "label": "nonbinder",
        "smiles": "CC(C(=O)O)N",
        "rationale": "Simple proteinogenic amino acid; despite carrying an amine, its small zwitterionic scaffold doesn't reproduce the cationic-amphiphile pattern MRGPRX2 agonists need (contrast with beta-alanine activating the paralog MRGPRX1).",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Cholesterol",
        "label": "nonbinder",
        "smiles": "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
        "rationale": "Bulky lipophilic steroid scaffold with a single hydroxyl and no basic nitrogen; hydrophobic ring system alone isn't enough to trigger MRGPRX2 without a cationic amine.",
        "source": "curated",
        "structural_reference": True,
    },
    {
        "name": "Aspirin",
        "label": "nonbinder",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "rationale": "Small aromatic drug dominated by an anionic carboxylic acid rather than a basic amine, so it doesn't fit the cationic-pocket binding pattern shared by MRGPRX2 agonists and antagonists.",
        "source": "curated",
        "structural_reference": True,
    },
    # --- Name-only references (no compact small-molecule SMILES yet) ---
    # These participate in name-based lookups/explanations but are skipped
    # by the Tanimoto similarity search until a representative structure is
    # curated. Add a "smiles" + set structural_reference=True once available.
    {
        "name": "VIP (Vasoactive Intestinal Peptide)",
        "label": "agonist",
        "smiles": None,
        "rationale": "PACAP/VIP-family neuropeptide agonist reported alongside PACAP 1-26/1-27 to directly activate MRGPRX2 on mast cells.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "PAMP-20",
        "label": "agonist",
        "smiles": None,
        "rationale": "Longer proadrenomedullin N-terminal peptide fragment; cationic/aromatic-rich agonist in the same family as PAMP-12.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "Human beta-defensin-3",
        "label": "agonist",
        "smiles": None,
        "rationale": "Cationic antimicrobial peptide; reported MRGPRX2 agonist in the defensin family alongside human beta-defensin-2.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "CXCL14 (and cleavage products)",
        "label": "agonist",
        "smiles": None,
        "rationale": "Chemokine and its proteolytic fragments reported to directly activate MRGPRX2, linking chemokine signalling to mast-cell degranulation.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "(R)-ZINC-3573",
        "label": "agonist",
        "smiles": None,
        "rationale": "Virtual-screening hit reported as a small-molecule MRGPRX2 agonist probe compound.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "C9",
        "label": "antagonist",
        "smiles": None,
        "rationale": "Reported small-molecule MRGPRX2 antagonist scaffold used as a pharmacological tool compound to block degranulation.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "C9-6",
        "label": "antagonist",
        "smiles": None,
        "rationale": "Optimised analogue of antagonist tool compound C9 with improved MRGPRX2 blocking potency.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "EVO756",
        "label": "antagonist",
        "smiles": None,
        "rationale": "Clinical-stage oral small-molecule MRGPRX2 antagonist studied for mast-cell-driven conditions (e.g. chronic urticaria); included by name as a confirmed literature/clinical reference because no resolvable small-molecule SMILES is publicly available for it - guessing a structure here would teach the classifier a fabricated structure-activity relationship.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "GSK Compound A",
        "label": "antagonist",
        "smiles": None,
        "rationale": "GSK MRGPRX2 antagonist series compound; reported sub-to-low-nanomolar potency tool compound for receptor blockade studies.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "GSK Compound B",
        "label": "antagonist",
        "smiles": None,
        "rationale": "GSK MRGPRX2 antagonist series compound reported with IC50 ~0.42 nM, among the most potent MRGPRX2 blockers described.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "PSB-172656",
        "label": "antagonist",
        "smiles": None,
        "rationale": "Subnanomolar-potency MRGPRX2 antagonist tool compound used to dissect receptor-selective pharmacology.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "Heterotricyclic compound 1",
        "label": "antagonist",
        "smiles": None,
        "rationale": "Heterotricyclic scaffold reported as an MRGPRX2 antagonist hit from medicinal-chemistry optimisation campaigns.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "KMH-45",
        "label": "antagonist",
        "smiles": None,
        "rationale": "Reported small-molecule MRGPRX2 antagonist tool compound used in mast-cell pharmacology studies.",
        "source": "curated",
        "structural_reference": False,
    },
    {
        "name": "Oat extract",
        "label": "antagonist",
        "smiles": None,
        "rationale": "Botanical extract reported to dampen MRGPRX2-mediated mast-cell degranulation, likely via its avenanthramide/polyphenol constituents.",
        "source": "curated",
        "structural_reference": False,
    },
]


def structural_reference_ligands() -> list[dict]:
    """Entries with a resolvable SMILES, suitable for fingerprint similarity search."""
    return [
        ligand
        for ligand in MRGPRX2_REFERENCE_LIGANDS
        if ligand.get("structural_reference") and ligand.get("smiles")
    ]


def all_reference_ligands() -> list[dict]:
    """Every curated entry, including name-only references for lookups/UI display."""
    return MRGPRX2_REFERENCE_LIGANDS
