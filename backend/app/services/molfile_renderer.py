from dataclasses import dataclass
import math


@dataclass
class Atom:
    x: float
    y: float
    element: str


@dataclass
class Bond:
    start: int
    end: int
    order: int


ELEMENT_COLORS = {
    "C": "#3a3a3c",
    "N": "#2563eb",
    "O": "#dc2626",
    "H": "#9ca3af",
    "S": "#b45309",
    "P": "#7c3aed",
    "F": "#0f766e",
    "CL": "#15803d",
    "BR": "#7c2d12",
}


def _parse_counts_line(line: str) -> tuple[int, int]:
    atom_count = int(line[0:3].strip())
    bond_count = int(line[3:6].strip())
    return atom_count, bond_count


def _parse_atom_line(line: str) -> Atom:
    return Atom(
        x=float(line[0:10].strip()),
        y=float(line[10:20].strip()),
        element=line[31:34].strip() or "C",
    )


def _parse_bond_line(line: str) -> Bond:
    return Bond(
        start=int(line[0:3].strip()) - 1,
        end=int(line[3:6].strip()) - 1,
        order=max(1, int(line[6:9].strip() or "1")),
    )


def molfile_to_highlight_svg(molfile: str) -> str:
    lines = [line.rstrip("\n") for line in molfile.splitlines()]
    if len(lines) < 5:
        raise ValueError("Invalid molfile: too few lines")

    atom_count, bond_count = _parse_counts_line(lines[3])
    atom_lines = lines[4 : 4 + atom_count]
    bond_lines = lines[4 + atom_count : 4 + atom_count + bond_count]

    atoms = [_parse_atom_line(line) for line in atom_lines]
    bonds = [_parse_bond_line(line) for line in bond_lines]

    visible_indices = [
        index for index, atom in enumerate(atoms) if atom.element.upper() != "H"
    ]
    if not visible_indices:
        visible_indices = list(range(len(atoms)))

    remap = {old_index: new_index for new_index, old_index in enumerate(visible_indices)}
    filtered_atoms = [atoms[index] for index in visible_indices]
    filtered_bonds = [
        Bond(
            start=remap[bond.start],
            end=remap[bond.end],
            order=bond.order,
        )
        for bond in bonds
        if bond.start in remap and bond.end in remap
    ]

    xs = [atom.x for atom in filtered_atoms]
    ys = [atom.y for atom in filtered_atoms]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    width = 360.0
    height = 320.0
    padding = 36.0
    span_x = max(max_x - min_x, 1.0)
    span_y = max(max_y - min_y, 1.0)
    scale = min((width - padding * 2) / span_x, (height - padding * 2) / span_y)

    def scale_point(x: float, y: float) -> tuple[float, float]:
        sx = (x - min_x) * scale + padding
        sy = height - ((y - min_y) * scale + padding)
        return sx, sy

    positions = [scale_point(atom.x, atom.y) for atom in filtered_atoms]
    degrees = [0 for _ in filtered_atoms]
    for bond in filtered_bonds:
        degrees[bond.start] += 1
        degrees[bond.end] += 1

    step = max(math.ceil(len(filtered_atoms) / 10), 4)
    parts = [
        '<svg class="compound-structure-graph" viewBox="0 0 360 320" xmlns="http://www.w3.org/2000/svg" fill="none">',
        '<rect width="360" height="320" rx="28" fill="#fbfbfd"/>',
    ]

    for bond in filtered_bonds:
        x1, y1 = positions[bond.start]
        x2, y2 = positions[bond.end]
        dashed = ' stroke-dasharray="5 5"' if bond.order > 1 else ""
        parts.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="#c9ccd3" stroke-width="2.6"{dashed} stroke-linecap="round" />'
        )

    for index, atom in enumerate(filtered_atoms):
        x, y = positions[index]
        element = atom.element.upper()
        fill = ELEMENT_COLORS.get(element, "#3a3a3c")
        highlight = (
            element in {"N", "O", "S", "P"}
            or (element == "C" and degrees[index] >= 3 and index % step == 0)
        )

        if highlight:
            ring_color = "#60a5fa" if element in {"N", "C"} else "#fb7185"
            parts.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="20" fill="{ring_color}" opacity="0.10" />'
            )
            parts.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="18" stroke="{ring_color}" stroke-width="2" opacity="0.75" />'
            )

        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="11.5" fill="#ffffff" opacity="0.94" />'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{y + 5:.1f}" text-anchor="middle" font-family="ui-sans-serif, system-ui, sans-serif" '
            f'font-size="14" font-weight="700" fill="{fill}">{element}</text>'
        )

    parts.append("</svg>")
    return "".join(parts)
