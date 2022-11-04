from glyles import Glycan
from PIL import Image, ImageDraw, ImageFont


class Tree:
    def __init__(self, glycan):
        self.nx_tree = glycan.parse_tree
        self.root = Node(0, 0, self.nx_tree)

    def assign_coords(self):
        self.root.assign_coords()
        return self

    def adjust_coords(self, offset):
        self.root.adjust_coords(offset)

    def get_bounds(self):
        return 0, 0, self.root.get_max_depth(), self.root.get_lower_bound()[-1][1]

    def draw_edges(self, img, **params):
        self.root.draw_edges(img, **params)

    def draw_nodes(self, img, glycan, **params):
        self.root.draw_nodes(img, glycan.parse_tree, **params)


class Node:
    def __init__(self, idx, depth, parse_tree):
        self.idx = idx
        self.depth = depth
        self.children = []
        self.coord = None

        for _, cidx in parse_tree.edges(self.idx):
            self.children.append((
                Node(cidx, self.depth + 1, parse_tree),
                parse_tree.get_edge_data(self.idx, cidx)["type"]
            ))

        self.children = list(
            sorted(
                self.children,
                key=lambda child: (child[0].get_max_depth(), child[0].get_num_children()),
                reverse=True
            )
        )

    def assign_coords(self):
        self.coord = 0
        if len(self.children) == 0:
            return 0

        for child, _ in self.children:
            child.assign_coords()

        for i in range(1, len(self.children)):
            lower = self.children[i - 1][0].get_lower_bound()
            upper = self.children[i][0].get_upper_bound()
            offset = max(l - u for (_, l), (_, u) in zip(lower, upper)) + 1
            self.children[i][0].adjust_coords(offset)

        self.coord = sum(child.coord for child, _ in self.children) / len(self.children)

        return self.coord

    def adjust_coords(self, offset):
        self.coord += offset
        for child, _ in self.children:
            child.adjust_coords(offset)

    def get_lower_bound(self):
        if len(self.children) != 0:
            return [(self.depth, self.coord)] + self.children[-1][0].get_lower_bound()
        return [(self.depth, self.coord)]

    def get_upper_bound(self):
        if len(self.children) != 0:
            return [(self.depth, self.coord)] + self.children[0][0].get_upper_bound()
        return [(self.depth, self.coord)]

    def get_num_children(self):
        return len(self.children)

    def get_max_depth(self):
        if len(self.children) == 0:
            return self.depth
        return max(child.get_max_depth() for child, _ in self.children)

    def get_coords(self):
        return self.depth, self.coord

    def draw_edges(self, img, **params):
        x, y = self.get_coords()
        k = 0.3
        l = 0.3
        for child, info in self.children:
            c_x, c_y = child.get_coords()
            img.line([
                ((x + 0.5) * params["width"], (y + 0.5) * params["height"]),
                ((c_x + 0.5) * params["width"], (c_y + 0.5) * params["height"])
            ], fill="black", width=params["line"])
            child.draw_edges(img, **params)
            img.text(
                ((x + k + 0.48) * params["width"], (c_y + (y - c_y) * (1 - l) + 0.35) * params["width"]),
                text=info[1],
                fill="black",
                anchor="lb",
                # font=ImageFont.truetype("Arial.ttf", size=15)
            )
            img.text(
                ((x - k + 1.48) * params["width"], (c_y + (y - c_y) * l + 0.35) * params["width"]),
                text=info[-2],
                fill="black",
                anchor="lb",
                # font=ImageFont.truetype("Arial.ttf", size=15),
            )

    def draw_nodes(self, img, tree, **params):
        draw_glycan(img, *(self.get_coords()), *snfg_fy(tree.nodes[self.idx]["type"]), **params)
        for child, _ in self.children:
            child.draw_nodes(img, tree, **params)


colors = {
    "Hexose": {  # same for Hexose, HexNAc, Hexosamine, Hexuronate, Unknown
        "Glc": "blue",
        "Man": "green",
        "Gal": "yellow",
        "Gul": "orange",
        "Alt": "lightpink",
        "All": "BlueViolet",
        "Tal": "LightSkyBlue",
        "Ido": "sienna",  # brown
    },
    "Deoxyhexose": {
        "Qui": "blue",
        "Rha": "green",
        "6dGul": "orange",
        "6dAlt": "lightpink",
        "6dTal": "LightSkyBlue",
        "Fuc": "red",  # brown
    },
    "DeoxyhexNAc": {
        "QuiNAc": "blue",
        "RhaNAc": "green",
        "6dAltNAc": "lightpink",
        "6dTalNAc": "LightSkyBlue",
        "FucNAc": "red",  # brown
    },
    "Di-deoxyhexose": {
        "Oli": "blue",
        "Tyv": "green",
        "Abe": "orange",
        "Par": "lightpink",
        "Dig": "BlueViolet",
        "Col": "LightSkyBlue",
    },
    "Pentose": {
        "Ara": "green",
        "Lyx": "yellow",
        "Xyl": "orange",
        "Rib": "lightpink",
    },
    "3-dideoxy-nunolusonic acids": {
        "Kdn": "green",
        "Neu5Ac": "BlueViolet",
        "neu5Gc": "LightSkyBlue",
        "Neu": "sienna",  # brown
        "Sia": "red",
    },
    "3,9-dideoxy-nunolusonic acids": {
        "Pse": "green",
        "Leg": "yellow",
        "Aci": "lightpink",
        "4eLeg": "LightSkyBlue",
    },
    "Unknown": {  # same for Hexose, HexNAc, Hexosamine, Hexuronate, Unknown
        "Bac": "blue",
        "LDmanHep": "green",
        "Kdo": "yellow",
        "Dha": "orange",
        "DDmanHep": "lightpink",
        "MurNAc": "BlueViolet",
        "MurNGc": "LightSkyBlue",
        "Mur": "sienna",  # brown
    },
    "Assigned": {  # same for Hexose, HexNAc, Hexosamine, Hexuronate, Unknown
        "Api": "blue",
        "Fru": "green",
        "Tag": "yellow",
        "Sor": "orange",
        "Psi": "lightpink",
    },
}


classes = {
    "Hexose": set(colors["Hexose"].keys()),
    "HexNAc": set([x + "NAc" for x in colors["Hexose"].keys()]),
    "Hexosamine": set([x + "N" for x in colors["Hexose"].keys()]),
    "Hexuronate": set([x + "A" for x in colors["Hexose"].keys()]),
    "Deoxyhexose": {"Qui", "Rha", "6dGul", "6dAlt", "6dTal", "Fuc"},
    "DeoxyhexNAc": {"QuiNAc", "RhaNAc", "6dAltNAc", "6dTalNAc", "FucNAc"},
    "Di-deoxyhexose": {"Oli", "Tyv", "Abe", "Par", "Dig", "Col"},
    "Pentose": {"Ara", "Lyx", "Xyl", "Rib"},
    "3-deoxy-nonulosonic acids": {"Kdn", "Neu5Ac", "Neu5Gc", "Neu", "Sia"},
    "3,9-dideoxy-nonulosonic acids": {"Pse", "Leg", "Aci", "4eLeg"},
    "Unknown": {"Bac", "LDmanHep", "Kdo", "Dha", "DDmanHep", "MurNAc", "MurNGc", "Mur"},
    "Assigned": {"Api", "Fru", "Tag", "Sor", "Psi"},
}


def snfg_fy(glycan):
    name = glycan.get_name()
    names = set(list(zip(*(glycan.get_recipe())))[0])
    if name in {
        "Oli", "Tyv", "Abe", "Par", "Dig", "Col", "Ara", "Lyx", "Xyl", "Rib", "Kdn", "Sia", "Pse", "Aci",
        "Bac", "Kdo", "Dha", "Api", "Fru", "Tag", "Sor", "Psi"
    }:
        return name, names.difference({"NAc"})
    if name in classes["Hexose"] and ("NAc" in names or "N" in names or "A" in names):
        if "NAc" in names and not (name in {"Alt", "Tal"} and "6d" in names):
            return name + "NAc", names.difference({"NAc", name})
        if "N" in names:
            return name + "N", names.difference({"N", name})
        if "A" in names:
            return name + "A", names.difference({"A", name})
    if name == "Gul" and "6d" in names:
        return "6dGul", names.difference({"6d", name})
    if name in {"Alt", "Tal"} and ("NAc" in names or "6d" in names):
        if "NAc" in names and "6d" in names:
            return "6d" + name + "NAc", names.difference({"6d", "NAc", name})
        if "NAc" in names:
            return name + "NAc", names.difference({"NAc", name})
        if "6d" in names:
            return "6d" + name, names.difference({"6d", name})
    if name in {"Qui", "Rha", "Fuc", "Mur"} and "NAc" in names:
        return name + "NAc", names.difference({"NAc", name})
    if name == "Neu" and "5Ac" in names or "Ac" in names:
        return "Neu5Ac", names.difference({"5Ac", name})
    if name == "Neu" and "5Gc" in names or "Gc" in names:
        return "Neu5Gc", names.difference({"5Gc", name})
    if name == "Mur" and "NGc" in names:
        return "Mur5Gc", names.difference({"5Gc", name})
    if name == "Leg" and "4e" in names:
        return "4eLeg", names.difference({"4e", name})
    # TODO: LDmanHep and DDmanHep
    return name, names.difference({name})


def draw_glycan(img, x, y, name, groups, **params):
    if name in classes["Hexose"]:
        draw_hexose(img, x, y, name, **params)
    elif name in classes["HexNAc"]:
        draw_hexnac(img, x, y, name, **params)
    elif name in classes["Hexosamine"]:
        draw_hexosamine(img, x, y, name, **params)
    elif name in classes["Hexuronate"]:
        draw_hexuronate(img, x, y, name, **params)
    elif name in classes["Deoxyhexose"]:
        draw_deoxyhexose(img, x, y, name, **params)
    elif name in classes["DeoxyhexNAc"]:
        draw_deoxyhexnac(img, x, y, name, **params)
    elif name in classes["Di-deoxyhexose"]:
        draw_di_deoxyhexose(img, x, y, name, **params)
    elif name in classes["Pentose"]:
        draw_pentose(img, x, y, name, **params)
    elif name in classes["3-deoxy-nonulosonic acids"]:
        draw_3_deoxy_nonulosonic_acids(img, x, y, name, **params)
    elif name in classes["3,9-dideoxy-nonulosonic acids"]:
        draw_3_9_dideocy_nonulosonic_acids(img, x, y, name, **params)
    elif name in classes["Unknown"]:
        draw_unknown(img, x, y, name, **params)
    elif name in classes["Assigned"]:
        draw_assigned(img, x, y, name, **params)
    else:
        draw_other(img, x, y, name, **params)
    text, width = arange_node_groups(groups, img, **params)
    img.multiline_text(
        ((x + 0.5) * params["width"] - width * 0.5, (y + 0.8) * params["height"]),
        text,
        fill="black",
        anchor="ma",
        align="center"
    )


def arange_node_groups(groups, img, **params):
    output = ""
    length = 0
    threshold = 0.5 * params["width"]
    lines = []
    for g in groups:
        if length == 0:
            output += g
            tl = img.textlength(g)
            if tl > threshold:
                output += "\n"
                lines.append(tl)
            else:
                length = tl
        else:
            tl = img.textlength(g)
            if length + tl > threshold:
                output += "\n" + g
                lines.append(length)
                length = tl
            else:
                output += g
                length += tl
    lines.append(length)
    return output, max(lines)


def draw_hexose(img, x, y, name, **params):
    img.ellipse([
        ((x + 0.25) * params["width"], (y + 0.25) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.75) * params["height"])
    ], fill=colors["Hexose"][name], width=params["stroke"], outline="black")


def draw_hexnac(img, x, y, name, **params):
    img.rectangle([
        ((x + 0.25) * params["width"], (y + 0.25) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.75) * params["height"])
    ], fill=colors["Hexose"][name[:3]], width=params["stroke"], outline="black")


def draw_hexosamine(img, x, y, name, **params):
    img.polygon([
        ((x + 0.25) * params["width"], (y + 0.25) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.75) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.25) * params["height"])
    ], fill=colors["Hexose"][name[:3]], width=params["stroke"], outline="black")
    img.polygon([
        ((x + 0.25) * params["width"], (y + 0.25) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.75) * params["height"]),
        ((x + 0.25) * params["width"], (y + 0.75) * params["height"])
    ], fill="white", width=params["stroke"], outline="black")


def draw_hexuronate(img, x, y, name, **params):
    if name in {"IdoA", "AltA"}:
        cols = ("white", colors["Hexose"][name[:3]])
    else:
        cols = (colors["Hexose"][name[:3]], "white")
    img.polygon([
        ((x + 0.25) * params["width"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.25) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.5) * params["height"])
    ], fill=cols[0], width=params["stroke"], outline="black")
    img.polygon([
        ((x + 0.25) * params["width"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.75) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.5) * params["height"])
    ], fill=cols[1], width=params["stroke"], outline="black")


def draw_deoxyhexose(img, x, y, name, **params):
    img.polygon([
        ((x + 0.25) * params["width"], (y + 0.75) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.25) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.75) * params["height"])
    ], fill=colors["Deoxyhexose"][name], width=params["stroke"], outline="black")


def draw_deoxyhexnac(img, x, y, name, **params):
    img.polygon([
        ((x + 0.25) * params["width"], (y + 0.75) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.25) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.75) * params["height"])
    ], fill="white", width=params["stroke"], outline="black")
    img.polygon([
        ((x + 0.5) * params["width"], (y + 0.75) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.25) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.75) * params["height"])
    ], fill=colors["DeoxyhexNAc"][name], width=params["stroke"], outline="black")


def draw_di_deoxyhexose(img, x, y, name, **params):
    img.rectangle([
        ((x + 0.25) * params["width"], (y + 0.375) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.625) * params["height"])
    ], fill=colors["Di-deoxyhexose"][name], width=params["stroke"], outline="black")


def draw_pentose(img, x, y, name, **params):
    img.polygon([
        ((x + 0.5) * params["width"], (y + 0.25) * params["height"]),

        ((x + 0.55902) * params["width"], (y + 0.43164) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.43164) * params["height"]),

        ((x + 0.5955) * params["width"], (y + 0.5439) * params["height"]),
        ((x + 0.66877) * params["width"], (y + 0.75) * params["height"]),

        ((x + 0.5) * params["width"], (y + 0.61328) * params["height"]),

        ((x + 0.33123) * params["width"], (y + 0.75) * params["height"]),
        ((x + 0.4045) * params["width"], (y + 0.5439) * params["height"]),

        ((x + 0.25) * params["width"], (y + 0.43164) * params["height"]),
        ((x + 0.44098) * params["width"], (y + 0.43164) * params["height"])
    ], fill=colors["Pentose"][name], width=params["stroke"], outline="black")


def draw_3_deoxy_nonulosonic_acids(img, x, y, name, **params):
    img.polygon([
        ((x + 0.25) * params["width"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.25) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.75) * params["height"])
    ], fill=colors["3-dideoxy-nunolusonic acids"][name], width=params["stroke"], outline="black")


def draw_3_9_dideocy_nonulosonic_acids(img, x, y, name, **params):
    img.polygon([
        ((x + 0.25) * params["width"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.375) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.625) * params["height"])
    ], fill=colors["3,9-dideoxy-nunolusonic acids"][name], width=params["stroke"], outline="black")


def draw_unknown(img, x, y, name, **params):
    img.polygon([
        ((x + 0.25) * params["width"], (y + 0.5) * params["height"]),
        ((x + 0.375) * params["width"], (y + 0.375) * params["height"]),
        ((x + 0.625) * params["width"], (y + 0.375) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.5) * params["height"]),
        ((x + 0.625) * params["width"], (y + 0.625) * params["height"]),
        ((x + 0.375) * params["width"], (y + 0.625) * params["height"]),
    ], fill=colors["Unknown"][name], width=params["stroke"], outline="black")


def draw_assigned(img, x, y, name, **params):
    img.polygon([
        ((x + 0.5) * params["width"], (y + 0.25) * params["height"]),
        ((x + 0.75) * params["width"], (y + 0.43164) * params["height"]),
        ((x + 0.66877) * params["width"], (y + 0.75) * params["height"]),
        ((x + 0.33123) * params["width"], (y + 0.75) * params["height"]),
        ((x + 0.25) * params["width"], (y + 0.43164) * params["height"])
    ], fill=colors["Assigned"][name], width=params["stroke"], outline="black")


def draw_other(img, x, y, name, groups, **params):
    pass


def create_snfg_img(glycan, filepath, **kwargs):
    params = {
        "width": 100,
        "height": 100,
        "stroke": 2,
        "line": 3,
    }
    params.update(**kwargs)
    tree = Tree(glycan).assign_coords()
    width, height = tree.get_bounds()[2:]

    img = Image.new(
        mode="RGB",
        size=(int((width + 1) * params["width"]), int((height + 1) * params["height"])),
        color=(255, 255, 255)
    )
    draw_img = ImageDraw.Draw(img)
    tree.draw_edges(draw_img, **params)
    tree.draw_nodes(draw_img, glycan, **params)
    img.save(filepath)
    return img


def main():
    iupac = "Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)GlcNAc(b1-3)[Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)GlcNAc(b1-6)]" \
            "Gal(b1-3)[GlcNAc(a1-4)Gal(b1-4)GlcNAc6S(b1-6)]GalNAc"
    iupac2 = "IdoA2S(b1-4)Fuc2S(a1-4)[Par2S(a1-4)Xyl2S(a1-2)][All2S(a1-3)TalNAc3S(b1-6)GulN3S(a1-3)6dAltNAc3S(b1-3)]Sia6S(a1-4)Pse2S(a1-4)[Gal2S(a1-3)]Bac2S(a1-4)Tag2S"
    glycan = Glycan(iupac2, tree_only=True)
    create_snfg_img(glycan, **{
        "width": 100,
        "height": 100,
        "stroke": 2,
        "line": 3,
    })


if __name__ == '__main__':
    main()
