from PIL.ImageDraw import Draw


class Tree:
    """
    Tree class to storeing a graphical representation of glycans in SNFG form.
    """

    def __init__(self, glycan):
        """
        Initialize the tree with the glycan.

        Args:
            glycan (Glycan): Glycan object to be represented in this tree object
        """
        self.root = Node(0, 0, glycan.parse_tree)

    def assign_coords(self):
        """
        Assign coordinates to every node in this tree.

        Returns:
            It-self
        """
        self.root.assign_coords()
        return self

    def adjust_coords(self, offset):
        """
        Adjust the coordinates by adjusting every node's coordinate.

        Args:
            offset (float): Offset by which all nodes should be moved
        """
        self.root.adjust_coords(offset)

    def get_bounds(self):
        """
        Bounding box in terms of node coordinates on a graphical representation.

        Returns:
            Tuple of 4 floats representing upper left and bottom right corners
        """
        return 0, 0, self.root.get_max_depth(), self.root.get_lower_bound()[-1][1]

    def draw_edges(self, img, **params):
        """
        Draw the edges of the tree into the ImageDraw object.

        Args:
            img (Draw): ImageDraw object to draw on
            **params: parameters for visualization
        """
        self.root.draw_edges(img, **params)

    def draw_nodes(self, img, glycan, **params):
        """
        Draw the nodes of the tree representing the individual monosaccharides in a sugar polymer.

        Args:
            img (Draw): ImageDraw object to draw on
            glycan (Glycan): Glycan object holding some information on the sugar monomers
            **params: parameters for visualization
        """
        self.root.draw_nodes(img, glycan.parse_tree, **params)


class Node:
    """
    A single node in representing a monosaccharide in a graphical SNFG representation.
    """
    def __init__(self, idx, depth, parse_tree):
        """
        Initialize a node representing a monosaccharide.

        Args:
            idx (int): Id of the glycan that matches the id in the networkx tree of the glycan object in the parse_tree.
            depth (int): Depth of the node in the tree, root has depth 0, leaves have the highest depth-value
            parse_tree (networkx.DiGraph): networkx graph representation of the glycan computed while parsing
        """
        self.idx = idx
        self.depth = depth
        self.children = []
        self.coord = None

        # iterate over the children and recursively initialize them.
        for _, cidx in parse_tree.edges(self.idx):
            self.children.append((
                Node(cidx, self.depth + 1, parse_tree),
                parse_tree.get_edge_data(self.idx, cidx)["type"]
            ))

        # sort the children based on the depth of their most far away leaf and the number of children
        self.children = list(
            sorted(
                self.children,
                key=lambda child: (child[0].get_max_depth(), child[0].get_num_children()),
                reverse=True
            )
        )

    def assign_coords(self):
        """
        Assign coordinates of this node. It is only necessary to assign y-coordinates as x coordinates are set to depth.

        Returns:
            The own y-coordinate relative to the parent.
        """
        # initially the coordinate is 0 and can be returned in case this node is a leaf
        self.coord = 0
        if len(self.children) == 0:
            return 0

        # compute the coordinates of the children
        for child, _ in self.children:
            child.assign_coords()

        # iterate over the children and arrange them to not overlap by moving the lower child down based on bounds
        for i in range(1, len(self.children)):
            lower = self.children[i - 1][0].get_lower_bound()
            upper = self.children[i][0].get_upper_bound()
            offset = max(l - u for (_, l), (_, u) in zip(lower, upper)) + 1
            self.children[i][0].adjust_coords(offset)

        # set own coordinate to be in middle height of their children
        self.coord = sum(child.coord for child, _ in self.children) / len(self.children)

        return self.coord

    def adjust_coords(self, offset):
        """
        Adjust the y-coordinates of a node by adding the offset to the coordinate and propagating the offset to the
        children.

        Args:
            offset (float): Offset to be added to the subtrees coordinates
        """
        self.coord += offset
        for child, _ in self.children:
            child.adjust_coords(offset)

    def get_lower_bound(self):
        """
        Get the lower bound of this subtree in terms of node coordinates.

        Returns:
            A List of tuples of x and y coordinates (floats) representing the lower bound for this subtree.
        """
        # compose the lower bound as the own coordinate the lower bound of the lowest child
        if len(self.children) != 0:
            return [(self.depth, self.coord)] + self.children[-1][0].get_lower_bound()
        return [(self.depth, self.coord)]

    def get_upper_bound(self):
        """
        Get the upper bound of this subtree in terms of node coordinates.

        Returns:
            A List of tuples of x and y coordinates (floats) representing the upper bound for this subtree.
        """
        # compose the upper bound as the own coordinate the upper bound of the highest child
        if len(self.children) != 0:
            return [(self.depth, self.coord)] + self.children[0][0].get_upper_bound()
        return [(self.depth, self.coord)]

    def get_num_children(self):
        """
        Get the number of children of this node.

        Returns:
            Number of children
        """
        return len(self.children)

    def get_max_depth(self):
        """
        Get the maximal depth value of a leaf of the subtree of this node. The return value will be relative to the
        root not to this node.

        Returns:
            Depth of the deepest leaf of this subtree
        """
        if len(self.children) == 0:
            return self.depth
        return max(child.get_max_depth() for child, _ in self.children)

    def get_coords(self):
        """
        Get the coordinates pair of this node.

        Returns:
            Tuple of two floats representing the x and y coordinate of this node
        """
        return self.depth, self.coord

    def draw_edges(self, img, **params):
        """
        Draw the edge between two nodes.

        Args:
            img (Draw): ImageDraw object to draw on
            **params: parameters for visualization
        """
        x, y = self.get_coords()
        k = 0.3 * (params["box"] / params["width"])
        l = 0.3 * (params["box"] / params["width"])
        for child, info in self.children:
            # draw a line from the center of this node to the centers of the child nodes
            c_x, c_y = child.get_coords()
            img.line([
                ((x + 0.5) * params["width"], (y + 0.5) * params["height"]),
                ((c_x + 0.5) * params["width"], (c_y + 0.5) * params["height"])
            ], fill="black", width=params["line"])
            child.draw_edges(img, **params)

            # for each edge denote the additional information such as a/b conformation and binding carbon id
            img.text(
                ((x + k + 0.5) * params["width"], (c_y + (y - c_y) * (1 - l) + 0.38) * params["height"]),
                text=info[1],
                fill="black",
                anchor="lb",
            )
            img.text(
                ((x - k + 1.48) * params["width"], (c_y + (y - c_y) * l + 0.38) * params["height"]),
                text=info[-2],
                fill="black",
                anchor="lb",
            )

    def draw_nodes(self, img, tree, **params):
        """
        Draw the SNFG representation of this monomer to the ImageDraw object.

        Args:
            img (Draw): ImageDraw object to draw on
            tree (networkx.DiGraph): ParseTree from a glycan, computed during parsing of the IUPAC string
            **params: parameters for visualization
        """
        # draw this node and all the children recursively
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
    """
    Convert the description of a glycan into the SNFG name and a list of attached functional groups.

    Args:
        glycan (Glycan): Glycan holding the information to be extracted

    Returns:
        Name of the glycan in SNFG-conform term and a list of remaining functional groups.
    """
    name = glycan.get_name()
    names = set(list(zip(*(glycan.get_recipe())))[0])
    if name in {
        "Oli", "Tyv", "Abe", "Par", "Dig", "Col", "Ara", "Lyx", "Xyl", "Rib", "Kdn", "Sia", "Pse", "Aci",
        "Bac", "Kdo", "Dha", "Api", "Fru", "Tag", "Sor", "Psi"
    }:
        return name, names.difference({name})
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
    if name == "Man" and "LD" in names and "Hep" in names:
        return "LDManHep", names.difference({"LD", "Man", "Hep"})
    if name == "Man" and "DD" in names and "Hep" in names:
        return "DDManHep", names.difference({"DD", "Man", "Hep"})
    return name, names.difference({name})


def draw_glycan(img, x, y, name, groups, **params):
    """
    Draw a single glycan node to the image.

    Args:
        img (Draw): ImageDraw object to paint on
        x (float): x coordinate to be drawn
        y (float): y coordinate to be drawn
        name (str): name of the monomer to be drawn
        groups (List[str]): List of strings to be written beyond the
        **params: See description in viz.create_snfg_image
    """
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
        draw_other(img, x, y, **params)

    # finally write the text under the symbol
    text, width = arrange_node_groups(groups, img, **params)
    img.multiline_text(
        ((x + 0.5) * params["width"] - width * 0.5, (y + 0.8) * params["height"]),
        text,
        fill="black",
        anchor="ma",
        align="center"
    )


def arrange_node_groups(groups, img, **params):
    """
    Arrange the names of the attached functional groups in a way to not exceed the width of the sign.

    Args:
        groups (List[str]): List of descriptions of functional groups attached to a glycan
        img (Draw): ImageDraw object to estimate the width of each textbox
        **params: See description in Glycan.create_snfg_image

    Returns:
        Text with newlines and maximal width of the resulting textbox for better centering it.
    """
    # set some variable and define a threshold for the maximal width of a line
    output = ""
    length = 0
    threshold = 0.5 * params["box"]
    lines = []

    # iterate over all groups
    for g in groups:
        # if the line is empty, add the functional group as new line and store the length in case it's below threshold
        if length == 0:
            output += g
            tl = img.textlength(g)
            if tl > threshold:
                output += "\n"
                lines.append(tl)
            else:
                length = tl
        # otherwise, check if the line has enough space for another functional group and add a new line eventually
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
    """
    Draw a filled circle for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.ellipse([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"])
    ], fill=colors["Hexose"][name], width=params["stroke"], outline="black")


def draw_hexnac(img, x, y, name, **params):
    """
    Draw a filled square for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.rectangle([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
    ], fill=colors["Hexose"][name[:3]], width=params["stroke"], outline="black")


def draw_hexosamine(img, x, y, name, **params):
    """
    Draw a crossed square for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.polygon([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
    ], fill=colors["Hexose"][name[:3]], width=params["stroke"], outline="black")
    img.polygon([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
    ], fill="white", width=params["stroke"], outline="black")


def draw_hexuronate(img, x, y, name, **params):
    """
    Draw a divided diamond for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    if name in {"IdoA", "AltA"}:
        cols = ("white", colors["Hexose"][name[:3]])
    else:
        cols = (colors["Hexose"][name[:3]], "white")
    img.polygon([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"]),
    ], fill=cols[0], width=params["stroke"], outline="black")
    img.polygon([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"]),
    ], fill=cols[1], width=params["stroke"], outline="black")


def draw_deoxyhexose(img, x, y, name, **params):
    """
    Draw a filled triangle for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.polygon([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
    ], fill=colors["Deoxyhexose"][name], width=params["stroke"], outline="black")


def draw_deoxyhexnac(img, x, y, name, **params):
    """
    Draw a divided triangle for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.polygon([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
    ], fill="white", width=params["stroke"], outline="black")
    img.polygon([
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
    ], fill=colors["DeoxyhexNAc"][name], width=params["stroke"], outline="black")


def draw_di_deoxyhexose(img, x, y, name, **params):
    """
    Draw a flat rectangle for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.rectangle([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] - 0.125 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] + 0.125 * params["box"]),
    ], fill=colors["Di-deoxyhexose"][name], width=params["stroke"], outline="black")


def draw_pentose(img, x, y, name, **params):
    """
    Draw a 5-star for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.polygon([
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] - 0.25 * params["box"]),

        ((x + 0.5) * params["width"] + 0.05902 * params["box"], (y + 0.5) * params["height"] - 0.06836 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] - 0.06836 * params["box"]),

        ((x + 0.5) * params["width"] + 0.0955 * params["box"], (y + 0.5) * params["height"] + 0.0439 * params["box"]),
        ((x + 0.5) * params["width"] + 0.16877 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),

        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] + 0.11328 * params["box"]),

        ((x + 0.5) * params["width"] - 0.16877 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
        ((x + 0.5) * params["width"] - 0.0955 * params["box"], (y + 0.5) * params["height"] + 0.0439 * params["box"]),

        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] - 0.06836 * params["box"]),
        ((x + 0.5) * params["width"] - 0.05902 * params["box"], (y + 0.5) * params["height"] - 0.06836 * params["box"]),
    ], fill=colors["Pentose"][name], width=params["stroke"], outline="black")


def draw_3_deoxy_nonulosonic_acids(img, x, y, name, **params):
    """
    Draw a filled diamond for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.polygon([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
    ], fill=colors["3-dideoxy-nunolusonic acids"][name], width=params["stroke"], outline="black")


def draw_3_9_dideocy_nonulosonic_acids(img, x, y, name, **params):
    """
    Draw a diamond for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.polygon([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] - 0.125 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] + 0.125 * params["box"]),
    ], fill=colors["3,9-dideoxy-nunolusonic acids"][name], width=params["stroke"], outline="black")


def draw_unknown(img, x, y, name, **params):
    """
    Draw a flat hexagon for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.polygon([
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"] - 0.125 * params["box"], (y + 0.5) * params["height"] - 0.125 * params["box"]),
        ((x + 0.5) * params["width"] + 0.125 * params["box"], (y + 0.5) * params["height"] - 0.125 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"]),
        ((x + 0.5) * params["width"] + 0.125 * params["box"], (y + 0.5) * params["height"] + 0.125 * params["box"]),
        ((x + 0.5) * params["width"] - 0.125 * params["box"], (y + 0.5) * params["height"] + 0.125 * params["box"]),
    ], fill=colors["Unknown"][name], width=params["stroke"], outline="black")


def draw_assigned(img, x, y, name, **params):
    """
    Draw a pentagon for an assigned glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        name (str): name of the glycan to select the correct color
        **params: See description in Glycan.create_snfg_image
    """
    img.polygon([
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] - 0.06836 * params["box"]),
        ((x + 0.5) * params["width"] + 0.16877 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
        ((x + 0.5) * params["width"] - 0.16877 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] - 0.06836 * params["box"]),
    ], fill=colors["Assigned"][name], width=params["stroke"], outline="black")


def draw_other(img, x, y, **params):
    """
    Draw an empty pentagon for an unknown glycan.

    Args:
        img (Draw): ImageDraw object to use for drawing.
        x (float): x coordinate of the shape on the image before offsetting
        y (float): y coordinate of the shape on the image before offsetting
        **params: See description in Glycan.create_snfg_image
    """
    img.polygon([
        ((x + 0.5) * params["width"], (y + 0.5) * params["height"] - 0.25 * params["box"]),
        ((x + 0.5) * params["width"] + 0.25 * params["box"], (y + 0.5) * params["height"] - 0.06836 * params["box"]),
        ((x + 0.5) * params["width"] + 0.16877 * params["box"], (y + 0.5) * params["height"] + 0.75 * params["box"]),
        ((x + 0.5) * params["width"] - 0.16877 * params["box"], (y + 0.5) * params["height"] + 0.25 * params["box"]),
        ((x + 0.5) * params["width"] - 0.25 * params["box"], (y + 0.5) * params["height"] - 0.06836 * params["box"]),
    ], fill="white", width=params["stroke"], outline="black")
