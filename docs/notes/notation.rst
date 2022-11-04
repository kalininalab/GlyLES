Notation
========

This package is able to convert a special form of IUPAC strings for glycans into SMILES strings. To convert normal
IUPAC notation for organic molecules into SMILES there are several tools already published, e.g.,
`OPSIN <https://opsin.ch.cam.ac.uk/>`.

Short Overview
--------------

In the following, we will show four different ways to describe lactose, a disaccharide of D-Galactose and D-Glucose
connected with a 1-4 beta-O-glycosidic bond as shown in the image beyond.

.. image:: lactose.png
    :width: 600
    :alt: Lactose [Image taken from: https://de.wikipedia.org/wiki/Lactose]

* IUPAC:
    `beta-D-galacto-hexopyranosyl-(1->4)-beta-D-gluco-hexopyranose`
* IUPAC-condensed:
    `Gal(b1-4)Glc b`
* IUPAC-full:
    `Galp(b1-4)Glcpb`
* IUPAC-simplified:
    `Galb4Glc b`

All four notations describe the same molecule. For the standard IUPAC, there are published tools available. For the
other three, we developed this tool. For the sake of simplicity, we will only describe IUPAC-condensed here. The other
two IUPAC derivates can be derived by either giving all information on the glycans (full) or by reducing the bond
description to the absolute minimum (simplified).

IUPAC-condensed
---------------

The IUPAC-condensed notation of glycans is the most commonly used notation in publications, databases, and other
sources. Its wide usage is due to the compact descriptions of glycans and the human-readability of the notation.

This notation is a specialized form of the general IUPAC notation that is used to describe other organic molecules in a
standardized way. It is specifically suited for glycans and their variability in structure and attached functional
groups. It is capable of representing these structures in a compact and human-readable form.

Chains
^^^^^^

Glycans have a complex, yet regular, structure, especially when described using the IUPAC notation.
For example, `Man(a1-4)Man(a1-4)Man(a1-4)Man` is a simple chain of four mannopyranose monosaccharides connected with
1-4 alpha-O-glycosidic bonds. Important to note is that the root of the glycan tree is always the last monosaccharide
in an IUPAC formula. The chain grows from the end by prepending elements like \verb+Man(a1-4)+, which corresponds to
adding a new mannose to a leaf-monomer of the glycan.

Branches
^^^^^^^^

The branching of the trees is described by introducing the new branch in square brackets. For example,
`Man(a1-4)Man(a1-3)[Man(a1-4)Man(a1-4)]Man(a1-4)Man` is a tree of mannopyranoses that has two monosaccharides before
splitting into two chains with again two mannopyranoses each. The branching is put right in front of the monomer where
the root of the side branch is bound to the main strand.

Even more branches
^^^^^^^^^^^^^^^^^^

This can be extended to three chains of two monosaccharides each bound to a single mannopyranose
`Man(a1-4)Man(a1-2)[Man(a1-4)Man(a1-3)][Man(a1-4)Man(a1-4)]Man`

Attachments, a.k.a functional groups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Modifications in IUPAC are directly annotated at the modified monosaccharide. The annotation is often done by first
naming the position at which carbon atom a functional group is attached. Then, an abbreviation of that group is given.
So, `Man3S` describes a mannose with a sulfate attached to the third carbon atom. The number can also be dropped in
case of a functional group attached to a standard position, as in `GalNAc` where an acetamid-group is attached to the
second carbon atom of galactopyranose. There are also other possibilities to denote a functional group but that is
beyond the scope of this introduction.

Variety of Attachments
^^^^^^^^^^^^^^^^^^^^^^
Modifications of glycans include added sulfur groups (`Man3S`), added phosphate groups (`Man3P`), added amine groups
(`ManN`), and many more. Functional groups can be attached to any oxygen, nitrogen, or a carbon atom in the
monosaccharide (e.g., `Man4S`, `Man4P`, or `ManNS`). Additionally, one can combine them in any way (e.g., `Man3S4P`).
The only restriction is that one cannot put two modifications to the same atom (e.g., `Man2S2P`).
