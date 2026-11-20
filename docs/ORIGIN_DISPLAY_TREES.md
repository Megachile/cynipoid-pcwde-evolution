# Origin-tree display views

Figure 2 and Figure S1 use fixed tip sets from the previous display files, recorded in `config/origin_display_tips.tsv`. Each row maps one display label to a complete-tree label and marks the outgroup tip used to root that view. The four removed CB-derived *Bassettia pallida* tips are excluded under the accepted removal record; no other display tips are changed.

`scripts/trees/prune_origin_display_trees.R` reads the post-removal complete origin trees in `repository_inputs/deposit/base/trees/`, checks their SHA256 values against `config/logical_inputs.tsv`, roots each tree on its view's outgroup tip, and prunes to the listed tips. The pinned tree SHA256 values are GH5_2 `6d1d92bfddf1e9cf570ac325e6ab282d585d19eb79799502273e820a246bf3a2`, PL1_8 `249673331254859fb25328607722c8d776ffe6a17f4f531c5bbfc923ccca45cc`, and PL4_4 `089725874cf12f1580d9d6d1a2e084f78acd20fc3c33a32340d9687618e07375`.

Ultrafast-bootstrap values belong to the complete-tree branches. The producer transfers them across rerooting by bipartition and gives every display node the value at the complete-tree most recent common ancestor of its kept descendants. It rereads each output tree and checks that rule for every internal node. `build/trees/origin_key_nodes.tsv` records the cynipoid clade and joining node for each view, their support, and the sister lineage. Figure 2 always labels those two nodes; it labels other internal nodes at support 70 or greater, except two-tip cherries. Figure S1 labels internal nodes at support 70 or greater.

The renderer consumes only the `build/trees` intermediates. `scripts/release/check_origin_key_nodes.py` checks that the Figure 2 labels and `build/numbers.json` agree with the producer's key-node table. The displayed branch lengths and supports come from the complete-tree inference, not a new inference on the reduced tip set.

### Figure S1 tip order

The GH5_2 Figure S1 cherry containing *Euceroptres primus* and *Belonocnema kinseyi* follows the order of those two rows in `config/origin_display_tips.tsv` (Euceroptres, then Belonocnema). This is a presentation setting inherited from the frozen S1 view. The producer swaps only the two sibling edge rows and their matching branch lengths when needed; the support at their parent and every clade remain unchanged. The frozen tree is not read during the build.
