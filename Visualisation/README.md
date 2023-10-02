# Visualisation of Metabolic Networks and related Data

### Cytoscape and IDARE2 cytoscape app
Cytoscape:

https://cytoscape.org/

https://github.com/cytoscape/cytoscape-tutorials/wiki

##### IDARE2:

https://www.mdpi.com/2218-1989/11/5/300

https://github.com/sysbiolux/IDARE

We recommend to start here:

https://github.com/sysbiolux/IDARE-QuickStart

Examples:
- Recon3D (Recon3D_cyto_split_v2, with link & reaction Label Font Size = 5)
- Blautia hydrogenotropica

##### Layout
"Prefuse Force Directed layouot" with Default Spring Length = 20 and Default Node Mass = 1 (Settings) gives aften nice results. But feel free to explore other settings and layouts.

##### Data import (fluxes) to Recon3D
- addfluxes2.m: generates table with reaction name (shared name) and data
- import to cytoscape via: File -> Import -> Table from File -> (Select file name). Then: Where to Import Table Data:To a Network Collection; Import Data as: Node Table Columns; Key Column For Network: shared name. Advanced Options: Delimiter: comma. OK.
- map flux data e.g. via Style -> Size -> (select data column) -> Continuous Mapping. Double-click on Mapping -> Add (to add additional base point) -> Node Size = 30 (for high & low flux values) & = 5 (for 0 flux values)
- See: Recon3D_cyto_split_v2_dataMapped

##### Keep only reactions and metabolites of interest
- Generate list of reactions to be removed (with a flux value between -10 and 10 only)
- Cytoscape: Select -> Nodes -> from ID List File -> (File name of above list)
- Delete (these reactions – and metabolites -)
- Remove orphan metabolites: Select all edges. Select nodes connected by selected edges. Invert node selection. Deselect all edges. Delete.
- Prefuse Force Directed Layout.
- Apps-> Create Subnetworks -> … (confirm) … -> Column to determine subnetwork: Transport, lysosomal (=subsystem); Select all (Subnetworks)
- See: Recon3D_cyto_v2_dataMapped_abs10fluxes; Recon3D_cyto_v2_dataMapped_abs10fluxesC6-3

### Escher (for small networks and pathways)
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004321

https://escher.readthedocs.io/en/latest/

Online: https://escher.github.io/#/

Python: https://escher.readthedocs.io/en/latest/escher-python.html
