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
"Prefuse Force Directed layouot" with Default Spring Length = 20 and Default Node Mass = 1 (Settings) often gives nice results. But feel free to explore other settings and layouts.

##### Data import (fluxes) to Recon3D - Approch I
- base Recon3D model: RECON3D_consistent_model.mat; RECON3D_consistent_model.xml; RECON3D_consistent_model.cys
- RECON3D model with pre-defined size and shape mapping: RECON3D_061123.cys
- metsFluxSum_v3FBA_061123.m: generates table (metsFluxSum_log2.txt) with metabolite/reaction name (sbml id), node size (based on flux sum per metabolite or absolute flux per reaction), shape and keep/remove flag
- shape is "Elipse" for metabolites above cutoff and "none" for metabolites below cutoff
- shape is "Triangle"/"V"/"none" for reactions with positive/negative/zero flux
- "none" shapes get "remove" flag, others "keep"
- Cytoscape -> open RECON3D_061123.cys
- import to cytoscape via: File -> Import -> Table from File -> (Select file name: metsFluxSum_log2.txt). Then: Where to Import Table Data:To a Network Collection; Import Data as: Node Table Columns; Key Column For Network: shared name. Advanced Options: Delimiter: comma. OK.
- size and shape mapping should automatically update network (RECON3D_061123_dataMapped.cys)
- Keep onnly active metabolites and reactions: Node table -> search "keep" -> invert node selection -> delete (RECON3D_061123_dataMapped_keepOnly.cys)
- Remove cofactors: Select -> Nodes -> from ID List File -> (File name of above list: metsCofactors.txt) -> delete (RECON3D_061123_dataMapped_keepOnly_removeCo.cys)

##### Data import (fluxes) to Recon3D - Approch II
- addfluxes2.m: generates table with reaction name (shared name) and data
- import to cytoscape via: File -> Import -> Table from File -> (Select file name). Then: Where to Import Table Data:To a Network Collection; Import Data as: Node Table Columns; Key Column For Network: shared name. Advanced Options: Delimiter: comma. OK.
- map flux data e.g. via Style -> Size -> (select data column) -> Continuous Mapping. Double-click on Mapping -> Add (to add additional base point) -> Node Size = 30 (for high & low flux values) & = 5 (for 0 flux values)
- See: Recon3D_cyto_split_v2_dataMapped

###### Keep only reactions and metabolites of interest - Approch II
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
