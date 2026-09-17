# Tissue geometries

Both files are static legacy VTK `POLYDATA` meshes in centimetres. They contain no membrane-voltage or other time-dependent data.

- `fiber_1d.vtk`: a 4.98 x 0.06 cm, two-element-thick representation of a one-dimensional cable; 501 nodes and 332 quadrilateral elements.
- `tissue_2d.vtk`: a 4.98 x 4.98 cm sheet; 27,889 nodes and 27,556 quadrilateral elements.

The spatial resolution is 0.03 cm (300 um). Original one-based node and element identifiers are stored in `source_node_id` and `source_element_id`; VTK connectivity is zero-based.

Cell arrays contain `material_id`, `cell_model_id` and `fiber_orientation`. Point arrays contain `S1_nodes`; the tissue mesh also contains `S2_nodes`. Masks use 1 for stimulated nodes and 0 otherwise. Detailed parameters are in `tissue_1d.yaml` and `tissue_2d.yaml` in this directory.
