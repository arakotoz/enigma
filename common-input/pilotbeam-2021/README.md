# Description of the ROOT files

## Align params files

Alignment parameters were computed w.r.t. a given geometry:

corrections (align params) | w.r.t geometry
`mftprealignment.root` | `o2sim_geometry.root`
`pass1_mft_alignment.root` | `o2sim_geometry-prealigned.root`
`pass1_wrt_ideal_mft_alignment.root` | `o2sim_geometry.root`
`pass1b_mft_alignment.root` | `pass1_o2sim_geometry-aligned.root`
`pass2_mft_alignment.root` | `o2sim_geometry.root`

`pass1_wrt_ideal_mft_alignment.root` was obtained by convoluting the pre-alignment parameters `mftprealignment.root` with alignment parameters from MillePede `pass1_mft_alignment.root`, using the macro `exploreGeom.C` .
## Geometry files

Alignment parameters were computed w.r.t initial geometry and integrated into the final geometry:

initial geometry | corrections (align params) | final geometry
 ---------------- | -------------------------- | ---------------
 `o2sim_geometry.root` | `mftprealignment.root` | `o2sim_geometry-prealigned.root`
 `o2sim_geometry-prealigned.root` | `pass1_mft_alignment.root` | `pass1_o2sim_geometry-aligned.root`
 `pass1_o2sim_geometry-aligned.root` | `pass1b_mft_alignment.root` | `pass1b_o2sim_geometry-aligned.root`
 `o2sim_geometry.root` | `pass2_mft_alignment.root` | `pass2_o2sim_geometry-aligned.root`
 `o2sim_geometry.root`| `pass1_wrt_ideal_mft_alignment.root` | `new-pass1_o2sim_geometry-aligned.root`

### Ideal geometry file

- `o2sim_geometry.root`

### Pre-aligned geometry file

- `o2sim_geometry-prealigned.root`

### Aligned geometry files from tests of Millipede

- `pass1_o2sim_geometry-aligned.root`
- `new-pass1_o2sim_geometry-aligned.root`
- `pass1b_o2sim_geometry-aligned.root`
- `pass2_o2sim_geometry-aligned.root`

## Other files

- `o2sim_grpMagField.root` can be used as `GLO/Config/GRPMagField/snapshot.root` for a reconstruction workflow with zero magnetic field
- `ctf_dictionary.root`, `MFTdictionary.bin` are old dictionaries used before [Jira # O2-2967](https://alice.its.cern.ch/jira/browse/O2-2967)