# Description of the files

## Align params

- `mftprealignment.root`
- `pass1_mft_alignment.root`
- `pass1b_mft_alignment.root`
- `pass2_mft_alignment.root`

## Geometry

initial geometry | corrections (align params) | final geometry
 ---------------- | -------------------------- | ---------------
 `o2sim_geometry.root` | `mftprealignment.root` | `o2sim_geometry-prealigned.root`
 `o2sim_geometry-prealigned.root` | `pass1_mft_alignment.root` | `pass1_o2sim_geometry-aligned.root`
 `pass1_o2sim_geometry-aligned.root` | `pass1b_mft_alignment.root` | `pass1b_o2sim_geometry-aligned.root`
 `o2sim_geometry.root` | `pass2_mft_alignment.root` | `pass2_o2sim_geometry-aligned.root`

### Ideal

- `o2sim_geometry.root`

### Pre-alignment

- `o2sim_geometry-prealigned.root`

### Tests of Millipede

- `pass1_o2sim_geometry-aligned.root`
- `pass1b_o2sim_geometry-aligned.root`
- `pass2_o2sim_geometry-aligned.root`

## Other files

- `o2sim_grpMagField.root` can be used as `GLO/Config/GRPMagField/snapshot.root` for a reconstruction workflow with zero magnetic field
- `ctf_dictionary.root`, `MFTdictionary.bin` are old dictionaries used before [Jira # O2-2967](https://alice.its.cern.ch/jira/browse/O2-2967)