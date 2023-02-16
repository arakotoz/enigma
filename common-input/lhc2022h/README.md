# Description of the ROOT files

## Align params files

Alignment parameters were computed w.r.t. a given geometry:

corrections (align params) | w.r.t geometry
-------------------------- | ---------------
`pass2_mft_alignment.root` | `o2sim_geometry.root`

`pass2_mft_alignment.root` was obtained with Mille records built with 12812194 tracks from two B Off runs:
- 4923670 tracks from pilot beam 2021 run 505713 (regenerated CTFs with cluster dictionary from data)
- and 1981840+5906684 tracks from LHC22h 

MillePede2 solver output from iteration 1 and iteration 10:
- Obtained 281060922 equations from 12812193 records (0       rejected). Fixed 80   globals
- Obtained 174846585 equations from 7931371 records (4880822 rejected). Fixed 80   globals

## Geometry files

Alignment parameters were computed w.r.t initial geometry and integrated into the final geometry:

initial geometry | corrections (align params) | final geometry
 ---------------- | -------------------------- | ---------------
 `o2sim_geometry.root` | `pass2_mft_alignment.root` | `pass2_o2sim_geometry-aligned.root`
 
### Ideal geometry file

- `o2sim_geometry.root`


### Aligned geometry files from tests of Millipede

- `pass2_o2sim_geometry-aligned.root`

## Other files

- `o2__itsmft__TopologyDictionary_1653153860732.root` is the new cluster dictionary made from real data, as detailed in [Jira # O2-2967](https://alice.its.cern.ch/jira/browse/O2-2967)