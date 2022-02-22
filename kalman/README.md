# Test track models in Kalman fitter

Toy model for testing track models in a Kalman fitter for tracks propagating without a magnetic field, with of without Multiple Coulomb Scattering.

Presented by R. Pezzi at MFT Software and Physics meeting 2021, Aug. 31 [indico](https://indico.cern.ch/event/1071622/)

## Track parametrisations

- KF_01: the track parameters are (x, y, slopex, slopey)
- KF_02: the track parameters are (x, y, phi, tanl)

## How to use

```shell
$ root.exe
.x KF_01.C+
```

or

```shell
$ root.exe
.x KF_02.C+
```

## The histograms

The histograms show the pull distributions for the four parameters and the total chi-square distribution.
