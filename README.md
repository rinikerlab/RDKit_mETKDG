# RDKit_mETKDG

This Docker image provides a Debian distribution with Python 3 and all other dependencies installed via the Debian package manager. A modified version of RDKit from the master branch (as of April 2019) is compiled from source.

In particular, the modifications are for RDKit's experimental torsion knowledge distance geometry ([ETKDG](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00654)) conformer generator in order to improve sampling of macrocyclic molecules.

## Build Image
`docker build -t rdkit .`

## Run Container
To bring up a interactive terminal:
`docker run -it --rm rdkit:latest`

Or invoke the python 3 interpreter via:
`docker run -it --rm rdkit:latest python`

## Examples
Example python scripts are available in the [examples folder](./examples/) with README.


## Maintainer
Shuzhe Wang

## Acknowledgement
This docker image is created with https://github.com/mcs07/docker-rdkit as a template.
