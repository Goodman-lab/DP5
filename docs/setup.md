# Setup

## Creating collateral files

### Input files

One should have MacroModel or TINKER and NWChem or Gaussian.

### Configuration files

The DP5 folder can contain `settings.cfg`, where the locations for the installed software packages
can be given, so that DP5 can launch them and use in the workflow.

The `settings.example` file provides a template for the configuration file
used by the tool. It should be renamed to `settings.cfg`, then its contents
changed as required. This can be done with the following command:

```bash
mv DP5/settings.example DP5/settings.cfg
$EDITOR DP5/settings.cfg
```

## Setting up clustered compute

To run calculations on a computational cluster, a password-less ssh connection
should be set up in both directions:

- Desktop -> Cluster
- Cluster -> Desktop

In most cases the modification of the relevant functions in `Gaussian.py`
or `NWChem.py` will be required to fit your situation.

## Utility scripts

The `scripts/` directory contains utility scripts are not part of
the tool, but are related to using it.

### `sdf2tinkerxyz`

This executable converts to and from TINKER nonstandard xyz files.

You may also want to add this to your `PATH` variable to use it
anywhere as well as the main tool.

### `Defaultslurm.sh`
 
This bash script helps run Slurm jobs on high performance compute.