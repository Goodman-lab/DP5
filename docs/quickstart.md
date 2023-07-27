# Quickstart


!!! note WARNING

    **This tool should be run on a Linux or MacOS system, using a Python version
    $\geq$ 3.9 !**

    Using Windows or earlier versions of python will require tweaking -- do so at your
    own risk...


The following code can be copied and pasted into any bash terminal,
and should quickstart you running the tool.

```bash
# Download the source code from GitHub, and navigate into its directory
git clone https://github.com/Goodman-lab/DP5
cd DP5

# Run the make installation target.
# If make is not installed on your system, either install it or follow page on
# installation on this website
make install

# Activate the virtual environment containing the installed dependencies.
# Poetry is a python tool like pip which is installed by the make target if
# if isn't installed already
poetry shell

# Run the tool to see its manual. More instructions can then be found on the
# usage page on this website
pydp4 --help
```

More detailed instructions for installing, setting up, and running the tool can be
found on the [installation](installation.md), [setup](setup.md) and [usage](usage.md)
pages respectively.
