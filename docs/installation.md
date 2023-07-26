# Installation

!!! note

    This is a more in-depth explanation of the process run automatically by the [quickstart guide](quickstart.md).

    If you have already followed that, you don't need to run these again!

!!! warning

    **This tool should be run on a Linux or MacOS system, using a Python version
    $\geq$ 3.9 !**

    Using Windows or earlier versions of python will require tweaking -- do so at your
    own risk...


## Downloading the repository

The source code for the repository is hosted on GitHub. It can be downloaded with the
`git clone` command:

```bash
git clone https://github.com/Goodman-lab/DP5
cd DP5
```

## Installing the tool


### Installing poetry

DP5 uses a number of different Python libraries, which it manages with a tool
called `poetry`, which is slightly similar to `pip`. You can install it with
the following command:

```bash
pip install poetry
```

### Creating the virtual environment

In Python, a virtual environment is a place where dependencies can be
installed without interfering with the rest of the system. It is considered
best practice to use them, and can avoid a lot of issues when working on
multiple projects which use different (possibly mutually exclusive) libraries

To create our virtual environment, we can run the following two commands:

```bash
poetry install --without=qml
.venv/bin/pip install qml
```

??? note

    The observant reader might notice this is slightly more complicated than
    normal. The reason for omitting the QML library, then installing it with
    `pip` separately is that the library is packaged incorrectly, and has build
    time dependencies which are not satisfied in a fresh virtual environment.

    This is tracked in a [GitHub issue](https://github.com/qmlcode/qml/issues/145),
    so might be resolved in the future.

### Entering the virtual environment

Once you have created the virtual environment, you need to switch to use
it instead of your global installation of Python. This can be done with
the following command:

```bash
poetry shell
```

!!! warning

    Entering a virtual environment is tied to your terminal session,
    so this command must be re-run **for every new session**.

    If an error about dependencies not being installed occurs when
    running the tool, try re-running it to make sure it is activated
    correctly.

### Adding to the `PATH` variable


