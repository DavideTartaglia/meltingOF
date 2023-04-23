# OpenFOAM-nrg

An [OpenFOAM](https://openfoam.org) library for simulation of nuclear reactors. Developed and tested with [OpenFOAM-6](https://github.com/OpenFOAM/OpenFOAM-6).

Steps for installation:

1.  Clone the repository to any installation directory outside any OpenFOAM project tree (e.g. `$HOME/OpenFOAM`):

    ```
    cd $HOME/OpenFOAM
    git clone <this-repository>
    ```

2.  Add the following line to the user's `$HOME/.bashrc` file:

    ```
    source $HOME/OpenFOAM/OpenFOAM-nrg-git/etc/bashrc
    ```
    This assumes the installation directory to be `$HOME/OpenFOAM`. If not, modify the path accordingly.

3.  Compile the library by running the `Allwmake` script:

    ```
    cd OpenFOAM-nrg-git
    ./Allwmake
    ```
    This will compile all solvers in `applications/solvers/` and libraries in `src/`.

