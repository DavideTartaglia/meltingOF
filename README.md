# OpenFOAM-nrg

A git repository of a paper I am trying to publish, in collaboration with NRG POLIMI (https://www.nuclearenergy.polimi.it). OpenFoam solver for fluids with different thermophysical properties. It is a temperature-based fixed-grid method. Developed and tested with [OpenFOAM-6](https://github.com/OpenFOAM/OpenFOAM-6). Compatible with OpenFOAM-nrg.

Steps for installation:

1.  Clone the repository to any installation directory outside any OpenFOAM project tree (e.g. `$HOME/OpenFOAM`):

    ```
    cd $HOME/OpenFOAM
    git clone <this-repository>
    ```

2.  Add the following line to the user's `$HOME/.bashrc` file:

    ```
    source $HOME/OpenFOAM/<this-repository>/etc/bashrc
    ```
    This assumes the installation directory to be `$HOME/OpenFOAM`. If not, modify the path accordingly.

3.  Compile the library by running the `Allwmake` script:

    ```
    cd <this-repository>
    ./Allwmake
    ```
    This will compile all solvers in `applications/solvers/` and libraries in `src/`.

