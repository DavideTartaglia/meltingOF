# meltingOF

OpenFoam melting solver for fluids with different thermophysical properties from solid to liquid. It is a temperature-based fixed-grid method. Developed and tested with [OpenFOAM-6](https://github.com/OpenFOAM/OpenFOAM-6). Developed in collaboration with NRG @polimi (https://www.nuclearenergy.polimi.it/about-us/). Compatible with OpenFOAM-nrg repository.

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

## Reference :
https://www.tandfonline.com/doi/full/10.1080/00295639.2023.2229576
