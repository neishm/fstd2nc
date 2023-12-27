# Conda packaging recipe   

- for each supported version of python, create a conda environement with the following python 3.X conda-verify conda-build boa anaconda-client   

```shell   
mamba create -n builder38 python 3.8 conda-verify conda-build boa anaconda-client   
```   

- activate the environement   
- from the fstd2nc directory run the following   
```shell
. activate builderXX 
conda mambabuild conda.recipe [-c <your private conda channel>]
```
- login to anaconda if you have a private channel
```shell
anaconda login
```
- upload your package   
```shell
anaconda upload <path to tar.bz2 file>
```


# Notes  
- removed pygeode dependency

