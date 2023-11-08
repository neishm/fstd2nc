[English version](README.md)

Aperçu
========
Ce module fournit un mécanisme de conversion entre les formats de fichiers FSTD et netCDF, via Python ou la ligne de commande.


Installation
==========

La façon plus simple d'installer est d'utiliser [pip](https://pip.pypa.io/en/stable):
```
pip install fstd2nc
```


Utilisation de base
===========

À partir de la ligne de commande
---------------------
```
python -m fstd2nc [options] <fichier(s) d'entrée> <fichier de sortie>

arguments optionnel:
  -h, --help            afficher ce message d'aide et quitter
  --version             afficher le numéro de version du programme et quitter
  --no-progress         Désactiver la barre de progression.
  --serial              Désactiver le multithreading/multitraitement. Utile
                        pour les machines à ressources limitées.
  --minimal-metadata    Ne pas inclure les attributs d'enregistrement internes
                        et d'autres informations internes dans les métadonnées
                        de sortie. Ceci est la comportement par défaut.
  --internal-metadata, --rpnstd-metadata
                        Inclure tous les attributs d'enregistrement internes
                        dans les métadonnées de sortie.
  --metadata-list nomvar,..., --rpnstd-metadata-list nomvar,...
                        Spécifier un ensemble minimal d'attributs
                        d'enregistrement internes à inclure dans le fichier de
                        sortie.
  --ignore-typvar       Indique au convertisseur d'ignorer la typvar lorsqu'il
                        décide si deux enregistrements font partie du même
                        champ. Le défaut est de diviser la variable sur
                        différents typvars.
  --ignore-etiket       Indique au convertisseur d'ignorer l'etiket lorsqu'il
                        décide si deux enregistrements font partie du même
                        champ. Le défaut est de diviser la variable sur
                        différents etikets.
  --vars VAR1,VAR2,...  Liste de variables séparées par des virgules à
                        convertir. Par défaut, toutes les variables sont
                        converties.
  --fill-value FILL_VALUE
                        La valeur de remplissage à utiliser pour les données
                        masquées (manquantes). Enregistré en tant qu'attribut
                        '_FillValue' dans les métadonnées. Le défaut est
                        '1e+30'.
  --datev, --squash-forecasts
                        Utiliser la date de validité pour l'axe "time". Ceci
                        est le défaut.
  --dateo, --forecast-axis
                        Utiliser la date d'analyse d'origine pour l'axe
                        temporel, et placer les heurs de prévision dans un axe
                        "forecast" séparé.
  --accum-vars NOM,NOM,...
                        Spécifier les champs à traiter comme des quantités
                        cumulées (en utilisant IP3 comme période de cumul).
  --ensembles           Rassembler différentes etikets pour la même variable
                        ensemble dans une axe "ensemble".
  --profile-momentum-vars VAR1,VAR2,...
                        Liste de variables séparées par des virgules qui
                        utilisent des niveaux thermodynamiques.
  --profile-thermodynamic-vars VAR1,VAR2,...
                        Liste de variables séparées par des virgules qui
                        utilisent des niveaux momentum.
  --missing-bottom-profile-level
                        Supposer que le niveau final des données de profil est
                        manquant.
  --vardict VARDICT     Utiliser les métadonnées du dictionnaire de variables
                        spécifié (format XML).
  --opdict              Similaire à ce qui précéde, mais utilise le
                        dictionnaire opérationnel standard CMC-RPN.
  --sfc-agg-vars NOM,NOM,...
                        Définir des champs d'agrégats de surface additionnel.
  --soil-depths SOIL_DEPTHS
                        Définir des profondeurs personnalisée pour les champs
                        de sol (WSOL,ISOL). Les défauts sont
                        0.05,0.1,0.2,0.4,1.0,2.0,3.0.
  --strict-vcoord-match
                        Exiger que les paramètres IP1/IP2/IP3 de la coordonnée
                        verticale correspondent aux paramètres IG1/IG2/IG3 du
                        champ afin d'être utilisés. La comportement par défaut
                        est d'utiliser l'enregistrement vertical de toute
                        façon s'il est le seul dans le fichier.
  --diag-as-model-level
                        Traiter les données diagnostic (près de la surface)
                        comme niveau de modèle '1.0'. Ceci est la comportement
                        par défaut.
  --split-diag-level    Placer les données diagnostic (près de la surface)
                        dans une variable séparée, loin de la sortie du modèle
                        3D.
  --ignore-diag-level   Ignorer les données sur la hauteur diagnostique (près
                        de la surface).
  --only-diag-level     Utiliser seulement la hauteur diagnostique (près de la
                        surface), en ignorant les autres niveaux
                        atmosphériques.
  --thermodynamic-levels, --tlev
                        Convertir seulement les données qui se trouvent sur
                        des niveaux verticaux 'thermodynamiques'.
  --momentum-levels, --mlev
                        Convertir seulement les données qui se trouvent sur
                        des niveaux verticaux 'momentum'.
  --vertical-velocity-levels, --wlev
                        Convertir seulement les données qui se trouvent sur
                        des niveaux 'vitesse verticale'.
  --subgrid-axis        Pour les données sur les 'supergrids', divisez les
                        grilles sur une axe "subgrid". Le défaut est garder
                        les sous-grilles ensemble comme elles sont dans le
                        fichier RPN.
  --keep-LA-LO          Inclure les enregstrement LA et LO, même s'ils
                        semblent redondants.
  --no-adjust-rlon      Pour les grilles tournées, n'ajustez PAS la coordonnée
                        rlon pour conserver la plage entre -180..180. Utilise
                        la plage que vient de librmn.
  --bounds              Inclure les limites des cellules de grille dans la
                        sortie.
  --filter CONDITION    Sous-ensemble d'enregistrements fichier standard RPN
                        en utilisant les critères donnés. Par exemple, pour
                        convertir seulement les prévisions sur 24 heures, vous
                        pouvez utiliser --filter ip2==24
  --exclude NOM,NOM,...
                        Exclure quelque axes, attributs, ou variables dérivées
                        de la sortie. Par exemple, l'exclusion de
                        'leadtime,reftime' peut aider les outils netCDF qui ne
                        reconnaissent pas leadtime et reftime comme des
                        coordonnées valides. Notez que les axes ne seront
                        exclus que s'ils ont une longueur de 1.
  --yin                 Sélectionner la première sous-grille d'une super
                        grille.
  --yang                Sélectionner la deuxième sous-grille d'une super
                        grille.
  --crop-to-smallest-grid
                        Couper les grilles au domaine plus petit (noyau
                        interne) pours les sorties LAM.
  --metadata-file METADATA_FILE
                        Utiliser les métadonnées du fichier spécifié.
  --rename OLDNAME=NEWNAME,...
                        Appliquer les changes de nom spécifiés aux variables.
  --conventions CONVENTIONS
                        Définir l'attribut "Conventions" pour le fichier
                        netCDF. Le défaut est "CF-1.6". Notez que cela n'a
                        aucun effet sur la structure du fichier.
  --no-conventions      Omettre entièrement l'attribut "Conventions" du
                        fichier netCDF.
  --time-units {seconds,minutes,hours,days}
                        Les unités de l'axe du temps de sortie. Le défaut est
                        hours.
  --reference-date AAAA-MM-JJ
                        La date de référence pour l'axe de temps de sortie. Le
                        défaut est la date de début dans le fichier standard
                        RPN.
  --fstd-compat         Ajoute une couche de compatibilité au fichier de
                        sortie netCDF, qu'il puisse également fonctionner
                        comme un fichier FSTD valide. EXPÉRIMENTAL.
  --msglvl {0,DEBUG,2,INFORM,4,WARNIN,6,ERRORS,8,FATALE,10,SYSTEM,CATAST}
                        Combien d'informations à imprimer sur stdout pendant
                        de la conversion. Le défaut est WARNIN
  --nc-format {NETCDF4,NETCDF4_CLASSIC,NETCDF3_CLASSIC,NETCDF3_64BIT_OFFSET,NETCDF3_64BIT_DATA}
                        Quelle variante de netCDF écrire. Le défaut est
                        NETCDF4.
  --zlib                Activer la compression pour le fichier netCDF. Ne
                        Fonctionne que pour les formats NETCDF4 et
                        NETCDF4_CLASSIC.
  --compression COMPRESSION
                        Niveau de compression du fichier netCDF. Ne utilisé
                        que si --zlib est specifié. Le défaut est 4.
  -f, --force           Remplacer le fichier de sortie s'il existe déjà.
  --no-history          Ne placer pas l'invocation de la ligne de commande
                        dans les métadonnées netCDF
  -q, --quiet           N'afficher pas quelque information moins les messages
                        d'erreur critiques. Implique --no-progress
```

Utilisation dans un script Python
========================

Conversion simple
--------------------------------------
```python
import fstd2nc
data = fstd2nc.Buffer("myfile.fst")
data.to_netcdf("myfile.nc")
```

Vous pouvez contrôler `fstd2nc.Buffer` en utilisant des paramètres similaires aux arguments de la ligne de commande.  La convention est que *--arg-name* de la ligne de commande serait passé comme *arg_name* de Python.

Par exemple:
```python
import fstd2nc
# Sélectionner uniquement les variables TT,HU.
data = fstd2nc.Buffer("myfile.fst", vars=['TT','HU'])
# Définir la date de référence au 1er Jan 2000 dans le fichier netCDF.
data.to_netcdf("myfile.nc", reference_date='2000-01-01')
```

Interfaçage avec xarray
---------------------------------------------------------------------------------

Pour les conversions plus compliquées, vous pouvez manipulez les données comme un objet [xarray.Dataset](http://xarray.pydata.org/en/stable/data-structures.html#dataset):
```python
import xarray as xr

# Ouvrir le fichier FSTD.
# Accéder les données comme un objet xarray.Dataset.
dataset = xr.open_dataset("myfile.fst", engine="fstd")
print (dataset)

# Convertir la pression de surface à Pa.
dataset['P0'] *= 100
dataset['P0'].attrs['units'] = 'Pa'

# (Pouvez encoure manipuler le dataset ici)
# ...

# Écrire le résultat final dans netCDF en utilisant xarray:
dataset.to_netcdf("myfile.nc")
```

Interfaçage avec iris
---------------------------------------------------------------------------------

Vous pouvez interfacer avec [iris](https://scitools.org.uk/iris/docs/latest/index.html) en utilisant la méthode `.to_iris()` (requiert iris au moins version 2.0).
Cela vous donnera un objet [iris.cube.CubeList](https://scitools.org.uk/iris/docs/latest/iris/iris/cube.html#iris.cube.CubeList):
```python
import fstd2nc
import iris.quickplot as qp
from matplotlib import pyplot as pl

# Ouvrir le fichier FSTD.
data = fstd2nc.Buffer("myfile.fst")

# Accéder les données comme un objet iris.cube.CubeList.
cubes = data.to_iris()
print (cubes)

# Tracer toutes les données (en supposant que nous avons des champs 2D)
for cube in cubes:
  qp.contourf(cube)
  pl.gca().coastlines()

pl.show()
```

Interfaçage avec pygeode
---------------------------------------------------------------------------------

Vous pouvez créer un objet [pygeode.Dataset](http://pygeode.github.io/dataset.html) en utilisant la méthode `.to_pygeode()` (requiert pygeode au moins version 1.2.2)
```python
import fstd2nc

# Ouvrir le fichier FSTD.
data = fstd2nc.Buffer("myfile.fst")

# Accéder les données comme un objet pygeode.Dataset.
dataset = data.to_pygeode()
print (dataset)
```

Interfaçage avec fstpy
---------------------------------------------------------------------------------

Vous pouvez charger des données à partir d'une table [fstpy](https://gitlab.science.gc.ca/CMDS/fstpy) en utilisant la méthode `.from_fstpy()` (requiert fstyp au moins version 2.1.9).
```python
import fstd2nc
import fstpy
table = fstpy.StandardFileReader('myfile.fst').to_pandas()
data = fstd2nc.Buffer.from_fstpy(table)
```
Vous pouvez aussi exporter vers une table fstpy en utilisant la méthode `.to_fstpy()`:
```python
import fstd2nc
table = fstd2nc.Buffer('myfile.fst').to_fstpy()
```


Exigences
============

Exigences de base
--------------------

Ce paquet requiert [Python-RPN](https://github.com/meteokid/python-rpn) pour lire/écrire des fichiers FSTD, et [netcdf4-python](https://github.com/Unidata/netcdf4-python) pour lire/écrire des fichiers netCDF.

Dépendances optionnelles
---------------------

Un dictionnaire de variables utile pour l'option `--vardict` et disponible [ici](https://collaboration.cmc.ec.gc.ca/cmc/CMOI/VariableDictionary/).

Pour lire un grand nombre de fichiers d'entrée (>100), cet utilitaire peut utiliser [pandas](https://github.com/pandas-dev/pandas) pour traiter rapidement l'en-tête d'enregistrement de FSTD.

La méthode Python `.to_xarray()` requiert les packages [xarray](https://github.com/pydata/xarray) et [dask](https://github.com/dask/dask).

La méthode Python `.to_iris()` requiert le paquet [iris](https://scitools.org.uk/iris/docs/latest/index.html), ainsi que les dépendances `.to_xarray()`.

La méthode Python `.to_pygeode()` requiert le paquet [pygeode](https://github.com/pygeode/pygeode), ainsi que les dépendances `.to_xarray()`.

La méthode Python `.to_fstpy()` requiert le paquet [fstpy](https://gitlab.science.gc.ca/CMDS/fstpy) (*lien interne*).
