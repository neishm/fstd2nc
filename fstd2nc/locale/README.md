How to contribute a translation
===============================

Any contribution is greatly appreciated.  The main file of interest is ```fr_CA.po```, which contains the current French translations.  If this file doesn't exist yet, you can create it by running
```
  make
```
This will create a new ```fr_CA.po``` file with all the original (English) strings, and stubs for inserting the French translations.

Once some translations are filled in or updated for this file, you can test them by re-running
```
  make
```
This will compile the translations into ```fstd2nc/locale/fr_CA/LC_MESSAGES/fstd2nc.mo```.  You should now see your translations when running the ```fstd2nc.py``` script from your local directory, with either ```LANGUAGE=fr_CA``` or ```CMCLNG=francais```, e.g.
```
export CMCLNG=francais
python -m fstd2nc --help
```

To contribute your translations back to the project, please create a commit with your ```fr_CA.po``` file, and submit a pull request.

