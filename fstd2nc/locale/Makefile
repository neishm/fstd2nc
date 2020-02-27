# Makefile for handling translation files.
# See README.md for more information.

# Target language for the translation.
LANGUAGE=fr_CA
INSTALL_FILE=$(LANGUAGE)/LC_MESSAGES/fstd2nc.mo

# Rule for compiling the final (binary) translation file.
$(INSTALL_FILE): $(LANGUAGE).po
	mkdir -p $(dir $@)
	msgfmt $< -o $@

# Rule for updating a translation file when new strings are added to the
# program.
$(LANGUAGE).po: fstd2nc.pot
	[ -e $@ ] || msginit -i fstd2nc.pot -l $(LANGUAGE) --no-translator
	msgmerge --update $(LANGUAGE).po fstd2nc.pot
	touch $@

# Rule for generating a template file.
# Also include strings from argparse.
ARGPARSE_PYC = $(shell python -c "import argparse; print (argparse.__file__)")
ARGPARSE_PY = $(basename $(ARGPARSE_PYC)).py
fstd2nc.pot: ../*.py ../mixins/*.py
	xgettext $^ $(ARGPARSE_PY) -d fstd2nc --package-name=fstd2nc --no-location
	mv fstd2nc.po fstd2nc.pot

.DELETE_ON_ERROR: table

# Generate a table of translations in mediawiki format.
# Could be pasted into a wiki page, to solicit translations from users.
table: fr_CA.po
	echo '{| width="80%" cellspacing="1" cellpadding="10" border="1"' > $@
	echo '! style="width:50%;" |English' >> $@
	echo '! style="width:50%;" |Français' >> $@
	sed -E '1,/^$$/d;s/^msgid "+/|-\n|/;s/^msgstr /|/;/^#/d;s/""//g;s/^"//;s/"$$//;' $< >> $@
	echo '|}' >> $@

.PHONY: clean

# Rule for cleaning up the files.
# Removes the auto-generated template and binary files, but leaves the
# translation file (.po).
clean:
	rm -f fstd2nc.pot $(INSTALL_FILE) table

