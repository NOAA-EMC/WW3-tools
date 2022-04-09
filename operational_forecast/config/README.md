This folder contains files named `renames.conf` and `compress.conf` meant to be used by a rename and compress function
These files are simple text files that must have the following structure:

```
# lines with hashtag or that start with space are ignored
# NOTE: variable names cannot contain spaces

# RENAME FORMAT
old_name1 => new_name1
old_name2 => new_name2


#COMPRESS FORMAT
varname_after_renaming_1 => #n_digits_to_keep
varname_after_renaming_2 => #n_digits_to_keep
```
*spaces are ignored*

These files should be inside the folders with the same names as the scripts outputs i.e. the file `config/gdas/renames.conf` is only used to rename files inside the `gdas` prefix
