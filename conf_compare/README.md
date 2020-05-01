# 2D4::conf_compare
A dislocation configuration comperator tool for the toolset 2D4.

## Short description
It compares configuration files and tell if they are the same, similar or different. Two configuration
* are the same if all dislocations' position both in *x* and *y*, and the Burger's vector are the same,
* similar, if the difference in *x* is smaller than a threshold τ_s for all dislocations and and average difference is also smaller than another threshold τ_a, and 
* different, if the difference in *x* is larger than the threshold τ_s or the average difference is larger than threshold τ_a
* and incompatible, if the number of dislocations or their *y* values are different.

## Dependencies
* boost-program-options

## Usage
Both individual files and lists of files also can be compared. The first two (so-called positional) argument must be either 
1. the filenames of the individual dislocation configuration files (endings must be `.dconf`), or
2. the filenames of the files containing the list of filenames to compare (endings of the file lists must be `.ini` and the listed filenames' endings must be `.dconf`). Important program option is the [`find-to-compare`](#find_to_compare).

Example for the 1st case: the first two input arguments of the program are the two files to compare:
```
conf_compare inputf_file1.dconf input_file2.dconf
```

## Further options:
* `--help`: shows a brief description of the parameters available
* `find-to-compare`: <a label="find_to_compare"></a> if set, then if two filenames A and B are given to compare lists of files, then for each filename from file A will be compared with a file listed in file B such a way that the time difference between the two is the smallest along all the files listed in B. Otherwise, files will be compared in order, A_1 with B_1, A_2 with B_2, as long as A has elements. Exceed elemnts from B are neglected.
* `--hide-copyright` or `-c`: hides the copyright text
* `--sort arg` or `-s arg`: the program orders the dislocations before the comparison. The order is descending, first in Burger's vector then in *y* value so that *b * y* is descending. Default argument is `y`, telling that order in *y* is expected. Use any other, or `u` to use unsorted.
* `--individual-tolerance arg`: if all dislocations' *x* coordinate difference absolute value is smaller than this, the two configuration are considered the *same*.
* `--similarity-tolerance arg`: if the average difference absolute value of the *x* coordinates is smaller than this, the two configuration considered to be *similar*.