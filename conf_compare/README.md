# 2D4::conf_compare
A dislocation configuration comperator tool for the toolset 2D4.

## Short description
It compares configuration files and tell if they are the same, similar or different. Two configuration is the same if all dislocations' position both in *x* and *y*, and the Burger's vector are the same. Similar, if there is a little difference in *x*, different, if the difference is *x* is larger, and incompatible, if the number of dislocations or their *y* values are different.

## Dependencies
* boost-program-options

## Usage
The first two input argument of the program is the two files to compare. Example:
`conf_compare inputf_file1.dconf input_file2.dconf`

Further options:
* `--help`: shows a brief description of the parameters available
* `--hide-copyright` or `-c`: hides the copyright text
* `--sort arg` or `-s arg`: the program orders the dislocations before the comparison. The order is descending, first in Burger's vector then in *y* value so that *b * y* is descending. Default argument is `y`, telling that order in *y* is expected. Use any other, or `u` to use unsorted.
* `--individual-tolerance arg`: if all dislocations' *x* coordinate difference absolute value is smaller than this, the two configuration are considered the *same*.
* `--similarity-tolerance arg`: if the average difference absolute value of the *x* coordinates is smaller than this, the two configuration considered to be *similar*.