# Onetep Parser

This is the parser for [Onetep](http://www.Onetep.org/).
It is part of the [NOMAD Laboratory](http://nomad-lab.eu).
The official version lives at

    git@gitlab.mpcdf.mpg.de:nomad-lab/parser-Onetep.git

you can browse it at

    https://gitlab.mpcdf.mpg.de/nomad-lab/parser-Onetep

It relies on having the nomad-meta-info and the python common repositories one level higher.
The simplest way to have this is to check out nomad-lab-base recursively:

    git clone --recursive git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-lab-base.git

then this will be in parsers/Onetep.


# TEST OUTPUT FILES 

Few output files where to test the parser are provided in the directory test examples.

        FILE NAME     |              FILE DESCRIPTION
    __________________|___________________________________________________
    "single_point.out" --> Single Point Calculation (minimum verbosity)
    