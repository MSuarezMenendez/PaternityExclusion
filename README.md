
# **Baleen whale microsatellite panel for individual identification and parentage assignment in *Mysticeti***

## Overview

Python scripts used on the manuscript "Baleen whale microsatellite panel for individual identification and parentage assignment in *Mysticeti*".

## Scripts

**Create_PO.py**: Outputs every combination of two individuals from a file.

Run example: ```python3 Create_PO.py```

**PO_check.py**: Parent-offspring exclusion.

Run example: ```python3 ./PO_check.py -m 30 -i Microsatellites.tsv -c PO_list.tsv```

+ **-m** Minimum number of markers that need to match to consider a PO pair
+ **-i** Input file with the microsatellite data (format example [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/Microsatellites.tsv))
+ **-c** Name of file with a list of PO pairs to check (format example [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/PO_list.tsv))


**Create_Calf_Mother.py**: Outputs every combination of calf-mother pairs, based on a list of samples and their sex (file names coded in script, examples [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/Samples.txt) and [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/Sexes.tsv)).

Run example: ```python3 Create_Calf_Mother.py```

**PaternityExclusion.py**: Paternity exclusion, it searches for a putative father for a given calf-mother pair.

Run example: ```python3 ./PaternityExclusion.py -m 30 -i Microsatellites.tsv -s Males.txt -c Calf_mother_list.tsv```

+ **-m** Minimum number of markers that need to match to consider a PO pair
+ **-i** Input file with the microsatellite data (format example [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/Microsatellites.tsv))
+ **-s** Name of file with a list of potential sires (format example [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/Males.txt))
+ **-c** Name of file with a list of Calf-Mother pairs to check (format example [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/Calf_mother_list.tsv))

**MultiLoci.py**: Increases the number of loci of a dataset by sampling with replacement from the input microsatellite dataset (microsatellite file input and number of loci to duplicate coded in the script).

Run example: ```python3 MultiLoci.py```

**Parentage_simulation.py**: Simulation of paternity assignment (Coded in script two input files with a list of females and males in the dataset, [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/Females.txt) and [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/Males.txt)).

Run example: ```python3 Parentage_simulation.py -r 1500 -i Extra_loci.tsv -o Simulation.pdf```

+ **-r** Number of replicates
+ **-i** Input file with the microsatellite data (output from MultiLoci.py, example [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/Extra_loci.tsv))
+ **-o** Output name for results figure

**Dataset_micro_sampling.py**: Creates subsampled datasets, starting from eight loci, adding an additional locus at a time till all loci are included. It creates 50 replicate datasets per set of loci (microsatellite file input coded in the script)

Run example: ```python3 Dataset_micro_sampling.py```

**Pat_Excl_Cal.py**: Calculates the combined non-exclusion probability of the first parent for the given dataset (input format example [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files/Microsatellites.tsv)).

Run example: ```python3 Pat_Excl_Cal.py Microsatellites.tsv```

**Examples of the input files for the scripts can be found** [**here**](https://github.com/MSuarezMenendez/PaternityExclusion/tree/main/Example_files)

## Preprint

M. Suárez-Menéndez, M. Bérubé, L. Bachmann, P. Best, M. P. Heide-Jørgensen, V. Lesage, T. Oosting, R. Prieto, C. Ramp, J. Robbins, R. Sears, M. A. Silva, M. Tollis, E. Vermeulen, G. A. Víkingsson, Ø. Wiig, and P. J. Palsbøll. Universal baleen whale microsatellite panel for individual identification and power to detect parentage (2023) bioRxiv 2023.04.12.536337; doi: https://doi.org/10.1101/2023.04.12.536337

## License

GNU General Public License v3.0

Copyright (C) 2023 Marcos Suarez

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (version 3 of the License)

This program is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose. See the
GNU General Public License for more details.

## Contact

Scripts were written by Marcos Suárez-Menéndez (marcos.sume@gmail.com) and Per Palsbøll
