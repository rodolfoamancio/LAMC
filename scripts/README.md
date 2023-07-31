This folder contains additional scripts used for running the simulation on parallel or post-processing data.

# Scripts usage

# run_multiprocess.py
This script is used to run the simulation on multiple threads.

For using it, initially `cd` into the folder with the simulation `.inp` files, *e.g*:

```
example_1.inp
example_2.inp
...
example_n.inp
```
Then, export the path to the `LAMC` root folder:

```
$ export LAMC_DIR=path/to/LAMC
```
Finally, run the script:
```
$ python path/to/run_multiprocess.py
```

# get_results_summary.py

This script is used for post-processing properties results. As default if post-process each file which matches the pattern `*_properties_data.out.dat`. For instance, in a folder with:
```
example_1_properties_data.out.dat
example_2_properties_data.out.dat
...
example_3_properties_data.out.dat
```
The script can be executed with:
```
$ python path/to/get_results_summary.py
```
Optionally, the post-processing can be run on parallel informing the number of threads with:
```
$ python path/to/get_results_summary.py --cpus=<number_of_threads>
```
Finally, the results will be generated from the original name:
```
example_1_summary.csv
example_2_summary.csv
...
example_3_summary.csv
```

# combine_summaries.py

This script consolidates the result of each `summary.csv` result into a single table. It is necessary to run it on the folder where the files are stores, it will combine any file which matches the pattern `*summary.csv`. It requires a name for the output file. Usage:
```
$ python path/to/combine_summaries.py
```

# get_profiles results.py

This script post-process the raw `raw_profile.out.dat` files generated during the simulation. It requires the number of bins and number of threads. It should be executed in the folder where the referred `raw_profile.out.dat` are stored and will process any file which maches such pattern. Usage:

```
$ python path/to/get_profiles_results.py --number_bins=<number_bins> --cpus=<number_threads>
```
