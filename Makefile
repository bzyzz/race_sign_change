# Build the C programs and scan libraries.
# PRIMESIEVE = prefix of a primesieve (>= 12) install, with include/primesieve.h and lib/libprimesieve.so
# (e.g. a conda environment with `conda install -c conda-forge primesieve`).
CC        ?= gcc
PRIMESIEVE ?= $(CONDA_PREFIX)
PS_FLAGS   = -I$(PRIMESIEVE)/include -L$(PRIMESIEVE)/lib -Wl,-rpath,$(PRIMESIEVE)/lib -lprimesieve

.PHONY: all count primesieve scan

all: count scan primesieve

count:                       # own segmented sieves (no dependencies)
	$(CC) -O3 -march=native -fopenmp census/team_count.c -o census/team_count -lm
	$(CC) -O3 -march=native census/count_races.c -o census/count_races -lm

primesieve:                  # independent count and lead-region listing with primesieve
	$(CC) -O2 -fopenmp census/check_primesieve.c $(PS_FLAGS) -o census/check_primesieve
	$(CC) -O2 census/zones.c $(PS_FLAGS) -o census/zones

scan:                        # explicit-formula scan kernel used by the Python scripts
	$(CC) -O3 -fopenmp -shared -fPIC census/scan.c -o census/libscan.so -lm
	$(CC) -O3 -fopenmp -shared -fPIC certify/scan.c -o certify/libscan.so -lm
