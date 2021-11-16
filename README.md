# dendritic_distance
Compute the distances between spines using spine ROIs and the average projection tiff


## Installation instructions

1. DD requires a  a recent Python distribution.

2. Open a terminal and type:

```bash
cd ..\documents\code #or similar as relevant for your machine
git clone git@github.com:GreggHeller1/dendritic_distance.git
cd phy
conda env create -f environment.yml
Conda activate dendritic_distance
pip install read_roi

```

3.  should now be installed and the dendritic distance environment should still be activated. 
4. to run the algorithm, change the path in config.py to a data directory containing a tiff and corresponding ROIs, or contianing nested directories with tiffs and ROIs
```bash
cd path/to/dendritic_distance_installation
python dendritic_distance.py
```

### Dependencies

For your information, here are the Python dependencies of dendritic_distance  (as found in `environment.yml`):

```

```

## Usage instructions
Change the path in 
