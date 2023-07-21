# HighFluxSPL
Compensating for dead times in high-flux single-photon lidar 

The code included here was used to generate the figures in the paper
“High-flux single-photon lidar,” *Optica* 2021
[Joshua Rapp](joshuarapp.org), [Yanting Ma](https://www.merl.com/people/yma), Robin Dawson, [Vivek K Goyal](https://www.vivekgoyal.org/) 

contact: jrapp@bu.edu

## Getting started
**Code**: Clone the repository to your local machine.

**Data**: [Download the data](https://drive.google.com/file/d/1TwwtDn4BheJ9Ij25taWhk2_VLMedrtry/view?usp=sharing) and expand the zipfile into the root directory of the repo.

## Running the Code

The code is written in MATLAB and has been verified with MATLAB R2021b.

Each “main” file generates a different figure from the paper:

- main_compare_td_te.m generates Fig. 4, showing simulations that verify the Markov chain limiting distribution correctly predicts the detection time density;

- main_ranging_comparison.m generates Fig. 5, showing experimental results, validating the Markov chain modeling with real SPAD lidar data;

- main_histogram_correction.m generates Fig. 6, showing a comparison between our proposed Markov chain-based histogram correction (MCHC) method with the one of [Isbaner et al.](http://projects.gwdg.de/projects/deadtimecorrectiontcspc);

- main_depth_image_reconstruction.m generates Fig. 7, showing the point cloud results of a lidar imaging experiment with a mannequin.

## Citation
### Methods with detector and electronics dead times
```
@article{rapp-dead-optica-2021,
	title = {High-flux single-photon lidar},
	volume = {8},
	url = {https://doi.org/10.1364/OPTICA.403190},
	doi = {10.1364/optica.403190},
	number = {1},
	journal = {Optica},
	author = {Rapp, Joshua and Ma, Yanting and Dawson, Robin and Goyal, Vivek},
	month = jan,
	year = {2021},
	pages = {30--39}
}
```
### Methods with only detector dead times
```

@article{rapp-dead-tsp-2019,
	title = {Dead time compensation for high-flux ranging},
	url = {https://ieeexplore.ieee.org/document/8705308/},
	doi = {10.1109/TSP.2019.2914891},
	journal = {IEEE Transactions on Signal Processing},
	author = {Rapp, Joshua and Ma, Yanting and Dawson, Robin M. A. and Goyal, Vivek K.},
	month = jul,
	year = {2019},
	pages = {3471--3486},
}

```
