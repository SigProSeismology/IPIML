![IPIML](https://github.com/SigProSeismology/IPIML/raw/main/IPIML_workflow.png)




#  IPIML: A deep-scan earthquake detection and location workflow Integrating Pair-Input deep learning model and  Migration Location method
IPIML is primarily developed and tested on Debian-based Linux OS systems. Therefore, we suggest using IPIML on such environments for the best experience. While it's possible to use IPIML on Windows and macOS, there may be challenges during compiling and running the workflow due to potential compatibility issues.

We greatly value community contributions and are steadfastly committed to continuously addressing and resolving any bugs that arise in the repository. Should you encounter any issues, please don't hesitate to contact us.

We implement the IPIML workflow in six steps, using two conda environments: ESR and MIL. Each environment serves a different purpose within the workflow, allowing for a clean and organized development process.

## Installation
The installation guides for these environments are provided below:

# ESR envoronment:
Create and activate a conda environment, ESR for detecting the primary events:
If you want to process with CPU:
```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n ESR python=3.7 tensorflow=1.14 keras=2.3.1 h5py=2.10 obspy spyder==5.0.3 pygmt matplotlib=3.2 pyyaml cudatoolkit cudnn pandas tqdm pyproj jupyter notebook basemap six~=1.15.0 numpy~=1.19.2 protobuf'<3.20,>=3.9.2'
pip install keras-rectified-adam
```
It you want to process with GPU:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n ESR python=3.7 tensorflow-gpu=1.14 keras-gpu=2.3.1 h5py=2.10 obspy spyder==5.0.3 pygmt matplotlib=3.2 pyyaml pandas tqdm pyproj jupyter notebook basemap six~=1.15.0 numpy~=1.19.2 protobuf'<3.20,>=3.9.2'
pip install keras-rectified-adam
```
# MIL envoronment:
Create and activate a conda environment, MIL for locating the primary events:
```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n MIL python=3.9 pygmt
conda activate MIL
pip install seisbench
```

### Install loki (GNU gcc compiler and openmp required)
```bash
git clone https://github.com/speedshi/LOKI.git
cd WHERE_LOKI_IS_STORED
pip install .
```

### Install MALMI 
```bash
git clone https://github.com/speedshi/MALMI.git
```

### Install NonLinLoc if you want to generate travetime tables in MALMI (optional)
Currently only *Vel2Grid* and *Grid2Time* programs are used and remember to put them in a executable path after compiling NonLinLoc. 
There are two ways to install NonLinLoc:

1. Through [NonLinLoc GitHub Page](https://github.com/alomax/NonLinLoc) (Recomended) 
Install example:
```bash
git clone https://github.com/alomax/NonLinLoc.git
cd NonLinLoc/src
mkdir bin   # bin/ is a subdirectory of src/
cmake .
make
echo 'export PATH="WHERE_CODE_IS_STORED/NonLinLoc/src/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

2. Follow [NonLinLoc Home Page](http://alomax.free.fr/nlloc/) for installing the NonLinLoc software. 
Install example:
```bash
wget http://alomax.free.fr/nlloc/soft7.00/tar/NLL7.00_src.tgz
mv NLL7.00_src.tgz WHERE_CODE_IS_STORED
cd WHERE_CODE_IS_STORED
mkdir NLL
tar -xzvf NLL7.00_src.tgz -C ./NLL
cd ./NLL/src
make -R distrib
echo 'export PATH="WHERE_CODE_IS_STORED/NLL/src:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

## Usage 
Follow the [user manual page](https://github.com/speedshi/MALMI/blob/main/user_manual.md) to use MALMI. 

## Reference 
Please cite the following paper in your documents if you use MALMI in your work. 
Peidong Shi, Francesco Grigoli, Federica Lanza, Gregory C. Beroza, Luca Scarabello, Stefan Wiemer; MALMI: An Automated Earthquake Detection and Location Workflow Based on Machine Learning and Waveform Migration. Seismological Research Letters 2022; doi: [https://doi.org/10.1785/0220220071](https://doi.org/10.1785/0220220071)

BibTex:
```
@article{xiao2021siamese,
  title={Siamese earthquake transformer: A pair-input deep-learning model for earthquake detection and phase picking on a seismic array},
  author={Xiao, Zhuowei and Wang, Jian and Liu, Chang and Li, Juan and Zhao, Liang and Yao, Zhenxing},
  journal={Journal of Geophysical Research: Solid Earth},
  volume={126},
  number={5},
  pages={e2020JB021444},
  year={2021},
  publisher={Wiley Online Library}
}
@article{10.1785/0220220071,
    author = {Shi, Peidong and Grigoli, Francesco and Lanza, Federica and Beroza, Gregory C. and Scarabello, Luca and Wiemer, Stefan},
    title = "{MALMI: An Automated Earthquake Detection and Location Workflow Based on Machine Learning and Waveform Migration}",
    journal = {Seismological Research Letters},
    year = {2022},
    month = {05},
    issn = {0895-0695},
    doi = {10.1785/0220220071},
    url = {https://doi.org/10.1785/0220220071},
    eprint = {https://pubs.geoscienceworld.org/ssa/srl/article-pdf/doi/10.1785/0220220071/5602568/srl-2022071.1.pdf},
}



```

## License 
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. For more details, see in the license file.

## Contributing
If you would like to contribute to the project or have any suggestions about the code, please feel free to create Pull Requests, raise issues and contact me. 
If you have any questions about the usage of this package or find bugs in the code, please also feel free to contact me.

## Contact information 
Copyright(C) 2023 Hamzeh Mohammadigheymasi 
Author: Hamzeh Mohammadigheymasi (hamzeh@ubi.pt), Peidong Shi (peidong.shi@sed.ethz.ch), and Xiao Zhuowei  (xiaozhuowei@mail.iggcas.ac.cn)
Email: , , 






