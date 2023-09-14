![IPIML](https://github.com/SigProSeismology/IPIML/raw/main/IPIML_workflow.png)




#  IPIML: A deep-scan earthquake detection and location workflow Integrating Pair-Input deep learning model and  Migration Location method
IPIML is primarily developed and tested on Debian-based Linux OS systems. Therefore, we suggest using IPIML on such environments for the best experience. While it's possible to use IPIML on Windows and macOS, there may be challenges during compiling and running the workflow due to potential compatibility issues.

We greatly value community contributions and are steadfastly committed to continuously addressing and resolving any bugs that arise in the repository. Should you encounter any issues, please don't hesitate to contact us.

We implement the IPIML workflow in six steps, using an IPIML conda environments:

## Installation
The installation guides for these environments are provided below:

# IPIML envoronment:
Create and activate a conda environment, IPIML for detecting the primary events:
If you want to process with CPU:
```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n IPIMLL python=3.9 tensorflow==2.11.1 keras==2.11.0 h5py obspy spyder pygmt matplotlib pyyaml pandas tqdm pyproj jupyter notebook basemap six numpy protobuf
conda activate IPIMLL
pip install keras-rectified-adam seisbench
```
It you want to process with GPU:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n IPIML python=3.9 tensorflow-gpu==2.11.1 keras-gpu==2.11.0 h5py obspy spyder pygmt matplotlib pyyaml cudatoolkit cudnn pandas tqdm pyproj jupyter notebook basemap six numpy protobuf
conda activate IPIML
pip install keras-rectified-adam seisbench
```

### Install loki (GNU gcc compiler and openmp required)
```bash
git clone https://github.com/speedshi/LOKI.git
conda activate IPIML
cd WHERE_LOKI_IS_STORED
pip install .
```

### Install IPIML 
```bash
git clone https://github.com/SigProSeismology/IPIML.git](https://github.com/SigProSeismology/IPIML.git
```

## Usage 

Optimized Deep Learning (DL)-based workflows can improve the efficiency and accuracy of earthquake detection and location processes. IPIML is a six-step automated event detection, phase association, and earthquake location workflow, which integrates the state-of-the-art Pair-Input DL model and waveform Migration Location methods (IPIML). 

 

## To cite: 
Please cite the following paper in your documents if you use MALMI in your work. 

In case you utilize IPIML for processing your data, it would be appreciated if you cite the following paper(s):


BibTex:
```
@article{mohammadigheymasi2023ipiml,
  title={IPIML: A deep-scan earthquake detection and location workflow Integrating Pair-Input deep learning model and Migration Location method},
  author={Mohammadigheymasi, Hamzeh and Shi, Peidong and Tavakolizadeh, Nasrin and Xiao, Zhuowei and Mousavi, S. Mostafa and Matias, Luis and Pourvahab, Mehran and Fernandes, Rui},
  journal={IEEE Transactions on Geoscience and Remote Sensing},
  volume={XX},
  number={XX},
  pages={XX--XX},
  year={2023},
  doi={XX.XXXX/XXXXXXX}
}

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

@article{mousavi2020earthquake,
  title={Earthquake transformer—an attentive deep-learning model for simultaneous earthquake detection and phase picking},
  author={Mousavi, S Mostafa and Ellsworth, William L and Zhu, Weiqiang and Chuang, Lindsay Y and Beroza, Gregory C},
  journal={Nature communications},
  volume={11},
  number={1},
  pages={3952},
  year={2020},
  publisher={Nature Publishing Group UK London}
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



