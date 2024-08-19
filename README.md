# HarmoDecon

### **Overview**

![Fig 1](https://github.com/raywzr/HarmoDecon/assets/81131673/bfc20f57-6258-4865-89b1-05ca1c9c8d84)

## **Description**

HarmoDecon is a deconvolution tool for spatially resolved transcriptomics data.

## **Installation**

Clone the project and create a virtual environment.

```bash
# clone the project
git clone https://github.com/raywzr/HarmoDecon.git
cd HarmoDecon

# create a conda environment
conda create -n HarmoDecon python=3.9.2
conda activate HarmoDecon
```

Install common package.

```bash
pip install --upgrade pip
pip install -r ./requirements.txt
```

##### For CPU implementation

Our model is lightweight. It can be run fast on a CPU.

```bash
pip install dgl==1.1.2
pip install torch==2.0.1+cpu torchvision==0.15.2+cpu
# or
pip install torch==2.0.1 torchvision==0.15.2 --index-url https://download.pytorch.org/whl/cpu
```

##### For CPU implementation

Look up the websites of DGL and Pytorch for command lines with specified CUDA settings.

https://www.dgl.ai/pages/start.html

https://pytorch.org/get-started/locally/

An example of CUDA 11.8:

```bash
pip install dgl==1.1.2 -f https://data.dgl.ai/wheels/cu118/repo.html
pip install torch==2.0.1 torchvision==0.15.2 --index-url https://download.pytorch.org/whl/cu118
```

##### Quick start with Docker

For convenience, we built a docker image of HarmoDecon (CPU version) to  implement our algorithm.

https://hub.docker.com/r/zrwangbrazy/harmodecon

```bash
docker pull zrwangbrazy/harmodecon
```

To run the docker:

```bash
# set the memory more than 4g. 
docker run --memory 5g harmodecon:v1
```

## **Input Format**

##### Single-cell RNA sequencing data: 

Should be in .h5ad format, with gene names annotated in "var", cell type annotated in "obs"(the key name should be "celltype").

![image-20240328153728086](https://github.com/raywzr/HarmoDecon/assets/81131673/20f44027-b28d-4a69-bb4b-337fae0f87fa)

![image-20240328153749220](https://github.com/raywzr/HarmoDecon/assets/81131673/e7ff0add-0896-429e-ada6-6cb24462cbbc)

##### Spatial transciptomics data: 

The gene Expression matrix should be in .h5ad format, with gene names annotated in "var".

![image-20240328153811950](https://github.com/raywzr/HarmoDecon/assets/81131673/bdb8bf7e-f00d-42d0-9810-c9392f82e4d6)

The location metadata should be in .txt or .tsv format, separated by '\t', including two columns: x and y.

![image-20240328145627548](https://github.com/raywzr/HarmoDecon/assets/81131673/db4c990d-bf18-4eb2-b5d1-a0fdf37007fe)

 Examples of data preprocessing(to get the above format) are included in ./script folder.

##### Pseudo spots data: 

In our method, scRNA-seq data are used to generate/synthesize pseudo spots. Pseudo-spot data includes cell type proportion profiles, which can be utilized to instruct the model.

![image-20240328153835689](https://github.com/raywzr/HarmoDecon/assets/81131673/289f0d60-0ec7-47a0-a484-eb8b8db6a2f2)

The input of the model can either be scRNA-seq data or pseudo spots. If the scRNA-seq is used, the model will synthesize spots from scratch. If pseudo spots are available, they can be reused for training to save time.

## **Example**

To start deconvolution on ST data. A JSON file with directories and hyperparameters is needed. For detailed information of setting hyperparameters, please refer to main.py.

To run the toy example on the osmFISH dataset[1, 2]:

```bash
python main.py --input ./configs/osm.json
```

To evaluate and visualize the result.

```bash
# repeat the input configs and specified the epoch and random seed of the model
python test.py --input ./configs/osm.json --epoch 20 --seed 1
```

![osm_L4_5 IT CTX](https://github.com/raywzr/HarmoDecon/assets/81131673/fae2dd0c-3f56-4934-9d36-971d8ae64483)

![osm_pie_plot](https://github.com/user-attachments/assets/415ba550-4b79-43e4-90a7-0ef0e270abec)


## **Reference**

[1] Codeluppi, S., Borm, L.E., Zeisel, A. *et al.* Spatial organization of the somatosensory cortex revealed by osmFISH. *Nat Methods* **15**, 932–935 (2018). https://doi.org/10.1038/s41592-018-0175-z

[2] Chen, J., Liu, W., Luo, T., Yu, Z., Jiang, M., Wen, J., Gupta, G. P., Giusti, P., Zhu, H., Yang, Y., & Li, Y. (2022). A comprehensive comparison on cell-type composition inference for spatial transcriptomics data. *Briefings in bioinformatics*, *23*(4), bbac245. https://doi.org/10.1093/bib/bbac245

## **Contact**

WANG Zirui: cszrwang@cse.hkbu.edu.hk zrwang18@gmail.com
