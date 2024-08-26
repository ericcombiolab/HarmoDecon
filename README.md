# HarmoDecon

### **Overview**

![HarmoDecon 7 4](https://github.com/user-attachments/assets/48790181-b18d-47b9-ae4e-ca6191a665f3)


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
# Set the memory more than 4g. 
docker run -t -d --memory 5g harmodecon:v2
cd home
cd HarmoDecon
# Run the training script to reproduce the result of the osmFISH dataset.
python main.py --input ./configs/osm.json
# Plot pie charts and heat maps. Specify the random seed and which epoch you want to fetch the pre-trained model.
python test.py --input ./configs/osm.json --seed 1 --epoch 20
```

## **Input Format**

##### Single-cell RNA sequencing data: 

Should be in .h5ad format, with gene names annotated in "var", cell type annotated in "obs"(the key name should be "celltype").

![317642930-20f44027-b28d-4a69-bb4b-337fae0f87fa](https://github.com/ericcombiolab/HarmoDecon/assets/81131673/198188b3-f853-495a-aa35-1cabb716bd78)

![317642955-e7ff0add-0896-429e-ada6-6cb24462cbbc](https://github.com/ericcombiolab/HarmoDecon/assets/81131673/d0bea230-72e1-4b29-8625-0aaa60248f03)

##### Spatial transciptomics data: 

The gene Expression matrix should be in .h5ad format, with gene names annotated in "var".

![317643005-bdb8bf7e-f00d-42d0-9810-c9392f82e4d6](https://github.com/ericcombiolab/HarmoDecon/assets/81131673/fcf095e9-346f-4f76-9526-75f4b4f48c09)

The location metadata should be in .txt or .tsv format, separated by '\t', including two columns: x and y.

![317643095-db4c990d-bf18-4eb2-b5d1-a0fdf37007fe](https://github.com/ericcombiolab/HarmoDecon/assets/81131673/812025e6-aa3c-4d27-91ac-779c1db4a085)

 Examples of data preprocessing(to get the above format) are included in ./script folder.

##### Pseudo spots data: 

In our method, scRNA-seq data are used to generate/synthesize pseudo spots. Pseudo-spot data includes cell type proportion profiles, which can be utilized to instruct the model.

![317643154-289f0d60-0ec7-47a0-a484-eb8b8db6a2f2](https://github.com/ericcombiolab/HarmoDecon/assets/81131673/ff2e51e3-8cfa-40a9-b18a-0bc111ef30a5)

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

![317643376-fae2dd0c-3f56-4934-9d36-971d8ae64483](https://github.com/ericcombiolab/HarmoDecon/assets/81131673/af776b42-c9f1-4f98-b022-2cfe89934b05)

<<<<<<< HEAD
![osm_pie_plot](https://github.com/user-attachments/assets/415ba550-4b79-43e4-90a7-0ef0e270abec)
=======
![osm_pie_plot](https://github.com/user-attachments/assets/cc1cfb25-c857-4209-bce1-85f7e078b678)
>>>>>>> 0807dc96d0c81bc18218cb25b35cd08628e41954


## **Reference**

[1] Codeluppi, S., Borm, L.E., Zeisel, A. *et al.* Spatial organization of the somatosensory cortex revealed by osmFISH. *Nat Methods* **15**, 932–935 (2018). https://doi.org/10.1038/s41592-018-0175-z

[2] Chen, J., Liu, W., Luo, T., Yu, Z., Jiang, M., Wen, J., Gupta, G. P., Giusti, P., Zhu, H., Yang, Y., & Li, Y. (2022). A comprehensive comparison on cell-type composition inference for spatial transcriptomics data. *Briefings in bioinformatics*, *23*(4), bbac245. https://doi.org/10.1093/bib/bbac245

## **Contact**

If any questions, please feel free to contact me or my supervisor, Prof. Eric Zhang.

WANG Zirui: cszrwang@cse.hkbu.edu.hk zrwang18@gmail.com

Supervisor, Dr.Eric Zhang: ericluzhang@hkbu.edu.hk
