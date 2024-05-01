# CACTUS 🌵 | Chemistry Agent Connecting Tool Usage to Science
<img width="1000" alt="Cactus_header" src="workflow_diagram_V2__1_.png"> 

# Preprint 

## Available here: [link to preprint]

# Introduction 

CACTUS is an innovative tool-augmented language model designed to assist researchers and chemists in various chemistry-related tasks. By integrating state-of-the-art language models with a suite of powerful cheminformatics tools, CACTUS provides an intelligent and efficient solution for exploring chemical space, predicting molecular properties, and accelerating drug discovery.


## Running Cactus 🏃

Getting started with Cactus is as simple as:

```python
from cactus.agent import Cactus

Model = Cactus(model_name="google/gemma7b", model_type="vllm")
Model.run("What is the molecular weight of the smiles: OCC1OC(O)C(C(C1O)O)O")
```

## Installation 💻

First, we must install the version of `pytorch` that is useable by your system. We can do this by following the `torch` install instructions [here](https://pytorch.org/get-started/locally/).
Note that for some older versions of CUDA, you'll need to build from the compiled wheels [here](https://pytorch.org/get-started/previous-versions/)

For example, using `pip` to install PyTorch 2.0.1 with CUDA 11.7 from the pre-compiled wheel:

```bash
pip install torch==2.0.1+cu117 torchvision==0.15.2+cu117 torchaudio==2.0.2 --index-url https://download.pytorch.org/whl/cu117
```

Once you have the appropriate version of `torch`, we can continue with the install.

To install `cactus`, we can build from source:

```bash
git clone https://gitlab.pnnl.gov/computational_data_science/cactus.git
cd cactus
python -m pip install -e .
```

## Benchmarking 📊

We provide scripts for generating lists of benchmarking questions to evaluate the performance of the CACTUS agent.

These scripts are located in the `benchmark` directory.

To build the dataset used in the paper, we can run:

```bash
python benchmark_creation.py
```

This will generate a readable dataset named `QuestionsChem.csv` for use with the `Cactus` agent.

## Future Directions

We are continuously working on improving CACTUS and expanding its capabilities for molecular discovery. Some of our planned features include:

🧬 Integration with physics-based models for 3D structure prediction and analysis
🔧 Support for advanced machine learning techniques (e.g., graph neural networks)
🎯 Enhanced tools for target identification and virtual screening    
📜 Improved interpretability and explainability of the model's reasoning process

We welcome contributions from the community and are excited to collaborate with researchers and developers to further advance the field of AI-driven drug discovery.

## Citation 

If you use CACTUS in your research, please cite our preprint:
