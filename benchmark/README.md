# Benchmarking Models

In this benchmarking study, we evaluate the performance of CACTUS across a diverse set of state-of-the-art language models. By testing CACTUS with these models, we aim to assess its robustness, versatility, and effectiveness in solving chemistry-related tasks.

## Models Tested

For this application we are benchmarking the following models:

| Model        | model_name                         |
|--------------|------------------------------------|
| `llama2-7b`  | `meta-llama/Llama-2-7b-hf`         |
| `llama3-8b`  | `meta-llama/Meta-Llama-3-8B`       |
| `mistral-7b` | `mistralai/Mistral-7B-v0.1`        |
| `gemma-7b`   | `google/gemma-7b-it`               |
| `falcon-7b`  | `tiiuae/falcon-7b`                 |
| `MPT-7b`     | `mosaicml/mpt-7b`                  |
| `Phi-2`      | `microsoft/phi-2`                  |
| `Phi-3`      | `microsoft/Phi-3-mini-4k-instruct` |
| `OLMo-1b`    | `allenai/OLMo-1B`                  |

These models were selected based on their strong performance in natural language tasks and their potential for adaptation to domain-specific applications.

## Benchmark Dataset

To evaluate the performance of CACTUS with each model, we have created a comprehensive benchmark dataset consisting of a wide range of chemistry-related questions. This dataset can be generated using the `benchmark_creation.py` script provided in this folder. For the preprint we use the following benchmark datasets that are included here for reproducibility:

| File                                      | Description                                                                             |
|-------------------------------------------|-----------------------------------------------------------------------------------------|
| `SingleStepQuestionList_Qualitative.csv`  | The 500 __Qualitative__ questions used in the manuscript                                |
| `SingleStepQuestionList_Quantitative.csv` | The 500 __Quantitative__ questions used in the manuscript                               |
| `SingleStepQuestionList_Combined.csv`     | A concatenation of the Qualitative and Quantitative questions (1000 questions in total) |

The other files contained in this folder are:

| File                                  | Description                                                                                                                |
|---------------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| `benchmark_creation.py`               | The script used to generate the above benchmark question lists                                                             |
| `compound_list.csv`                   | A collection of 1000 molecules found on PubChem for benchmark creation                                                     |
| `run_benchmark.py`                    | A script to run a benchmark by using CACTUS                                                                                |
| `plot_creation.py`                    | A series of methods for generating the plots used in the manuscript                                                        |
| `plot_runtime.py`                     | A dedicated python script to generate the accuracy vs time plot in the manuscript                                          |
| `benchmark_files/`                    | A directory of all the resulting benchmark files used in the manuscript                                                    |
| `benchmark_files/Data_Analysis.ipynb` | A notebook detailing information about the dataset used, as well as how to calculate the expected answers from the dataset |

To run the benchmark using the `run_benchmark.py` script, you can do the following:

```shell
python run_benchmark.py --model_name "google/gemma-7b" --model-type "vllm" --input-csv "SingleStepQuestionList_Combined.csv" --output-csv "output.csv"
```

Optionally you can include a `--cache-dir` and a `--log-file` but these default to `None` and `my_log.txt` respectfully.

## Bechmarking Results

Comparison of model performance among 7B parameter models using minimal and domain-specific prompts. The Gemma-7b and Mistral-7b models demonstrate strong performance and adaptability across prompting strategies, highlighting their potential for widespread applicability in various computational settings, from high-performance clusters to more modest research setups.

<img width="1000" alt="Benchmarking Models" src="benchmark_files/Combined.png"> 

