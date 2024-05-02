# Benchmarking Models

In this benchmarking study, we evaluate the performance of CACTUS across a diverse set of state-of-the-art language models. By testing CACTUS with these models, we aim to assess its robustness, versatility, and effectiveness in solving chemistry-related tasks.

## Models Tested

For this application we are benchmarking the following models:

| Model        | model_name                  |
|--------------|-----------------------------|
| `llama2-7b`  | `meta-llama/Llama-2-7b-hf`  |
| `mistral-7b` | `mistralai/Mistral-7B-v0.1` |
| `gemma-7b`   | `google/gemma-7b-it`        |
| `falcon-7b`  | `tiiuae/falcon-7b`          |
| `MPT-7b`     | `mosaicml/mpt-7b`           |
| `Phi-2`      | `microsoft/phi-2`           |
| `OLMo-1b`    | `allenai/OLMo-1B`           |

These models were selected based on their strong performance in natural language tasks and their potential for adaptation to domain-specific applications.

## Benchmark Dataset

To evaluate the performance of CACTUS with each model, we have created a comprehensive benchmark dataset consisting of a wide range of chemistry-related questions. This dataset, named 'QuestionsChem.csv', can be generated using the 'benchmark_creation.py' script provided in this folder.

## Bechmarking Results

Comparison of model performance among 7B parameter models using minimal and domain-specific prompts. The Gemma-7b and Mistral-7b models demonstrate strong performance and adaptability across prompting strategies, highlighting their potential for widespread applicability in various computational settings, from high-performance clusters to more modest research setups.

<img width="1000" alt="Benchmarking Models" src="benchmark_files/Combined.png"> 

