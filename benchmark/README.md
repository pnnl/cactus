# Benchmarking Models

For this application we are benchmarking the following models:

| Model                 | model_name                                                      | Needs HF Token? | Notes                        |
|-----------------------|-----------------------------------------------------------------|-----------------|------------------------------|
| `opt13b`              | `facebook/opt-13b`                                              | N               |                              |
| `alpaca13b`           | `chavinlo/alpaca-13b`                                           | N               |                              |
| `galpacainstruct`     | `GeorgiaTechResearchInstitute/galactica-6.7b-evol-instruct-70k` | N               | An instruct fine-tuned model |
| `wizardvicuna13b`     | `ehartford/Wizard-Vicuna-13B-Uncensored`                        | N               | Uncensored Model             |
| `platypus30b`         | `lilloukas/Platypus-30B`                                        | N               |                              |
| `openorcaplatypus13b` | `Open-Orca/OpenOrca-Platypus2-13B`                              | N               |                              |
| `llama213b`           | `meta-llama/Llama-2-13b-hf`                                     | Y               |                              |
| `mixtral7b`           | `mistralai/Mixtral-8x7B-v0.1`                                   | N               | Mixture of Experts Model     |
| `gemma7b`             | `google/gemma-7b-it`                                            | Y               | An instruct fine-tuned model |

These models can be used by the `Cactus` agent by loading the model with the text in the first column (`Model`)

```python
from cactus.agent import Cactus

# Loading a model that does not require token
Model = Cactus(model_name="wizardvicuna13b")
```
For loading a model that requires a token:
```python
HF_TOKEN = "hf_..."

Model = Cactus(model_name="llama213b", token=HF_TOKEN)
```
Or if you set Token as env variable named: `$HF_Token` so `os.getenv("HF_Token")` returns the token

```python
Model = Cactus(model_name="llama213b")
```

