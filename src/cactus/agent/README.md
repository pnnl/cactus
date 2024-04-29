# The CACTUS Agent

## How to Load Different Model Types

To use `Cactus`, you first must choose the source the model is being loaded from.

There are currently 3 main sources for LLMs in `Cactus`:
1. API
    - OpenAI API
    - Gemini API
2. HuggingFace Pipelines
    - Locally stored models downloaded from HuggingFace
    - Implemented through `Pipelines`
3. VLLM
    - Locally stored models downloaded from HuggingFace
    - Implemented through `VLLM`

When loading in `Cactus`, you will need to specify which model type you are trying to load:

```python
from cactus.agent import Cactus
api_model = Cactus(model_name="gpt-4", model_type="api")
pipelines_model = Cactus(model_name="wizardvicuna13b", model_type="pipelines")
vllm_model = Cactus(model_name="google/gemma7b", model_type="vllm")
```
### Loading Local Models from Cache

By default `Pipelines` and `VLLM` save tokenizers and model weights for the models you use to where the package is saved on your system. For example, if you use `conda`, the weights are saved to your `/envs` folder. This can be inconvenient if you want to save and load the models from a specific location (a location with more storage space for instance). To address this, `Pipelines` and `VLLM` allow for a path argument to tell the package where to look for the model. This is called the `cache_dir` in `Pipelines` and `download_dir` for `VLLM`. To make things convenient, we can simply pass a `cache_dir` argument to the `Cactus` agent and it will save/load the desired model to your specified path whether you are using `Pipelines` or `VLLM`.

```python
CACHE = "Path/To/Model/Weights"

pipelines_model = Cactus(model_name="gemma7b", model_type="pipelines", cache_dir=CACHE)
vllm_model = Cactus(model_name="mosaicml/mpt-7b", model_type="vllm", cache_dir=CACHE)
```

This will load in the tokenizer and model weights from the given path. 

🔔 Note: If this is the first time you are using a particular model, this command will download the tokenizer and model weights to the path first, and then load them in.

## How to Prompt Engineer

A key consideration to make when using open-source models hosted on HuggingFace is that the prompting strategy can be wildly different. We tried to consider a range of possible prompts that could be used with common models, but your mileage may vary. To address this, we provide a quick example on how to rapidly prompt engineer a new model using a python shell.

If, for example, we load in a VLLM model in the python REPL:
```python
from cactus.agent import Cactus

model = Cactus(model_name="meta-llama/Llama2", model_type="vllm")
```

This will load in some default prompt template comprised of a `PREFIX`, `FORMAT_INSTRUCTIONS`, and `SUFFIX` to describe how the LLM agent should proceed with answering questions. If you installed the project in editable mode (`python -m pip install -e .`) you can modify the `prompts.py` file with you model's specific prompt style. Or if you want to do it on the fly while running the code, you can access the current prompt from the following:
```python
model.agent_executor.agent.llm_chain.prompt.template

# To preserve the formatting, use:
print(model.agent_executor.agent.llm_chain.prompt.template)

# To make changes to the prompt just overwrite the resulting string
model.agent_executor.agent.llm_chain.prompt.template = """This is a new prompt!"""
```

This method of rapid prompt engineering can then be saved as a new `PREFIX`, `FORMAT_INSTRUCTIONS`, and `SUFFIX` to add to the prompts file to be loaded by `cactus.py`.
