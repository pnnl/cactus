"""Agent File for Constructing Cactus"""
from typing import Optional

# from .prompts import ...
from .tools import make_tools
from .huggingface_model_loaders import HuggingFacePipelineFactory, pipeline_resolver
from .openai_model_loader import load_openai_model
from .palm_model_loader import load_google_model

from langchain.agents import AgentExecutor
from langchain.agents.mrkl.base import ZeroShotAgent

__all__ = ["Cactus"]


def _load_model(
    model_name: Optional[str] = None,
    cache_dir: Optional[str] = None,
    max_length: int = 2000,
    use_8bit: bool = True,
):
    """
    Load a HuggingFace LLM.

    This function loads a HuggingFace LLM specified by its name using a
    registered model loader.

    Args:
        model_name (str): The name of the model to load.
        cache_dir (str, optional): The directory to cache the model data.
            Default is None.
        max_length (int, optional): The maximum allowed length for the model.
            Default is 2000.
        use_8bit (bool, optional): Whether to use 8-bit quantization for the
            model. Default is True.

    Returns:
        Model: The loaded HuggingFace LLM.

    """
    if model_name in ["gpt-3.5-turbo", "gpt-4"]:
        return load_openai_model(model_name, temperature=0.7)

    if model_name in ["models/chat-bison-001", "models/text-bison-001"]:
        return load_google_model(model_name, temperature=0.7)

    model_loader: HuggingFacePipelineFactory = pipeline_resolver.make(
        model_name,
        cache_dir=cache_dir,
        max_length=max_length,
        use_8bit=use_8bit,
    )
    return model_loader.load_model()


class Cactus:  # pylint: disable=too-few-public-methods
    """LLM Agent for Answering Cheminformatics Questions"""

    def __init__(
        self,
        tools=None,
        model_name="alpaca13b",
        cache_dir=None,
        max_length=2000,
        use_8bit=True,
    ):
        try:
            llm = _load_model(
                model_name,
                cache_dir,
                max_length,
                use_8bit,
            )

        except AssertionError:
            print(f"Module name: {model_name}, not found!")

        if tools is None:
            tools = make_tools()
            self.agent_executor = AgentExecutor.from_agent_and_tools(
                tools=tools,
                agent=ZeroShotAgent.from_llm_and_tools(
                    llm,
                    tools,
                ),
                verbose=True,
            )

    def run(self, prompt):
        outputs = self.agent_executor({"input": prompt})
        return outputs["output"]
