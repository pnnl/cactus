"""Agent File for Constructing Cactus"""
from typing import Optional

# from .prompts import ...
# from .tools import ...

from .huggingface_model_loaders import HuggingFacePipelineFactory, pipeline_resolver

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
            pass

    def run(self, prompt):
        pass
