import langchain
#import rmrkl
from .prompts import ...
from .tools import ... #TODO: Rename dir to tools

from .huggingface_model_loaders import HuggingFacePipelineFactory, pipeline_resolver

def _load_model(
    model_name: Optional[str] = None,
    cache_dir: str = None,
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
        model_name, cache_dir=cache_dir, max_length=max_length, use_8bit=use_8bit
    )
    return model_loader.load_model()

class CACTUS:
    def __init__(
            self,
            tools=None,
            model_name='alpaca13b',
            cache_dir=None,
            max_length=2000,
            use_8bit=True,
    ):
        try:
            llm = _load_model(model_name, cache_dir, max_length, use_8bit)

        except:
            assert #TODO Finish Assert

        if tools is None:
            pass
