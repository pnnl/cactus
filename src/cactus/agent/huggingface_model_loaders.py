"""Classes defining HuggingFace Model imports."""

from langchain_community.llms import HuggingFacePipeline
from transformers import (
    AutoModelForCausalLM,
    AutoTokenizer,
    pipeline,
)


def pipelines_model(
    model_id: str,
    cache_dir: str | None = None,
    max_length: int = 2000,
    use_8bit: bool = False,
) -> HuggingFacePipeline:
    """Create a text generation pipeline using the specified pre-trained HuggingFace model.

    Parameters
    ----------
        model_id: The identifier of the pre-trained model.
        cache_dir: The directory to cache the pre-trained model files.
        max_length: The maximum length of the generated text.
        use_8bit: Whether to load the model using 8-bit quantization.

    Returns
    -------
        A HuggingFacePipeline instance for text generation.
    """
    tokenizer = AutoTokenizer.from_pretrained(model_id, cache_dir=cache_dir)

    model = AutoModelForCausalLM.from_pretrained(
        model_id, load_in_8bit=use_8bit, device_map="auto", cache_dir=cache_dir
    )

    model.tie_weights()

    pipe = pipeline("text-generation", model=model, tokenizer=tokenizer, max_length=max_length)

    llm = HuggingFacePipeline(pipeline=pipe)
    llm.model_id = model_id

    return llm
