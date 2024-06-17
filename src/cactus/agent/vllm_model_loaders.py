"""Loads a vLLM hosted model for optimized inference."""

from langchain.llms import VLLM


def vllm_model(model, cache_dir, *, trust_remote_code=True, temperature=0.7, dtype):
    """Load in a vllm compatible model.

    Parameters
    ----------
         model (str): Name of the model you want to load
         cache_dir (str): Path to the directory you want to save weights to
         trust_remote_code (Bool): Whether or not you want to trust remote code from HuggingFace
         temperature (Float): Hyperparameter dictating how creative the model should be when running inference

    Returns
    -------
         llm (BaseModel): A langchain compatible LLM ready for inference
    """
    llm = VLLM(
        model=model,
        download_dir=cache_dir,
        trust_remote_code=trust_remote_code,  # For HF Models
        temperature=temperature,
        dtype=dtype,
    )

    return llm
