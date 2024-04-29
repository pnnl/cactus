"""Loads a vLLM hosted model for optimized inference."""

from langchain.llms import VLLM


def vllm_model(model, cache_dir, *, trust_remote_code=True, temperature=0.7):
    llm = VLLM(
        model=model,
        download_dir=cache_dir,
        trust_remote_code=trust_remote_code,  # For HF Models
        temperature=temperature,
    )

    return llm
