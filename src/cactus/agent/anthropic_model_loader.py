"""Script to load Anthropic Models (Claude) into Cactus."""

import os

from langchain_anthropic import ChatAnthropic


def load_anthropic_model(model_name, *, temperature, api_key=None):
    """Load an Anthropic chat model into LangChain.

    Parameters
    ----------
      model_name: The name of the Anthropic chat model to load.
      temperature: The hyperparameter that controls the randomness of the generated text
      api_key: The Anthropic API key. If not provided, the environment variable
        ANTHROPIC_API_KEY will be used.

    Returns
    -------
      A LangChain ChatAnthropic model object.
    """
    if api_key is None:
        try:
            if api_key is None:
                # Attempt to get the API key from the environment variables
                api_key = os.getenv("ANTHROPIC_API_KEY")

            if not api_key:
                error_msg = "API key not provided or found"
                raise ValueError(error_msg)

        except ValueError as error:
            error_msg = f"Error: {error}"
            raise (error_msg) from error

    if model_name in [
        "claude-3-opus-20240229",
        "claude-3-sonnet-20240229",
        "claude-3-haiku-20240307",
    ]:
        llm = ChatAnthropic(
            model=model_name,
            temperature=temperature,
            api_key=api_key,
        )
    else:
        error_message = f"Model name: {model_name}, not found!"
        raise ValueError(error_message)

    return llm
