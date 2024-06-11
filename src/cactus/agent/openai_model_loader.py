"""Model loading scripts for loading Open AI models through LangChain."""

import os

from langchain_openai import ChatOpenAI, OpenAI


def load_openai_model(model_name, *, temperature, api_key=None):
    """Load an OpenAI chat model into LangChain.

    Parameters
    ----------
      model_name: The name of the OpenAI chat model to load.
      temperature: The hyperparameter that controls the randomness of the generated text
      api_key: The OpenAI API key. If not provided, the environment variable
        OPENAI_API_KEY will be used.

    Returns
    -------
      A LangChain ChatOpenAI model object.
    """
    if api_key is None:
        try:
            if api_key is None:
                # Attempt to get the API key from the environment variables
                api_key = os.getenv("OPENAI_API_KEY")

            if not api_key:
                raise ValueError("API key not provided or found")

        except ValueError as error:
            raise error

    if model_name in ["gpt-3.5-turbo", "gpt-4"]:
        llm = ChatOpenAI(
            model=model_name,
            temperature=temperature,
            api_key=api_key,
        )

    else:
        llm = OpenAI(
            model=model_name,
            temperature=temperature,
            api_key=api_key,
        )

    return llm
