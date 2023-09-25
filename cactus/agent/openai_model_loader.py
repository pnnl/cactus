"""Model loading scripts for loading Open AI models through LangChain"""

import os

from langchain.chat_models import ChatOpenAI
from langchain.llms import OpenAI


def load_openai_model(model_name, *, temperature, api_key=None):
    """Loads an OpenAI chat model into LangChain.

    Args:
      model_name: The name of the OpenAI chat model to load.
      temperature: The hyperparameter that controls the randomness of the generated text
      api_key: The OpenAI API key. If not provided, the environment variable
        OPENAI_API_KEY will be used.

    Returns:
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
            print(f"Error: {error}")

    if model_name in ["gpt-3.5-turbo", "gpt-4"]:
        llm = ChatOpenAI(
            model_name=model_name,
            temperature=temperature,
            openai_api_key=api_key,
        )

    else:
        llm = OpenAI(
            model_name=model_name,
            temperature=temperature,
            oppenai_api_key=api_key,
        )

    return llm
