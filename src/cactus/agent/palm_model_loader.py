""" Model loading scripts for loading Google Generative AI Models through LangChain"""

import os
from langchain.chat_models import ChatGooglePalm
from langchain.llms import GooglePalm


def load_google_model(model_name, *, temperature, api_key=None):
    """Loads a Google generative AI model into LangChain.

    Args:
      model_name: The name of the Google model to load.
      temperature: The hyperparameter that controls the randomness of the generated text.
      api_key: The Google PaLM API key. If not provided, the environment variable
        GOOGLE_API_KEY will be used.

    Returns:
      A LangChain ChatGooglePalm model object.
    """

    if api_key is None:
        try:
            if api_key is None:
                # Attempt to get the API key from the environment variables
                api_key = os.getenv("GOOGLE_API_KEY")

            if not api_key:
                raise ValueError("API key is not provided or found")

        except ValueError as error:
            print(f"Error: {error}")

    if model_name in ["models/chat-bison-001"]:
        llm = ChatGooglePalm(
            model_name=model_name,
            temperature=temperature,
            google_api_key=api_key,
        )

    else:
        llm = GooglePalm(
            model_name=model_name,
            temperature=temperature,
            google_api_key=api_key,
        )

    return llm
