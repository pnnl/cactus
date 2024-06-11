"""Model loading scripts for loading Google Generative AI Models through LangChain."""

import os

from langchain_google_genai import ChatGoogleGenerativeAI, GoogleGenerativeAI


def load_google_model(model_name, *, temperature, api_key=None):
    """Load a Google generative AI model into LangChain.

    Parameters
    ----------
      model_name: The name of the Google model to load.
      temperature: The hyperparameter that controls the randomness of the generated text.
      api_key: The Google Gemini API key. If not provided, the environment variable
        GOOGLE_API_KEY will be used.

    Returns
    -------
      A LangChain ChatGoogleGenerativeAI model object.
    """
    if api_key is None:
        try:
            if api_key is None:
                # Attempt to get the API key from the environment variables
                api_key = os.getenv("GOOGLE_API_KEY")

            if not api_key:
                raise ValueError("API key is not provided or found")

        except ValueError as error:
            raise error

    if model_name in ["gemini-pro"]:
        llm = ChatGoogleGenerativeAI(
            model=model_name,
            temperature=temperature,
            google_api_key=api_key,
        )

    else:
        llm = GoogleGenerativeAI(
            model_name=model_name,
            temperature=temperature,
            google_api_key=api_key,
        )

    return llm
