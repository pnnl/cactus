"""Agent File for Constructing Cactus."""

import platform

from langchain.agents import AgentExecutor
from langchain.agents.mrkl.base import ZeroShotAgent

from cactus.agent.anthropic_model_loader import load_anthropic_model
from cactus.agent.gemini_model_loader import load_google_model
from cactus.agent.huggingface_model_loaders import pipelines_model
from cactus.agent.openai_model_loader import load_openai_model
from cactus.agent.prompts import FORMAT_INSTRUCTIONS, PREFIX, SUFFIX
from cactus.agent.tools import make_tools

if platform.system() != "Darwin":
    from cactus.agent.vllm_model_loaders import vllm_model

__all__ = ["Cactus"]


def _load_model(
    model_name: str,
    model_type: str,
    cache_dir: str | None = None,
    max_length: int = 2000,
    *,
    use_8bit: bool = True,
    dtype=str,
):
    """Load a HuggingFace LLM.

    This function loads a HuggingFace LLM specified by its name using a
    registered model loader.

    Arguments
    ---------
        model_name (str): The name of the model to load.
            Default is "google/gemma-7b".
        model_type (str): The type of model being used.
            Default is "vllm".
        cache_dir (str, optional): The directory to cache the model data.
            Default is None.
        max_length (int, optional): The maximum allowed length for the model.
            Default is 2000.
        use_8bit (bool, optional): Whether to use 8-bit quantization for the
            model. Default is True.

    Returns
    -------
        Model: The loaded HuggingFace LLM.

    """
    if model_type == "api":
        if model_name in ["gpt-3.5-turbo", "gpt-4"]:
            return load_openai_model(model_name, temperature=0.7)

        if model_name in ["gemini-pro"]:
            return load_google_model(model_name, temperature=0.7)

        if model_name in [
            "claude-3-opus-20240229",
            "claude-3-sonnet-20240229",
            "claude-3-haiku-20240307",
        ]:
            return load_anthropic_model(model_name, temperature=0.7)

    elif model_type == "pipelines":
        return pipelines_model(
            model_name,
            cache_dir=cache_dir,
            max_length=max_length,
            use_8bit=use_8bit,
        )

    elif model_type == "vllm" and platform.system() != "Darwin":
        return vllm_model(
            model=model_name,
            cache_dir=cache_dir,
            dtype=dtype,
        )


class Cactus:  # pylint: disable=too-few-public-methods
    """LLM Agent for Answering Cheminformatics Questions."""

    def __init__(
        self,
        model_name="google/gemma-7b",
        model_type="vllm",
        cache_dir=None,
        max_length=2000,
        *,
        tools=None,
        use_8bit=True,
        dtype="float16"
    ):
        try:
            llm = _load_model(
                model_name,
                model_type,
                cache_dir,
                max_length,
                use_8bit=use_8bit,
                dtype=dtype,
            )

        except AssertionError as err:
            error_message = f"Module name: {model_name}, not found!"
            raise ValueError(error_message) from err

        if tools is None:
            tools = make_tools()
            self.tool_names = [tool.name for tool in tools]
            self.agent_executor = AgentExecutor.from_agent_and_tools(
                tools=tools,
                agent=ZeroShotAgent.from_llm_and_tools(
                    llm,
                    tools,
                    prefix=PREFIX,
                    format_instructions=FORMAT_INSTRUCTIONS,
                    suffix=SUFFIX,
                    input_variables=["input", "tool_names", "agent_scratchpad"],
                ),
                verbose=True,
                handle_parsing_errors=True,
            )

    def run(self, prompt):
        """Run the cactus agent over the question."""
        tool_names = self.tool_names
        outputs = self.agent_executor({"input": prompt, "tool_names": tool_names})
        return outputs["output"]
