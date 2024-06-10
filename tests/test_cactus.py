"""Test the Cactus Agent initialization."""

import pytest
from cactus.agent import Cactus


def test_cactus_init_default_model():
    """Tests Cactus class initialization with default model."""
    cactus = Cactus()
    assert cactus.agent_executor.agent.llm.model_name == "google/gemma-7b"
    assert cactus.agent_executor.agent.llm.model_type == "vllm"


def test_cactus_init_custom_model():
    """Tests Cactus class initialization with custom model."""
    custom_model_name = "google/gemma-2b"
    cactus = Cactus(model_name=custom_model_name)
    assert cactus.agent_executor.agent.llm.model_name == custom_model_name
    assert cactus.agent_executor.agent.llm.model_type == "vllm"


def test_cactus_init_invalid_model():
    """Tests Cactus class initialization with invalid model."""
    with pytest.raises(ValueError) as excinfo:
        Cactus(model_name="invalid_model_name")
    assert "Module name: invalid_model_name, not found!" in str(excinfo.value)


def test_cactus_run(mocker):
    """Tests Cactus class run method with mocked output."""
    prompt = "What is the molecular formula of water?"
    expected_output = "H2O"
    mocker.patch.object(Cactus.agent_executor, "run", return_value={"output": expected_output})
    cactus = Cactus()
    output = cactus.run(prompt)
    assert output == expected_output
