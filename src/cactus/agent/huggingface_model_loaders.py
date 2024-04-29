"""Classes defining HuggingFace Model imports for CACTUS"""

import os

from class_resolver import Resolver
from langchain.llms import HuggingFacePipeline
from transformers import (
    AutoModelForCausalLM,
    AutoTokenizer,
    LlamaForCausalLM,
    LlamaTokenizer,
    pipeline,
)


class HuggingFacePipelineFactory:
    """
    A factory class for creating text generation pipelines using HuggingFace's
    transformers library.

    Args:
        model_id (str): The identifier of the pre-trained model to be used for
            text generation.
        cache_dir (str or None, optional): The directory to cache the
            pre-trained model files. Default is None.
        max_length (int, optional): The maximum length of the generated text.
            Default is 2000.
        use_8bit (bool, optional): Whether to load the model using 8-bit
            quantization. Default is True.

    Methods:
        load_model(model_id):
            This method must be implemented in the derived class.
            It loads a specific pre-trained model identified by the given
            'model_id'.

        create_pipeline(): Creates a text generation pipeline using the
            specified pre-trained model and tokenizer.

    Attributes:
        model_id (str): The identifier of the pre-trained model used for
            text generation.
        cache_dir (str or None): The directory to cache the pre-trained model
            files.
        max_length (int): The maximum length of the generated text.
        use_8bit (bool): Whether the model was loaded using 8-bit quantization.

    Example:
        factory = HuggingFacePipelineFactory(model_id='gpt2',
                                             cache_dir='/path/to/cache',
                                             max_length=1000)
        pipeline_instance = factory.create_pipeline()
    """

    def __init__(self, model_id, cache_dir=None, max_length=2000, *, use_8bit=True, token=None):
        self.model_id = model_id
        self.cache_dir = cache_dir
        self.max_length = max_length
        self.use_8bit = use_8bit
        self.token = token

    def load_model(self) -> HuggingFacePipeline:
        """
        This method must be implemented in the derived class.
        It loads a specific pre-trained model identified by the 'model_id'
        attribute of the instance.

        Raises:
            NotImplementedError: This method is meant to be overridden in the
                derived class.

        Returns:
            object: The loaded pre-trained model.
        """
        return self.create_pipeline()

    def create_pipeline(self) -> HuggingFacePipeline:
        """
        Creates a text generation pipeline using the specified pre-trained
        model and tokenizer.

        Returns:
            HuggingFacePipeline: A text generation pipeline using HuggingFace
                 transformers.
        """
        tokenizer = AutoTokenizer.from_pretrained(
            self.model_id,
            cache_dir=self.cache_dir,
            token=self.token,
        )

        model = AutoModelForCausalLM.from_pretrained(
            self.model_id,
            load_in_8bit=self.use_8bit,
            device_map="auto",
            cache_dir=self.cache_dir,
            token=self.token,
            attn_implementation="flash_attention_2",
        )

        model.tie_weights()

        pipe = pipeline(
            "text-generation", model=model, tokenizer=tokenizer, max_length=self.max_length
        )

        llm = HuggingFacePipeline(pipeline=pipe)
        llm.model_id = self.model_id

        return llm


class Opt13bLoader(HuggingFacePipelineFactory):
    """
    A derived class that represents a text generation pipeline for the
    'facebook/opt-13b' model.

    This class inherits from the HuggingFacePipelineFactory, which is a factory
    class for creating text generation pipelines using Hugging Face's
    transformers library.

    Args:
        cache_dir (str or None, optional): The directory to cache the
            pre-trained model files. Default is None.
        max_length (int, optional): The maximum length of the generated text.
            Default is 2000.
        use_8bit (bool, optional): Whether to load the model using 8-bit
            quantization. Default is True.

    Attributes:
        model_id (str): The identifier of the pre-trained model used for text
            generation ('facebook/opt-13b').
        cache_dir (str or None): The directory to cache the pre-trained model
            files.
        max_length (int): The maximum length of the generated text.
        use_8bit (bool): Whether the model was loaded using 8-bit quantization.

    Example:
        opt_loader = Opt13bLoader(cache_dir='/path/to/cache', max_length=1000)
        pipeline_instance = opt_loader.load_model()
    """

    def __init__(self, cache_dir=None, max_length=2000, *, use_8bit=True):
        model_id = "facebook/opt-13b"  # Set the model_id
        super().__init__(model_id, cache_dir=cache_dir, max_length=max_length, use_8bit=use_8bit)


class Alpaca13bLoader(HuggingFacePipelineFactory):
    """
    A derived class that represents a text generation pipeline for the
    'chavinlo/alpaca-13b' model.

    This class inherits from the HuggingFacePipelineFactory, which is a factory
    class for creating text generation pipelines using Hugging Face's
    transformers library.

    Args:
        cache_dir (str or None, optional): The directory to cache the
            pre-trained model files. Default is None.
        max_length (int, optional): The maximum length of the generated text.
            Default is 2000.
        use_8bit (bool, optional): Whether to load the model using 8-bit
            quantization. Default is True.

    Attributes:
        model_id (str): The identifier of the pre-trained model used for text
            generation ('chavinlo/alpaca-13b').
        cache_dir (str or None): The directory to cache the pre-trained model
            files.
        max_length (int): The maximum length of the generated text.
        use_8bit (bool): Whether the model was loaded using 8-bit quantization.

    Example:
        alpaca_loader = Alpaca13bLoader(cache_dir='/path/to/cache',
                                        max_length=1000)
        pipeline_instance = alpaca_loader.load_model()
    """

    def __init__(self, cache_dir=None, max_length=2000, *, use_8bit=True):
        model_id = "chavinlo/alpaca-13b"
        super().__init__(model_id, cache_dir=cache_dir, max_length=max_length, use_8bit=use_8bit)

    def create_pipeline(self) -> HuggingFacePipeline:
        tokenizer = LlamaTokenizer.from_pretrained(self.model_id, cache_dir=self.cache_dir)

        model = LlamaForCausalLM.from_pretrained(
            self.model_id, load_in_8bit=self.use_8bit, device_map="auto", cache_dir=self.cache_dir
        )

        model.tie_weights()

        pipe = pipeline(
            "text-generation", model=model, tokenizer=tokenizer, max_length=self.max_length
        )

        llm = HuggingFacePipeline(pipeline=pipe)
        llm.model_id = self.model_id

        return llm


class GalpacaInstructLoader(HuggingFacePipelineFactory):
    """
    A derived class that represents a text generation pipeline for the
    'GeorgiaTechResearchInstitute/galactica-6.7b-evol-instruct-70k' model.

    This class inherits from the HuggingFacePipelineFactory, which is a factory
    class for creating text generation pipelines using Hugging Face's
    transformers library.

    Args:
        cache_dir (str or None, optional): The directory to cache the
            pre-trained model files. Default is None.
        max_length (int, optional): The maximum length of the generated text.
            Default is 2000.
        use_8bit (bool, optional): Whether to load the model using 8-bit
            quantization. Default is True.

    Attributes:
        model_id (str): The identifier of the pre-trained model used for text
            generation
            ('GeorgiaTechResearchInstitute/galactica-6.7b-evol-instruct-70k').
        cache_dir (str or None): The directory to cache the pre-trained model
            files.
        max_length (int): The maximum length of the generated text.
        use_8bit (bool): Whether the model was loaded using 8-bit quantization.

    Example:
        galpaca_loader = GalpacaInstructLoader(cache_dir='/path/to/cache',
                                               max_length=1000)
        pipeline_instance = galpaca_loader.load_model()
    """

    def __init__(self, cache_dir=None, max_length=2000, *, use_8bit=True):
        model_id = "GeorgiaTechResearchInstitute/galactica-6.7b-evol-instruct-70k"
        super().__init__(model_id, cache_dir=cache_dir, max_length=max_length, use_8bit=use_8bit)


class WizardVicuna13bLoader(HuggingFacePipelineFactory):
    """
    A derived class that represents a text generation pipeline for the
    'ehartford/Wizard-Vicuna-13B-Uncensored' model.

    This class inherits from the HuggingFacePipelineFactory, which is a factory
    class for creating text generation pipelines using Hugging Face's
    transformers library.

    Args:
        cache_dir (str or None, optional): The directory to cache the
            pre-trained model files. Default is None.
        max_length (int, optional): The maximum length of the generated text.
            Default is 2000.
        use_8bit (bool, optional): Whether to load the model using 8-bit
            quantization. Default is True.

    Attributes:
        model_id (str): The identifier of the pre-trained model used for text
            generation ('ehartford/Wizard-Vicuna-13B-Uncensored').
        cache_dir (str or None): The directory to cache the pre-trained model
            files.
        max_length (int): The maximum length of the generated text.
        use_8bit (bool): Whether the model was loaded using 8-bit quantization.

    Example:
        wiz_vic_loader = WizardVicunaLoader(cache_dir='/path/to/cache',
                                            max_length=1000)
        pipeline_instance = wiz_vic_loader.load_model()
    """

    def __init__(self, cache_dir=None, max_length=2000, *, use_8bit=True):
        model_id = "ehartford/Wizard-Vicuna-13B-Uncensored"
        super().__init__(model_id, cache_dir=cache_dir, max_length=max_length, use_8bit=use_8bit)


class Platypus30bLoader(HuggingFacePipelineFactory):
    """
    A derived class that represents a text generation pipeline for the
    'lilloukas/Platypus-30B' model.

    This class inherits from the HuggingFacePipelineFactory, which is a factory
    class for creating text generation pipelines using Hugging Face's
    transformers library.

    Args:
        cache_dir (str or None, optional): The directory to cache the
            pre-trained model files. Default is None.
        max_length (int, optional): The maximum length of the generated text.
            Default is 2000.
        use_8bit (bool, optional): Whether to load the model using 8-bit
            quantization. Default is True.

    Attributes:
        model_id (str): The identifier of the pre-trained model used for text
            generation ('lilloukas/Platypus-30B').
        cache_dir (str or None): The directory to cache the pre-trained model
            files.
        max_length (int): The maximum length of the generated text.
        use_8bit (bool): Whether the model was loaded using 8-bit quantization.

    Example:
        platypus_loader = Platypus30bLoader(cache_dir='/path/to/cache',
                                            max_length=1000)
        pipeline_instance = platypus_loader.load_model()
    """

    def __init__(self, cache_dir=None, max_length=2000, *, use_8bit=True):
        model_id = "lilloukas/Platypus-30B"
        super().__init__(model_id, cache_dir=cache_dir, max_length=max_length, use_8bit=use_8bit)


class OpenOrcaPlatypus13BLoader(HuggingFacePipelineFactory):
    """
    A derived class that represents a text generation pipeline for the
    'Open-Orca/OpenOrca-Platypus2-13B' model.

    This class inherits from the HuggingFacePipelineFactory, which is a factory
    class for creating text generation pipelines using Hugging Face's
    transformers library.

    Args:
        cache_dir (str or None, optional): The directory to cache the
            pre-trained model files. Default is None.
        max_length (int, optional): The maximum length of the generated text.
            Default is 2000.
        use_8bit (bool, optional): Whether to load the model using 8-bit
            quantization. Default is True.

    Attributes:
        model_id (str): The identifier of the pre-trained model used for text
            generation ('Open-Orca/OpenOrca-Platypus2-13B').
        cache_dir (str or None): The directory to cache the pre-trained model
            files.
        max_length (int): The maximum length of the generated text.
        use_8bit (bool): Whether the model was loaded using 8-bit quantization.

    Example:
        orca_platypus_loader = OpenOrcaPlatypus13BLoader(cache_dir='/path/to/cache',
                                            max_length=1000)
        pipeline_instance = orca_platypus_loader.load_model()
    """

    def __init__(self, cache_dir=None, max_length=2000, *, use_8bit=True):
        model_id = "Open-Orca/OpenOrca-Platypus2-13B"
        super().__init__(model_id, cache_dir=cache_dir, max_length=max_length, use_8bit=use_8bit)


class Llama213BLoader(HuggingFacePipelineFactory):
    """Load the Llama2 Model"""

    def __init__(self, cache_dir=None, max_length=2000, *, use_8bit=True):
        model_id = "meta-llama/Llama-2-13b-hf"
        token = os.getenv("HF_Token")
        super().__init__(
            model_id, cache_dir=cache_dir, max_length=max_length, use_8bit=use_8bit, token=token
        )


class Mixtral7BLoader(HuggingFacePipelineFactory):
    """
    A derived class that represents a text generation pipeline for the
    'Open-Orca/OpenOrca-Platypus2-13B' model.

    This class inherits from the HuggingFacePipelineFactory, which is a factory
    class for creating text generation pipelines using Hugging Face's
    transformers library.

    Args:
        cache_dir (str or None, optional): The directory to cache the
            pre-trained model files. Default is None.
        max_length (int, optional): The maximum length of the generated text.
            Default is 2000.
        use_8bit (bool, optional): Whether to load the model using 8-bit
            quantization. Default is True.

    Attributes:
        model_id (str): The identifier of the pre-trained model used for text
            generation ('Open-Orca/OpenOrca-Platypus2-13B').
        cache_dir (str or None): The directory to cache the pre-trained model
            files.
        max_length (int): The maximum length of the generated text.
        use_8bit (bool): Whether the model was loaded using 8-bit quantization.

    Example:
        orca_platypus_loader = OpenOrcaPlatypus13BLoader(cache_dir='/path/to/cache',
                                            max_length=1000)
        pipeline_instance = orca_platypus_loader.load_model()
    """

    def __init__(self, cache_dir=None, max_length=2000, *, use_8bit=True):
        model_id = "mistralai/Mixtral-8x7B-v0.1"
        super().__init__(model_id, cache_dir=cache_dir, max_length=max_length, use_8bit=use_8bit)


class Gemma7BLoader(HuggingFacePipelineFactory):
    """
    A derived class that represents a text generation pipeline for the
    'Open-Orca/OpenOrca-Platypus2-13B' model.

    This class inherits from the HuggingFacePipelineFactory, which is a factory
    class for creating text generation pipelines using Hugging Face's
    transformers library.

    Args:
        cache_dir (str or None, optional): The directory to cache the
            pre-trained model files. Default is None.
        max_length (int, optional): The maximum length of the generated text.
            Default is 2000.
        use_8bit (bool, optional): Whether to load the model using 8-bit
            quantization. Default is True.

    Attributes:
        model_id (str): The identifier of the pre-trained model used for text
            generation ('Open-Orca/OpenOrca-Platypus2-13B').
        cache_dir (str or None): The directory to cache the pre-trained model
            files.
        max_length (int): The maximum length of the generated text.
        use_8bit (bool): Whether the model was loaded using 8-bit quantization.

    Example:
        orca_platypus_loader = OpenOrcaPlatypus13BLoader(cache_dir='/path/to/cache',
                                            max_length=1000)
        pipeline_instance = orca_platypus_loader.load_model()
    """

    def __init__(self, cache_dir=None, max_length=2000, *, use_8bit=True):
        model_id = "google/gemma-7b-it"
        token = os.getenv("HF_Token")
        super().__init__(
            model_id,
            cache_dir=cache_dir,
            max_length=max_length,
            use_8bit=use_8bit,
            token=token,
        )


pipeline_resolver = Resolver.from_subclasses(
    HuggingFacePipelineFactory, suffix="Loader", default=Alpaca13bLoader
)
