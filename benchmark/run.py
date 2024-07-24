from cactus.agent import Cactus
import pandas as pd
import os
import logging
import time


def main(logger):
    df = pd.read_csv("SingleStepQuestionList_Quantitative.csv")

    column_data = df["Question"]
    print(column_data)
    os.environ["GOOGLE_API_KEY"] = "AIzaSyA2L4C-nRNmApSNzFB7CtCiQyD8nISzpmI"
    Model = Cactus(model_name="gemini-pro", model_type="api")
    results = []

    for value in column_data:
        logger.info("Question: {value}")
        try:
            result = Model.run(value)
        except Exception as e:
            result = f"Error: {str(e)}"

        results.append(result)
        logger.info("Result: {result}")

    df["result_column"] = results
    print(df)
    print(os.getcwd())
    df.to_csv("gemini_prompt_quantitative_a100.csv", index=False)


if __name__ == "__main__":
    logger = logging.getLogger("my_logger")
    file_handler = logging.FileHandler("my_log.txt")
    file_handler.setFormatter(logging.Formatter("%(asctime)s:%(name)s:%(levelname)s:%(message)s"))
    logger.addHandler(file_handler)
    start = time.time()
    main(logger)
    end = time.time()
    print(end - start)
