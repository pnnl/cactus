"""Script for running CACTUS on a set of benchmark questions."""

import argparse
import logging

import pandas as pd
from cactus.agent import Cactus


def main(logger, model_name, model_type, cache_dir, input_csv, output_csv):
    """Run the benchmark on a set of questions."""
    df = pd.read_csv(input_csv)

    column_data = df["Question"]
    logging.info(column_data)

    model = Cactus(
        model_name=model_name,
        model_type=model_type,
        cache_dir=cache_dir,
    )
    results = []

    for value in column_data:
        logger.info("Question: {value}")
        try:
            result = model.run(value)
        except Exception as e:
            result = f"Error: {str(e)}"

        results.append(result)
        logger.info("Result: {result}")

    df["result_column"] = results
    logging.info(df)
    df.to_csv(output_csv, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-name", "-m", required=True)
    parser.add_argument("--model-type", "-t", required=True)
    parser.add_argument("--cache-dir", required=False, default=None)
    parser.add_argument("--input-csv", "-i", required=True)
    parser.add_argument("--output-csv", "-o", required=True)
    parser.add_argument("--log-file", required=False, default="my_log.txt")
    args = parser.parse_args()

    logger = logging.getLogger("my_logger")
    logger.setLevel(logging.INFO)
    file_handler = logging.FileHandler(args.log_file)
    file_handler.setFormatter(logging.Formatter("%(asctime)s:%(name)s:%(levelname)s:%(message)s"))
    logger.addHandler(file_handler)

    main(
        logger,
        args.model_name,
        args.model_type,
        args.cache_dir,
        args.input_csv,
        args.output_csv,
    )
