#!/usr/bin/env python3
# pylint: disable=missing-module-docstring
#  See the NOTICE file distributed with this work for additional information
#  regarding copyright ownership.
#
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
"""Select assessed data from the database and create a report.
This script connects to a MySQL database, retrieves data related to sequencing runs,
applies quality checks using FastQC and STAR, and generates a report summarizing the results.
"""

import argparse
import re
import typing
from pathlib import Path
import pandas as pd
import pymysql
from sqlalchemy import create_engine


def connect_to_db(host:str, user:str, password:str, database:str, port=3306) \
    -> pymysql.connections.Connection:
    """
    Establishes a connection to the MySQL database.

    Args:
        host: The hostname or IP address of the MySQL server.
        user: The username to connect to the database.
        password: The password for the user.
        db: The name of the database to connect to.
        port: The port of the MySQL server (default is 3306).
    Returns:
        A pymysql connection object.
    """
    try:
        connection = pymysql.connect(host=host, user=user, password=password, database=database, port=port)
        return connection
    except pymysql.MySQLError as e:
        print(f"Error connecting to MySQL: {e}")
        return None


def fastqc_quality(row:pd.Series) -> bool:
    """Calculate FastQC quality based on the criteria.
    The function checks the FastQC quality criteria for each row and returns True if the criteria are met.


    Args:
        row (_type_): dataframe row

    Returns:
        bool: True if the criteria are met, False otherwise.
    """
    fastqc_criteria = [
        row["per_base_sequence_quality"] in ["PASS", "WARN"],
        row["overrepresented_sequences"] in ["PASS", "WARN"],
        row["per_base_n_content"] in ["PASS", "WARN"],
        row["per_sequence_quality_scores"] in ["PASS", "WARN"],
    ]
    return sum(fastqc_criteria) == 4


# Define STAR Quality classification
def star_quality(row:pd.Series) -> bool:
    """Calculate STAR quality based on the criteria.
    The function checks the STAR quality criteria for each row and returns True if the criteria are met.
    Args:
        row (_type_): dataframe row
    Returns:
        bool: True if the criteria are met, False otherwise.
    """
    return (
        row["uniquely_mapped_reads_percentage"] >= 50
        and row["percentage_reads_mapped_to_multiple_loci"] <= 30
        and row["percentage_reads_unmapped_too_short"] <= 20
    )


def check_fastqc_star_quality(df:pd.DataFrame) -> pd.DataFrame:
    """Apply FastQC and STAR quality to each row
    Args:
        df (pd.DataFrame): DataFrame containing the data to be processed.
    Returns:
        pd.DataFrame: DataFrame with additional columns for FastQC and STAR quality.
    """
    # Apply FastQC and STAR quality to each row
    df["fastqc_pass"] = df.apply(fastqc_quality, axis=1)
    df["star_quality"] = df.apply(star_quality, axis=1)

    # Assign final quality labels

    df["final_fastqc_status"] = df.groupby("run_accession")["fastqc_pass"].transform(
        lambda x: "Failed" if not all(x) else "Passed"
    )
    df["final_star_status"] = df.groupby("run_accession")["star_quality"].transform(
        lambda x: "Failed" if not all(x) else "Passed"
    )
    df["passed_both"] = (df["final_fastqc_status"] == "Passed") & (df["final_star_status"] == "Passed")
    # print(df.head(10))
    return df


def create_report(df:pd.DataFrame, tissue_report_file:str) -> pd.DataFrame:
    """Create a report based on the DataFrame.
    This function summarizes the FastQC and STAR quality results for each taxon_id and run_accession.
    It calculates the number of runs that passed both quality checks and the percentage of runs that passed.
    The report is printed to the console.
    1. Group the DataFrame by taxon_id and run_accession.
    2. Aggregate the final_fastqc_status and final_star_status columns.
    3. Create a new column passed_both to indicate if both quality checks passed.
    4. Group the DataFrame by taxon_id and aggregate the results.
    5. Calculate the percentage of runs that passed both quality checks.
    6. Print the report to the console.
    7. Return the summarized DataFrame.
    8. The report includes the following columns:
        - taxon_id: The taxon ID.
        - passed_fastqc: The number of runs that passed FastQC.
        - passed_star: The number of runs that passed STAR.
        - passed_both: The number of runs that passed both quality checks.
        - total_runs: The total number of runs for the taxon_id.
        - percentage_passed_both: The percentage of runs that passed both quality checks.
    Args:
        df (pd.DataFrame): DataFrame containing the data to be processed.
    Returns:
        pd.DataFrame: DataFrame with the summarized report.
    """
    # Collapse to one row per run_accession
    run_summary = (
        df.groupby(["taxon_id", "run_accession"])
        .agg(
            final_fastqc_status=("final_fastqc_status", "first"),
            final_star_status=("final_star_status", "first"),
        )
        .reset_index()
    )
    run_summary["passed_both"] = (run_summary["final_fastqc_status"] == "Passed") & (
        run_summary["final_star_status"] == "Passed"
    )
    # REPORT
    taxon_pass_summary = (
        run_summary.groupby("taxon_id")
        .agg(
            passed_fastqc=("final_fastqc_status", lambda x: (x == "Passed").sum()),
            passed_star=("final_star_status", lambda x: (x == "Passed").sum()),
            passed_both=("passed_both", "sum"),
            total_runs=("run_accession", "nunique"),
        )
        .reset_index()
    )

    # Compute percentage outside of .agg()
    taxon_pass_summary["percentage_passed_both"] = (
        taxon_pass_summary["passed_both"] / taxon_pass_summary["total_runs"] * 100
    ).round(2)
    taxon_pass_summary["passed_fastqc_only"] = (
        taxon_pass_summary["passed_fastqc"] - taxon_pass_summary["passed_both"]
    )
    taxon_pass_summary["passed_star_only"] = (
        taxon_pass_summary["passed_star"] - taxon_pass_summary["passed_both"]
    )

    # print(taxon_pass_summary.head())
    # print("\n=== Run Quality Summary by Taxon ===")
    for _, row in taxon_pass_summary.iterrows():
        print(
            f"Taxon ID: {int(row['taxon_id'])} | "
            f"Total Runs: {row['total_runs']} | "
            f"Passed BOTH: {row['passed_both']} | "
            # f"FastQC Passed: {row['passed_fastqc']} | "
            # f"STAR Passed: {row['passed_star']} | "
            f"🧪 FastQC Only: {row['passed_fastqc_only']} | "
            f"🚀 STAR Only: {row['passed_star_only']} | "
            f"✔️ Both %: {row['percentage_passed_both']}%"
        )
        # Prepare tissue-level summary
    if "tissue_prediction" in df.columns:
        run_tissue_summary = (
            df.groupby(["taxon_id", "tissue_prediction", "run_accession"])
            .agg(
                final_fastqc_status=("final_fastqc_status", "first"),
                final_star_status=("final_star_status", "first"),
            )
            .reset_index()
        )

        run_tissue_summary["passed_both"] = (run_tissue_summary["final_fastqc_status"] == "Passed") & (
            run_tissue_summary["final_star_status"] == "Passed"
        )

        tissue_pass_summary = (
            run_tissue_summary.groupby(["taxon_id", "tissue_prediction"])
            .agg(
                passed_fastqc=("final_fastqc_status", lambda x: (x == "Passed").sum()),
                passed_star=("final_star_status", lambda x: (x == "Passed").sum()),
                passed_both=("passed_both", "sum"),
                total_runs=("run_accession", "nunique"),
            )
            .reset_index()
        )

        tissue_pass_summary["percentage_passed_both"] = (
            tissue_pass_summary["passed_both"] / tissue_pass_summary["total_runs"] * 100
        ).round(2)

        tissue_pass_summary["passed_fastqc_only"] = (
            tissue_pass_summary["passed_fastqc"] - tissue_pass_summary["passed_both"]
        )
        tissue_pass_summary["passed_star_only"] = (
            tissue_pass_summary["passed_star"] - tissue_pass_summary["passed_both"]
        )

        # Save to text file
        with open(tissue_report_file, "w") as f:
            for _, row in tissue_pass_summary.iterrows():
                f.write(
                    f"Taxon ID: {int(row['taxon_id'])} | "
                    f"Tissue: {row['tissue_prediction']} | "
                    f"Total Runs: {row['total_runs']} | "
                    f"Passed BOTH: {row['passed_both']} | "
                    # f"FastQC Passed: {row['passed_fastqc']} | "
                    # f"STAR Passed: {row['passed_star']} | "
                    f"FastQC Only: {row['passed_fastqc_only']} | "
                    f"STAR Only: {row['passed_star_only']} | "
                    f"BOTH %: {row['percentage_passed_both']}%\n"
                )

    return df


def prioritise_tissues(df: pd.DataFrame, priority_tissues: typing.List[str]) -> pd.DataFrame:
    """
    Assign a numeric priority to tissue predictions based on a predefined list.

    This function adds a 'priority' column to the input DataFrame:
    - Rows where the 'tissue_prediction' matches any of the `priority_tissues`
    (case-insensitive, partial match) are assigned a priority of 1.
    - All other rows receive a priority of 2.

    Args:
        df (pd.DataFrame): DataFrame containing a 'tissue_prediction' column.
        priority_tissues (List[str]): List of tissue names to prioritize (e.g., ['brain', 'liver']).

    Returns:
        pd.DataFrame: Modified DataFrame with an added 'priority' column.
    """
    priority_regex = re.compile(r"|".join(map(re.escape, priority_tissues)), re.IGNORECASE)
    df["priority"] = df["tissue_prediction"].apply(
        lambda x: 1 if isinstance(x, str) and priority_regex.search(x) else 2
    )
    return df


def filter_data(df:pd.DataFrame) -> list:
    """
    Filter the DataFrame to select tissue-annotated samples based on quality criteria.

    Steps:
    1. Remove duplicate run_accessions.
    2. Filter out rows with null tissue predictions.
    3. Exclude rows with tumor/cancer-related tissue labels.
    4. Prioritize specific tissue types (e.g., brain, heart, lung).
    5. Group by tissue and select samples up to 250M total reads.
    6. If fewer than 500 samples selected, add more high-quality unannotated samples.

    Args:
        df (pd.DataFrame): Input DataFrame with tissue predictions and QC info.
        debug (bool): If True, write intermediate CSVs for debugging.

    Returns:
        List[str]: List of selected run_accessions.
    """
    df = df.drop_duplicates(subset="run_accession", keep="first")
    df_tissue = df[df["tissue_prediction"].notnull()].reset_index(drop=True)
    df_tissue = df_tissue[
        ~df_tissue["tissue_prediction"].str.contains(
            r"tumor|cancer|melanoma|cancer/disease", flags=re.IGNORECASE, regex=True
        )
    ].reset_index(drop=True)
    # Prioritise tissue types
    priority_tissues = ["heart", "lung", "brain", "ovary", "ovaries", "testes", "testis", "gonad", "gonads"]
    df_tissue = prioritise_tissues(df_tissue, priority_tissues)
    df_tissue.to_csv("tissueall.csv", index=False)
    df_tissue = df_tissue[
        ~(
            (~df_tissue["passed_both"])
            & (df_tissue["final_fastqc_status"] == "Failed")
            & (df_tissue["final_star_status"] == "Failed")
        )
    ]
    df_tissue.to_csv("tissueless.csv", index=False)
    df_no_tissue = df[df["tissue_prediction"].isna()].reset_index(
        drop=True
    )
    run_accessions = []
    # Filter tissue-annotated samples
    for _, group in df_tissue.groupby("tissue_prediction"):
        group = group.sort_values(
            by=["priority", "passed_both", "final_fastqc_status", "final_star_status"], ascending=False
        )
        group.to_csv("group.csv", index=False)
        total_length = 0
        for _, row in group.iterrows():
            if total_length + row["total_sequences"] <= 250 * 10**6:
                run_accessions.append(row["run_accession"])
                total_length += row["total_sequences"]

    df_tissue[df_tissue["run_accession"].isin(run_accessions)].to_csv("tissuedanio.csv", index=False)

    if len(run_accessions) < 500:
        for _, group in df_no_tissue.groupby("run_accession"):
            total_length = 0
            if group["passed_both"].any():
                for _, row in group[group["passed_both"]].iterrows():
                    run_accessions.append(row["run_accession"])

    return run_accessions


def clean_repeated_words(text:str)-> str:
    """Remove consecutive duplicate words from a string."

    Args:
        text (str): The input string from which to remove 
        consecutive duplicate words.

    Returns:
        str: A string with consecutive duplicate words removed,
        preserving the first occurrence of each word.
    """
    words = text.split()
    if not words:
        return ""
    # Keep only the first occurrence, skip consecutive duplicates
    cleaned_words = [words[0]]
    for w in words[1:]:
        if w != cleaned_words[-1]:
            cleaned_words.append(w)
    return " ".join(cleaned_words)


def main() -> None:
    """Module's entry-point."""
    parser = argparse.ArgumentParser(prog="llm_prediction.py", description="Predict tissue using LLMs")

    parser.add_argument("--taxon_id", default="10116", type=str, required=True, help="Taxonomy ID")
    parser.add_argument("--host", type=str, default="mysql-ens-genebuild-prod-1", required=False, help="Host")
    parser.add_argument("--user", type=str, default="ensadmin", required=False, help="User")
    parser.add_argument("--password", type=str, default="ensembl", required=False, help="Password")
    parser.add_argument(
        "--database",
        default="gb_transcriptomic_registry",
        type=str,
        required=False,
        help="Database",
    )
    parser.add_argument("--port", default=4527, type=int, required=False, help="Port")
    parser.add_argument("--file_name", type=str, required=False, help="Output file name")
    parser.add_argument(
        "--read_type",
        type=str,
        choices=["short", "long"],
        default="short",
        required=False,
        help="Read type short/long",
    )
    parser.add_argument(
        "--csv_for_main",
        action="store_true",
        required=False,
        help="if true will produce the csv we need for main pipeline otherwise \
            the one for the full alignment pipeline",
    )
    args = parser.parse_args()
    engine = create_engine(
        f"mysql+pymysql://{args.user}:{args.password}@{args.host}:{args.port}/{args.database}"
    )

    # Load the rows to fix (e.g., 150k rows)
    query = ""
    if args.read_type == "short":
        query = (
            "SELECT r.taxon_id, r.run_accession, r.sample_accession, r.platform, r.paired, \
            r.tissue_prediction, d.file_name, d.file_url,d.md5,d.per_base_sequence_quality, \
            d.per_sequence_quality_scores, d.per_base_n_content, d.overrepresented_sequences,\
            d.total_sequences, a.uniquely_mapped_reads_percentage, \
            a.percentage_reads_mapped_to_multiple_loci,a.percentage_reads_unmapped_too_short, \
            a.assembly_accession \
            FROM run r INNER JOIN data_files d ON r.run_id=d.run_id INNER JOIN align a ON r.run_id=a.run_id \
            WHERE  r.qc_status='ALIGNED' AND r.paired=1 AND r.platform='ILLUMINA' AND r.taxon_id = "
            + args.taxon_id
        )
    elif args.read_type == "long":
        query = (
            "SELECT r.taxon_id, r.run_accession, r.sample_accession, r.platform, r.paired, \
            r.tissue_prediction, d.file_name, d.file_url,d.md5,d.per_base_sequence_quality, \
            d.per_sequence_quality_scores, d.per_base_n_content, d.overrepresented_sequences, \
            d.total_sequences, a.uniquely_mapped_reads_percentage, \
            a.percentage_reads_mapped_to_multiple_loci,a.percentage_reads_unmapped_too_short \
            FROM run r INNER JOIN data_files d ON r.run_id=d.run_id INNER JOIN align a ON r.run_id=a.run_id \
            WHERE  r.qc_status='ALIGNED' AND r.paired=0 AND r.platform='PACBIO_SMRT' AND r.taxon_id = "
            + args.taxon_id
        )
    # Connect to the database
    # db_connection = connect_to_db(**db_config)
    # Clean the input text
    try:
        # if db_connection:
        # df = pd.read_sql(query, db_connection)
        df = pd.read_sql(query, engine)
        # print(f"Loaded {len(df)} rows.")
    except pymysql.MySQLError as e:
        print(f"Error connecting to MySQL: {e}")

    if not df.empty:
        df = check_fastqc_star_quality(df)

        #df["tissue_prediction"] = df["tissue_prediction"].apply(
        #    lambda x: re.search(r"^(.*)", str(x)).group(1) if re.search(r"^(.*)", str(x)) else x
        #)
        df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: str(x) if pd.notnull(x) else x
        )
        # df["tissue_prediction"] = df["tissue_prediction"].apply(
        #    lambda x: (m.group(1) if (m := re.search(r"^(.*)", str(x))) else x)
        # )
        df["tissue_prediction"] = (
            df["tissue_prediction"]
            .astype(str)
            .str.replace('"', "", regex=True)
            .replace("Answer: ", "", regex=True)
            .replace("\n", "", regex=True)
            .replace("Output: ", "", regex=True)
        )
        df["tissue_prediction"] = df["tissue_prediction"].apply(
            lambda x: re.sub(
                r"[, ]+",
                "_",  # replace space/comma with underscore
                re.sub(
                    r"\b\d+\b|\bday\b|\b[a-zA-Z]\b",
                    "",  # remove numbers, "day", and single letters
                    clean_repeated_words(str(x)),
                ),
            ).strip(
                "_"
            )  # remove leading/trailing underscores
        )
        df["tissue_prediction"] = df["tissue_prediction"].astype(str).str.lower()

        df.loc[
            df["tissue_prediction"].astype(str).str.contains("NONE", case=False, na=False), "tissue_prediction"
        ] = None
        df.loc[
            df["tissue_prediction"].astype(str).str.contains("nan", case=False, na=False), "tissue_prediction"
        ] = None

        create_report(df, f"{Path(args.file_name).parent}/{args.taxon_id}_tissue_summary.txt")
        # print(df.head())
        # Remove duplicates based on 'run_accession' while keeping the first row for each
        df_original = df.copy()
        run_accessions = filter_data(df)

        # Build a regex pattern from run_accession list
        selected_accessions = run_accessions[0:25000]
        # df_final = df_original[df_original["run_accession"].isin(run_accessions[0:250])].copy()
        df_final = df_original[df_original["run_accession"].isin(selected_accessions)].copy()
        df_final["run_accession"] = pd.Categorical(
            df_final["run_accession"], categories=selected_accessions, ordered=True
        )
        df_final = df_final.sort_values(by=["run_accession", "file_name"])  # "run_accession")

        df_final.drop_duplicates(inplace=True)
        df_final["predicted_tissue"] = df_final.apply(
            lambda row: (
                row["sample_accession"] if pd.isnull(row["tissue_prediction"]) else row["tissue_prediction"]
            ),
            axis=1,
        )
        # print(df_final.head())
        if args.csv_for_main:
            df_final.loc[:, "file_name"] = df_final["file_name"].astype(str) + ".fastq.gz"
            df_final.loc[:, "col1"] = 1
            df_final.loc[:, "col_1"] = -1
            df_final.loc[:, "col0"] = 0
            df_final.loc[:, "ENA"] = "ENA"

            with open(args.file_name, "w") as f:
                df_final.to_csv(
                    f,
                    sep="\t",
                    index=False,
                    columns=[
                        "predicted_tissue",
                        "run_accession",
                        "col1",
                        "file_name",
                        "col_1",
                        "col1",
                        "col0",
                        "ENA",
                        "platform",
                        "sample_accession",
                        "file_url",
                        "md5",
                    ],
                    header=False,
                )
        else:
            # taxon_id,gca,platform,paired,tissue,run_accession,pair1,md5_1,pair2,md5_2
            output_df = (
                df_final.groupby("run_accession")
                .apply(
                    lambda g: pd.Series(
                        {
                            "taxon_id": g["taxon_id"].iloc[0],
                            "assembly_accession": g.get("assembly_accession", pd.Series(["NA"])).iloc[
                                0
                            ],  # optional if gca exists
                            "platform": g["platform"].iloc[0],
                            "paired": bool(int(g["paired"].iloc[0])),
                            # "paired": g["paired"].iloc[0],
                            "predicted_tissue": g["predicted_tissue"].iloc[0],
                            "pair1": g[g["file_name"].str.contains("_1")]["file_url"].values[0],
                            "md5_1": g[g["file_name"].str.contains("_1")]["md5"].values[0],
                            "pair2": g[g["file_name"].str.contains("_2")]["file_url"].values[0],
                            "md5_2": g[g["file_name"].str.contains("_2")]["md5"].values[0],
                        }
                    )
                )
                .reset_index()
            )
            output_df = output_df.drop(columns=["tissue"], errors="ignore").rename(
                columns={"predicted_tissue": "tissue"}
            )

            with open(args.file_name, "w") as f:
                output_df.to_csv(
                    f,
                    sep=",",
                    index=False,
                    columns=[
                        "taxon_id",
                        "assembly_accession",
                        "platform",
                        "paired",
                        "tissue",
                        "run_accession",
                        "pair1",
                        "md5_1",
                        "pair2",
                        "md5_2",
                    ],
                    header=True,
                )

        # print("✅ CSV READY!!!!.")


if __name__ == "__main__":
    main()
