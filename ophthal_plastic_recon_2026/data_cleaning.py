# GOAL: Perform various data cleaning tasks
# - Check all authors and institutions,
#   and fill in missing instutions for authors with known affiliations

import pandas as pd
import unicodedata

# Pull in data
data = pd.read_excel("ophthalmology_pubmed_dataset_full.xlsx")

# Create df for all authors + affiliations
author_cols = [
    col
    for col in data.columns
    if "author" in col
    and "name" in col
    and col not in ["all_authors_ordered", "total_authors"]
]

affiliation_cols = [col.replace("name", "affiliation") for col in author_cols]

# Create a df of two columns, the first with author names, the second with their affiliations
# drop any obs where author_name is NA
authors_df = pd.DataFrame(
    {
        "author_name": data[author_cols].values.flatten(),
        "author_affiliation": data[affiliation_cols].values.flatten(),
    }
)
authors_df = (
    authors_df.dropna(subset=["author_name"]).reset_index(drop=True).drop_duplicates()
)

# For all duplicated names, keep the first non-missing affiliation
master_author_affiliations = authors_df.copy()
master_author_affiliations["author_affiliation"] = authors_df.groupby("author_name")[
    "author_affiliation"
].transform(lambda x: x.ffill().bfill())
master_author_affiliations = master_author_affiliations.drop_duplicates(
    subset=["author_name"]
).reset_index(drop=True)

no_missing_affil = master_author_affiliations[
    master_author_affiliations["author_affiliation"].notna()
]

# Check for missing affiliations
missing_affiliations = authors_df[authors_df["author_affiliation"].isna()]

# For each author column, fill in missing affiliations from master_author_affiliations
for author_col in author_cols:
    affil_col = author_col.replace("name", "affiliation")
    data[affil_col] = data.apply(
        lambda row: (
            master_author_affiliations.loc[
                master_author_affiliations["author_name"] == row[author_col],
                "author_affiliation",
            ].values[0]
            if pd.isna(row[affil_col])
            and row[author_col] in master_author_affiliations["author_name"].values
            else row[affil_col]
        ),
        axis=1,
    )


# Save files
data.to_csv("ophthalmology_pubmed_dataset_full_authors_filled.csv", index=False)
master_author_affiliations.to_csv("master_author_affiliations.csv", index=False)


# Get a subset of all first author affiliations & last author affiliations
df_filled = pd.read_csv("ophthalmology_pubmed_dataset_firstlast_auth_complete.csv")


# Apply some cleaning to the first and last author affiliations
# by removing any leading/trailing whitespace, single lowercase letters, and periods at the end
# as well as any special characters (e.g. accents) to get a cleaner list of unique affiliations
df_filled["first_author_affiliation"] = (
    df_filled["first_author_affiliation"]
    .str.strip()
    .str.replace(r"^[a-z]\s+", "", regex=True)
    .str.replace(r"\.$", "", regex=True)
    .apply(
        lambda x: (
            unicodedata.normalize("NFKD", x).encode("ascii", "ignore").decode("ascii")
            if pd.notna(x)
            else x
        )
    )
    .str.strip()
)
df_filled["last_author_affiliation"] = (
    df_filled["last_author_affiliation"]
    .str.strip()
    .str.replace(r"^[a-z]\s+", "", regex=True)
    .str.replace(r"\.$", "", regex=True)
    .apply(
        lambda x: (
            unicodedata.normalize("NFKD", x).encode("ascii", "ignore").decode("ascii")
            if pd.notna(x)
            else x
        )
    )
    .str.strip()
)


# Combine first and last author affiliations into a single df, drop duplicates, and save to csv
affiliation_df = pd.DataFrame(
    pd.concat(
        [df_filled["first_author_affiliation"], df_filled["last_author_affiliation"]],
        axis=0,
    )
    .drop_duplicates()
    .reset_index(drop=True)
)
affiliation_df.columns = ["affiliation"]

# Minor data cleaning
# - Remove any leading/trailing whitespace, single lowercase letters, and periods at the end
clean_affiliations = (
    affiliation_df["affiliation"]
    .str.strip()
    .str.replace(r"^[a-z]\s+", "", regex=True)
    .str.replace(r"\.$", "", regex=True)
    .str.strip()
    .drop_duplicates()
    .reset_index(drop=True)
)

clean_affiliations_no_special_letters = clean_affiliations.apply(
    lambda x: (
        unicodedata.normalize("NFKD", x).encode("ascii", "ignore").decode("ascii")
        if pd.notna(x)
        else x
    )
)

clean_affiliations_no_special_letters.to_csv("first_last_author_affiliations.csv", index=False)

# Pull in cleaned up affiliations + country + continent file
affil_country_continent = pd.read_csv("organization_country_continent_complete.csv").dropna(subset=["organization"]).reset_index(drop=True)


# Join in the country and continent info to the cleaned affiliations
df_filled.drop(columns=["first_author_continent", "last_author_continent"], inplace=True)

df_filled_complete = df_filled.merge(
    affil_country_continent.rename(
        columns={"organization": "first_author_affiliation",
                 "country": "first_author_country",
                 "continent": "first_author_continent"},
        inplace=False,
    )[["first_author_affiliation", "first_author_country", "first_author_continent"]],
    on="first_author_affiliation",
    how="left",
)

df_filled_complete = df_filled_complete.merge(
    affil_country_continent.rename(
        columns={"organization": "last_author_affiliation",
                 "country": "last_author_country",
                 "continent": "last_author_continent"},
        inplace=False,
    )[["last_author_affiliation", "last_author_country", "last_author_continent"]],
    on="last_author_affiliation",
    how="left",
)

# Checking for any missing countries after the merge
missing_first_author_countries = df_filled_complete[
    df_filled_complete["first_author_country"].isna() & df_filled_complete["first_author_affiliation"].notna()
]["first_author_affiliation"].unique()

# Checking for any missing countries after the merge
missing_last_author_countries = df_filled_complete[
    df_filled_complete["last_author_country"].isna() & df_filled_complete["last_author_affiliation"].notna()
]["last_author_affiliation"].unique()

# Save the final cleaned and merged dataset
df_filled_complete.to_csv("ophthalmology_pubmed_dataset_firstlast_attrib_filled.csv", index=False)