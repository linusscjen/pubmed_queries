# GOAL: Connect to PubMed and pull data based on requirements from Jenny

import os
import re
import time

import pandas as pd
from Bio import Entrez
from dotenv import find_dotenv, load_dotenv

# Load environment variables from .secrets.ini
load_dotenv(find_dotenv(".secrets.ini"))


# SETUP - put "email" and "pubmed_api_key" in .secrets.ini
Entrez.email = os.environ["email"]
Entrez.api_key = os.environ["pubmed_api_key"]


# -----------------------------
# CONFIG
# -----------------------------
JOURNALS = ["Ophthalmic Plast Reconstr Surg", "Orbit"]

ARTICLE_TYPES = [
    "Case Reports",
    "Clinical Study",
    "Clinical Trial",
    "Comparative Study",
    "Meta-Analysis",
    "Multicenter Study",
    "Observational Study",
    "Randomized Controlled Trial",
    "Review",
    "Systematic Review",
    "Video-Audio Media",
]


# -----------------------------
# COUNTRY LIST (simplified)
# -----------------------------
COUNTRIES = [
    "USA",
    "United States",
    "UK",
    "United Kingdom",
    "Canada",
    "Australia",
    "Germany",
    "France",
    "Italy",
    "Spain",
    "China",
    "Japan",
    "India",
    "Brazil",
    "Netherlands",
    "South Korea",
    "Sweden",
    "Switzerland",
]


# -----------------------------
# BUILD QUERY
# -----------------------------
def build_query(journal):
    # article_filter = " OR ".join([f'"{t}"[Publication Type]' for t in ARTICLE_TYPES])

    # return f"""
    # ("{journal}"[Journal])
    # AND ({article_filter})
    # AND ("2016/01/01"[PDAT] : "2025/12/31"[PDAT])
    # """
    return f"""
    ("{journal}"[Journal])
    AND ("2016/01/01"[PDAT] : "2025/12/31"[PDAT])
    """


# -----------------------------
# SEARCH
# -----------------------------
def search_pubmed(query, batch_size=100):
    ids = []
    retstart = 0

    while True:
        handle = Entrez.esearch(
            db="pubmed", term=query, retmax=batch_size, retstart=retstart
        )
        record = Entrez.read(handle)
        handle.close()

        batch_ids = record["IdList"]
        if not batch_ids:
            break

        ids.extend(batch_ids)
        retstart += batch_size

        time.sleep(0.3)

    return ids


# -----------------------------
# FETCH XML
# -----------------------------
def fetch_details(id_list):
    handle = Entrez.efetch(db="pubmed", id=",".join(id_list), retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records


# -----------------------------
# HELPERS
# -----------------------------
def normalize_name(author):
    """Return 'First Middle Last'"""
    first = author.get("ForeName", "")
    last = author.get("LastName", "")
    return f"{first} {last}".strip()


def extract_credentials(text):
    if not text:
        return None
    matches = re.findall(r"\b(MD|PhD|DO|MBBS|MS|FRCS)\b", text)
    return ", ".join(set(matches)) if matches else None


def extract_country(affiliation):
    if not affiliation:
        return None

    for country in COUNTRIES:
        if country.lower() in affiliation.lower():
            return country

    # fallback: last token
    parts = affiliation.split(",")
    if parts:
        return parts[-1].strip()

    return None


def extract_doi(article):
    try:
        id_list = article["PubmedData"]["ArticleIdList"]
        for id_obj in id_list:
            if id_obj.attributes.get("IdType") == "doi":
                return str(id_obj)
    except:
        pass
    return None


# -----------------------------
# PARSER
# -----------------------------
def parse_article(article):
    result = {}

    article_data = article.get("MedlineCitation", {}).get("Article", {})

    # ---- Title
    result["title"] = article_data.get("ArticleTitle", "")

    # ---- Journal
    result["journal"] = article_data.get("Journal", {}).get("Title", "")

    # ---- Date
    def extract_year_month(pub_date):
        if "MedlineDate" in pub_date:
            year = pub_date["MedlineDate"][:4]
            months = pub_date["MedlineDate"][5:8].strip()
        elif "Year" in pub_date and "Season" in pub_date:
            year = pub_date["Year"]
            months = pub_date["Season"][:3].strip()
        elif "Year" in pub_date and "Month" in pub_date:
            year = pub_date["Year"]
            months = pub_date["Month"]
        elif "Year" in pub_date:
            year = pub_date["Year"]
            months = ""
        else:
            raise ValueError(f"Unknown date format in PubDate, keys: {pub_date.keys()}")
        return year, months

    pub_date = (
        article_data.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
    )
    result["year"], result["month"] = extract_year_month(pub_date)
    # result["year"] = pub_date.get("Year", "")
    # result["month"] = pub_date.get("Month", "")

    # ---- DOI
    result["doi"] = extract_doi(article)

    # ---- Article type
    pub_types = article_data.get("PublicationTypeList", [])
    result["article_type"] = ", ".join(pub_types)

    # ---- Authors
    authors = article_data.get("AuthorList", [])

    names = []
    affiliations = []
    credentials = []
    countries = []

    for author in authors:
        if "LastName" in author:
            name = normalize_name(author)
            names.append(name)

            affil = ""
            if author.get("AffiliationInfo"):
                affil = author["AffiliationInfo"][0].get("Affiliation", "")

            affiliations.append(affil)

            cred = extract_credentials(affil)
            credentials.append(cred)

            country = extract_country(affil)
            countries.append(country)

    # ---- Aggregates
    result["total_authors"] = len(names)
    result["all_authors_ordered"] = " | ".join(names)

    # ---- First author
    if names:
        result["first_author_name"] = names[0]
        result["first_author_affiliation"] = affiliations[0]
        result["first_author_credentials"] = credentials[0]
        result["first_author_country"] = countries[0]

    # ---- Last author
    if names:
        result["last_author_name"] = names[-1]
        result["last_author_affiliation"] = affiliations[-1]
        result["last_author_credentials"] = credentials[-1]
        result["last_author_country"] = countries[-1]

    # ---- Middle authors (individual columns)
    if len(names) > 2:
        for idx, (name, affil, cred, country) in enumerate(
            zip(names[1:-1], affiliations[1:-1], credentials[1:-1], countries[1:-1]), 1
        ):
            result[f"middle_author_{idx}_name"] = name
            result[f"middle_author_{idx}_affiliation"] = affil
            result[f"middle_author_{idx}_credentials"] = cred
            result[f"middle_author_{idx}_country"] = country

    return result


# MAIN
def run_pipeline():
    all_results = []

    for journal in JOURNALS:
        print(f"Processing: {journal}")

        query = build_query(journal)
        ids = search_pubmed(query)

        print(f"Found {len(ids)} papers")

        for i in range(0, len(ids), 100):
            batch = ids[i : i + 100]
            records = fetch_details(batch)

            for article in records["PubmedArticle"]:
                parsed = parse_article(article)
                all_results.append(parsed)

            time.sleep(0.5)

    df = pd.DataFrame(all_results)

    # ---- Final cleanup
    # df.fillna("", inplace=True)

    df.to_csv("ophthalmology_pubmed_dataset.csv", index=False)

    print("Saved → ophthalmology_pubmed_dataset.csv")


if __name__ == "__main__":
    run_pipeline()
