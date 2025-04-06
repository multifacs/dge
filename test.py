import requests

# Search for GSE65194 in the GEO DataSets (gds) database
search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
params = {
    "db": "gds",
    "term": "GSE65194[Accession]",
    "retmode": "json"
}
response = requests.get(search_url, params=params).json()
uid = response["esearchresult"]["idlist"][0]  # Get the UID
print(uid)

from xml.etree import ElementTree as ET

# Fetch the record for the UID
fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
fetch_params = {
    "db": "gds",
    "id": uid,
    "retmode": "xml"
}
response = requests.get(fetch_url, params=fetch_params)
print(response.text)

# Parse the XML to find the platform (GPL)
root = ET.fromstring(response.content)
platform = None

for docsum in root.findall("DocSum"):
    for item in docsum.findall("Item"):
        if item.get("Name") == "GPL":
            platform = item.text
            break

print(f"Platform: {platform}")