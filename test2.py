import requests

# Step 1: Search for GSE65194 in GEO (db=geo)
search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
search_params = {
    "db": "geo",
    "term": "GSE65194[Accession]",
    "retmode": "json"
}
search_response = requests.get(search_url, params=search_params).json()
uid = search_response["esearchresult"]

print(uid)

exit()

# Step 2: Fetch summary in JSON format
summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
summary_params = {
    "db": "geo",
    "id": uid,
    "retmode": "json"
}
summary_response = requests.get(summary_url, params=summary_params).json()

# Extract platform from the JSON response
gpl = summary_response["result"][uid]["gpl"]
print(f"Platform for GSE65194: {gpl}")