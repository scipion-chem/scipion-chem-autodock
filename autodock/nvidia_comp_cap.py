import requests
import pandas as pd
from bs4 import BeautifulSoup
import json

url = 'https://developer.nvidia.com/cuda-gpus'
data = requests.get(url).text
soup = BeautifulSoup(data, 'html.parser')

tables = soup.find_all('table')

d = {}
for table in tables:
	for row in table.tbody.find_all('tr'):
		columns = row.find_all('td')
		if(columns != []):
			d[columns[0].text.strip().lower()] = columns[1].text.strip().replace('.', '')
			
with open("NVIDIA_ComputeCapabilities.json", "w") as f:
	json.dump(d, f)

