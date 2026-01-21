# https://bill1224.tistory.com/381
# https://covariants.org/
# https://docs.nextstrain.org/projects/ncov/en/latest/guides/data-prep/gisaid-full.html

import datetime
import pandas as pd
from bs4 import BeautifulSoup as beautifulsoup
import time
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from html_table_parser import parser_functions as parser

options = Options()
options.add_argument("--headless")
options.add_argument("window-size=1400,1500")
options.add_argument("--disable-gpu")
options.add_argument("--no-sandbox")
options.add_argument("start-maximized")
options.add_argument("enable-automation")
options.add_argument("--disable-infobars")
options.add_argument("--disable-dev-shm-usage")

url = "https://covariants.org/"

# 드라이버 연결
browser = webdriver.Chrome(options=options)
# 웹사이트 이동
browser.get(url)
WebDriverWait(browser, 60).until(EC.presence_of_element_located((By.TAG_NAME, 'table')))
soup = beautifulsoup(browser.page_source, 'html.parser')


## corona_name.txt 생성

table = soup.find('table')
tmp = parser.make2d(table)
name = pd.DataFrame(tmp[1:], columns=tmp[0])

for c in name.columns:
    name[c].replace(c, '', regex=True, inplace=True)
    if c == "WHO Label":
        name['tmp_c'] = name[c].str[2:]
        name[c] = name['tmp_c']
        name.drop('tmp_c', axis=1, inplace=True)

name.to_csv('/files/corona_name.txt', sep='\t', index=False)


## corona_vari 생성

href = []
tmp = soup.find('tbody').find_all('tr')
for i in tmp:
    a = i.find('a').attrs['href']
    href.append(url[:-1]+a)

nextstrain = list(name['Nextstrain Clade'])
tmp = []
for i in range(len(href)):
    v =[]
    tmp_url = href[i]
    browser = webdriver.Chrome(options=options)
    browser.get(tmp_url)
    WebDriverWait(browser, 10).until(EC.presence_of_all_elements_located((By.TAG_NAME, 'ul')))
    tmp_soup = beautifulsoup(browser.page_source, 'html.parser')
    tmp_div = tmp_soup.find('div', class_='my-2 col')
    li = tmp_div.find_all('li')
    for l in li:
        variant = l.text
        v.append(variant)
    #
    tmp.append(v)

m = 0
for i in tmp:
    if len(i) > m: m = len(i)

vari = pd.DataFrame(columns=nextstrain, index=range(0, m))
for i in range(len(tmp)):
    for j in range(len(tmp[i])):
        vari.iloc[j, i] = tmp[i][j]

vari.to_csv('/files/corona_vari.txt', sep='\t', index=False)
