import re
import time
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


def msearch(smiles):
    molecule_search_url = "https://rmg.mit.edu/molecule_search"
    options = Options()
    options.add_argument("--headless")
    browser = webdriver.Chrome(options=options)
    browser.get(molecule_search_url)
    elem = browser.find_element(By.NAME, "species_identifier")
    elem.send_keys(smiles + Keys.RETURN)
    time.sleep(10)
    button = browser.find_element(By.NAME, "thermo")
    button.click()
    # get contents
    contents = WebDriverWait(browser, 10).until(EC.presence_of_element_located((By.ID, 'contents')))
    # contents = browser.find_element(By.ID, "contents")
    all_text = contents.text
    res = re.search(r"CHEMKIN", all_text)
    if res is None:
        return ""
    chemcontent = all_text[res.span()[0]:]
    chemlines = []
    for l in chemcontent.split("\n"):
        chemlines.append(l)
        if l.endswith(" 4"):
            break
    chemkin_poly = "\n".join(chemlines)
    browser.close()
    return chemkin_poly
