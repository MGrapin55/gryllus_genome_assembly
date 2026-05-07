# Rules 

**TOTAL**
SUSPECT_PRODUCT_NAMES: 2524 product_names contain suspect phrases or characters  


## Fatal Rules
FATAL: SUSPECT_PRODUCT_NAMES: 6 features contain 'Low Quality Protein'
```python
# 4. FATAL: Low Quality Protein
    if re.search(r'(?i)low quality protein', name):
        name = re.sub(r'(?i)low quality protein', 'hypothetical protein', name)
        notes_to_add.append("Originally labeled as Low Quality Protein")
```

FATAL: SUSPECT_PRODUCT_NAMES: 28 features end with '.'
```python
# 3. Clean up formatting errors flagged by NCBI
    name = name.rstrip('.') # Trailing periods
```
FATAL: SUSPECT_PRODUCT_NAMES: 154 features contain unbalanced brackets or parentheses
```python
# 3. Clean up formatting errors flagged by NCBI
    name = name.replace('()', '').replace('[]', '') # Empty brackets
```
FATAL: SUSPECT_PRODUCT_NAMES: 8 features contain '. '
```python
# 3. Clean up formatting errors flagged by NCBI
    name = name.replace('. ', ' ') # Periods followed by space
```

FATAL: SUSPECT_PRODUCT_NAMES: 12 features do not contain letters in product name
```python
# Fallback
    if not name or not any(c.isalpha() for c in name):
        name = "hypothetical protein"
```


## Non Fatal Rules

SUSPECT_PRODUCT_NAMES: Putative Typo
SUSPECT_PRODUCT_NAMES: 106 features May contain plural
**Drop the plural endings changes**


SUSPECT_PRODUCT_NAMES: 2 features starts with 'gp'. May contain systematic gene product identifiers from phage  
```python
# 9. Putative Typos (Caps, 'gp', plural 'proteins')
    if name.isupper():
        name = name.lower()
    
    name = re.sub(r'(?i)^gp\s*', '', name)
```
**this error is not a  problem**


SUSPECT_PRODUCT_NAMES: 40 features are all capital letters  

> **These are all domains from what it looks like, we could names them something like (Putative <domain> domain containing protein)** 

SUSPECT_PRODUCT_NAMES: 2 features contains 'other'. Does the product name include a descriptive phrase?  
**This should just be putative sugar transporter protein**   

SUSPECT_PRODUCT_NAMES: Suspicious phrase; should this be nonfunctional?
SUSPECT_PRODUCT_NAMES: 320 features contain 'Fragment'  
**Remove (fragment) and move species to notes so its like "evidence from <species>"** 

SUSPECT_PRODUCT_NAMES: May contain database identifier more appropriate in note; remove from product name
SUSPECT_PRODUCT_NAMES: 944 features contains three or more numbers together that may be identifiers more appropriate in note
**We can leave these as they contain useful information**  

SUSPECT_PRODUCT_NAMES: 
**We can leave these for future domain queries**   

SUSPECT_PRODUCT_NAMES: 12 features contain '(TC'
```python 
# Need to remove the phrase `Belongs to the`

```


SUSPECT_PRODUCT_NAMES: Implies evolutionary relationship; change to -like protein  
SUSPECT_PRODUCT_NAMES: 122 features contain 'Homolog'  
SUSPECT_PRODUCT_NAMES: 2 features contain 'Homologue'  
```python 
# 5. Evolutionary relationship (Homolog/Homologue)
    name = re.sub(r'(?i)\bhomolog(ue)?\b', '-like protein', name)

```  


SUSPECT_PRODUCT_NAMES: Product name does not contain letters
**Need to change to hypothetical protein**

SUSPECT_PRODUCT_NAMES: Use short product name instead of descriptive phrase
SUSPECT_PRODUCT_NAMES: 196 features start with 'belongs'  
**change 'Belongs to the' to 'Putative' and then add 'protein'**


SUSPECT_PRODUCT_NAMES: 364 features Is longer than 100 characters. Remove descriptive phrases or synonyms from product names.
Keep valid long product names, eg long enzyme names  
**Remove '(<species>)' and add that the notes. evidence from <species>.**   

SUSPECT_PRODUCT_NAMES: 286 features end with 'activity'
**leave activity**

SUSPECT_PRODUCT_NAMES: use protein instead of gene as appropriate  
SUSPECT_PRODUCT_NAMES: 18 features contain 'gene'  
**I will check these manuelly** 


### New edits to make to handle fringe cases

* Any parathenesis () and that does't match 'Fragment' as contain move to notes.
* Any unopened brackets or parenthesis find the first space move remove that content otherwise make hypothetical protein if not spaces 
* Recheck to make sure there is is no (Fragment)
* double check to make sure no changes end with a .
* If it does not contain letter make it hypothetical protein 
* add a check that if the orginal matches the corrected_names the it just uses the orginal name. (I.e if its not different don't chnage it. )


### Next round of changes
* we have these 'ACWDOJ_007683' values that I think are from a wrong column or something 


### 4-2-2026 Changes 
* It appears that a lot of the discrepancies found by table2asm are being cutoff by some weird spacing or something 
* Manuelly inspecting these changes show that these changes are actually fine. I will just have to see with the NCBI team if they will accept them. 

* Need to run it twice. and need to tell it to treat [] differently

* looking at the output. Running it second time to only correct for a regrex (speices)
* and need to tell it to ignore []


## 4-3-2026
* Leaving Plurals
* Leaving GP
* Leaving Numbers 
* Leaving underscore 
* Leaving TC
* Leaving not contain letters

* Need specific regex to match species if posible. 
```python
# =========================================================
        # TIER 1: PHRASE & KEYWORD REPLACEMENTS
        # =========================================================
        # 1.1 'other' -> putative sugar transporter protein
        if re.search(r'(?i)\bother\b', name):
            name = "putative sugar transporter protein"

        # 1.2 'belongs to' -> Putative [name] protein
        if re.search(r'(?i)^belongs to (the )?', name):
            name = re.sub(r'(?i)^belongs to (the )?', 'Putative ', name)
            if not re.search(r'(?i)protein$', name):
                name = name.strip() + " protein"

        # 1.3 Evolutionary relationship (homolog -> -like)
        name = re.sub(r'(?i)\bhomolog(ue)?\b', '-like protein', name)

        # 1.4 FATAL: Low Quality Protein # Should I just add the orginal label
        if re.search(r'(?i)low quality protein', name):
            name = re.sub(r'(?i)low quality protein', 'hypothetical protein', name)
            notes_to_add.add("Originally labeled as Low Quality Protein")

        # =========================================================
        # TIER 2: PARENTHETICAL & SPECIES EXTRACTIONS
        # =========================================================
        paren_matches = re.findall(r'\(([^)]+)\)', name)
        for match in paren_matches:
            if match.lower() == 'fragment':
                pass 
            # Protect short gene symbols/numbers like the '2' in l(2)tid
            elif len(match) <= 2 and match.isalnum():
                continue 
            elif re.match(r'^[A-Z][a-z]+\s+[a-z]+(?:\s+.*)?$', match):
                notes_to_add.add(f"evidence from {match}")
            else:
                notes_to_add.add(f"Original note: {match}")
            
            name = name.replace(f"({match})", "")

        # =========================================================
        # TIER 3: STRUCTURAL ARTIFACT CLEANUP
        # =========================================================
        # 3.1 Strip 'fragment'
        name = re.sub(r'(?i)\bfragment\b', '', name)
        
        # 3.2 Clean empty or broken brackets
        name = name.replace('()', '').replace('[]', '')
        if any(char in name for char in '()[]'):
            if ' ' in name:
                words = name.split()
                clean_words = [w for w in words if not any(c in w for c in '()[]')]
                name = " ".join(clean_words)
            else:
                name = "hypothetical protein"

        # 3.3 Formatting (Periods and Whitespace)
        name = name.replace('. ', ' ') 
        name = re.sub(r'\s+', ' ', name).strip()
        name = name.rstrip('. ')

        # =========================================================
        # TIER 4: FINAL FORMATTING & FALLBACKS
        # =========================================================
        # 4.1 All caps handling (e.g., "ABCD" -> "Putative ABCD domain-containing protein")
        if name.isupper() and any(c.isalpha() for c in name) and not name.startswith("Putative"):
            name = f"Putative {name} domain-containing protein"

        # 4.2 Empty/Non-Alpha Fallback
        if not name or not any(c.isalpha() for c in name):
            name = "hypothetical protein"
```