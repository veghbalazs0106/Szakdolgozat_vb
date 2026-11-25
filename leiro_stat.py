# -*- coding: utf-8 -*-
"""
Cukorbetegség aluldetektálásának okai, következményei, megoldási lehetőségei – szakdolgozat (Végh Balázs, BCE, 2026)
Az elemzés teljes kódja: adatbeolvasás, adatfeldolgozás, leíró statisztika, regressziós modellek.

Minden vizsgálat (változó) ugyanazon a lépéssoron megy végig:
  1. beolvasás és előkészítés  -> hosszú formátum: Orszag | Ev | Value
  2. illesztés a népességhez   -> csak az országok maradnak bent (régiók, jövedelmi csoportok stb. kiesnek)
  3. évválasztás               -> vizsgálatonként egyetlen év: a 2022-2024 közötti legfrissebb év, amelyben az országok
                                  legalább 80%-áról van adat; ha egyikben sincs meg a 80%, a legjobb lefedettségű év
  4. nyers érték               -> leíró statisztika, dobozábra és hisztogram
  5. népességarányosítás       -> csak az esetszám / abszolút összeg jellegű mutatóknál; leíró statisztika és ábrák
  6. log transzformáció (log1p)-> az (arányosított) értékeken; leíró statisztika és ábrák
  7. javasolt modellváltozat   -> log csak akkor, ha az alapváltozat nem elfogadható (|ferdeség| vagy |csúcsosság| >= 1)
                                  és a log javít rajta

Regressziós rész:
  K modell (fő modell)  -> a kezeletlen cukorbetegek aránya (NCD-RisC, 2022): okcsoportok (hierarchikus regresszió),
                           változószelekció (LASSO-stabilitás, versengő modellek keresztvalidációval), robusztusság,
                           adatév-érzékenység, kihagyott változók ellenőrzése, modelltípus-választás, ellenőrzés az IDF
                           nem diagnosztizált arányával, megoldási forgatókönyvek
  Halálozási ellenőrzés -> B modell: halálesetek 1000 cukorbetegre ~ nem diagnosztizált arány; H1/H2: a cukorbetegek és a
                           halálesetek népességarányos számának különbsége, parciális korreláció az IDF-aránnyal

Kimenetek (a futtatás mappájában):
  leiro_statisztikak_osszesito.xlsx   minden leíró statisztika és a javasolt modellváltozat
  modell_adatok.xlsx                  a javasolt változatok országonként, a felhasznált évekkel
  regresszio_kezeletlen_arany.xlsx    a K modell minden eredménye
  regresszio_eredmenyek.xlsx          a B modell
  regresszio_halalozas_kulonbseg.xlsx a H1/H2 célváltozók és a parciális korrelációk
  regresszio_modelltipus.xlsx         a modelltípus-választás (B és K modell)
  *.png                               a leíró statisztikai és a regressziós ábrák
"""
# IMPORTOK ----------------------------------------------------------------------------------------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import csv
import re
import io
import unicodedata
import warnings
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes
import matplotlib.figure

# Az ábrák cím nélkül készülnek, a címet a dolgozatban az ábrafelirat adja.
ABRACIM = False
if not ABRACIM:
    matplotlib.axes.Axes.set_title = lambda self, *a, **k: None
    matplotlib.figure.Figure.suptitle = lambda self, *a, **k: None
    plt.title = lambda *a, **k: None
    plt.suptitle = lambda *a, **k: None

# BEÁLLÍTÁSOK -------------------------------------------------------------------------------------------------------------------------------------------------------------

ALSO_EV = 2022
FELSO_EV = 2024
LEFEDETTSEGI_KUSZOB = 0.8
NORMALITASI_KUSZOB = 1.0

# Egyetlen közös lista a Világbank-féle nem-ország sorokról.
nem_orszagok = [
    # Regionális és földrajzi aggregátumok
    'Africa Eastern and Southern',
    'Africa Western and Central',
    'Arab World',
    'Caribbean small states',
    'Central Europe and the Baltics',
    'East Asia & Pacific',
    'East Asia & Pacific (excluding high income)',
    'Euro area',
    'Europe & Central Asia',
    'Europe & Central Asia (excluding high income)',
    'European Union',
    'Latin America & Caribbean',
    'Latin America & Caribbean (excluding high income)',
    'Middle East & North Africa',
    'Middle East, North Africa, Afghanistan & Pakistan',
    'Middle East, North Africa, Afghanistan & Pakistan (excluding high income)',
    'North America',
    'Other small states',
    'Pacific island small states',
    'Small states',
    'South Asia',
    'Sub-Saharan Africa',
    'Sub-Saharan Africa (excluding high income)',
    'World',

    # Jövedelmi és gazdasági besorolások
    'Heavily indebted poor countries (HIPC)',
    'High income',
    'Least developed countries: UN classification',
    'Low & middle income',
    'Low income',
    'Lower middle income',
    'Middle income',
    'OECD members',
    'Upper middle income',
    'Not classified',

    # Demográfiai csoportok
    'Early-demographic dividend',
    'Late-demographic dividend',
    'Post-demographic dividend',
    'Pre-demographic dividend',

    # Világbank hitelezési csoportok (IDA & IBRD)
    'East Asia & Pacific (IDA & IBRD countries)',
    'Europe & Central Asia (IDA & IBRD countries)',
    'IBRD only',
    'IDA & IBRD total',
    'IDA blend',
    'IDA only',
    'IDA total',
    'Latin America & the Caribbean (IDA & IBRD countries)',
    'Middle East, North Africa, Afghanistan & Pakistan (IDA & IBRD)',
    'South Asia (IDA & IBRD)',
    'Sub-Saharan Africa (IDA & IBRD countries)',
]

# WHO- és OWID-féle régiós / összesítő sorok.
nem_orszagok_egyeb = [
    'Global', 'World', 'Africa', 'Americas', 'Europe', 'Eastern Mediterranean', 'South-East Asia', 'Western Pacific',
    'Asia', 'North America', 'South America', 'Oceania', 'European Union (27)',
    'High-income countries', 'Low-income countries', 'Lower-middle-income countries', 'Upper-middle-income countries',
]

# Ez a szótár tartalmazza azokat az elemeket, amik más megfogalmazásban de szerepelhetnek a két df-ben.
nev_szotar = {
    'United States of America': 'United States',
    'United Kingdom of Great Britain and Northern Ireland': 'United Kingdom',
    'China, Hong Kong SAR': 'Hong Kong SAR, China',
    'Hong Kong SAR': 'Hong Kong SAR, China',
    'China, Macao SAR': 'Macao SAR, China',
    'Macao SAR': 'Macao SAR, China',
    'Bolivia (Plurinational State of)': 'Bolivia',
    'Venezuela (Bolivarian Republic of)': 'Venezuela, RB',
    'Venezuela': 'Venezuela, RB',
    'Iran (Islamic Republic of)': 'Iran, Islamic Rep.',
    'Iran': 'Iran, Islamic Rep.',
    'Russia': 'Russian Federation',
    'Syria': 'Syrian Arab Republic',
    'Vietnam': 'Viet Nam',
    "Dem. People's Republic of Korea": "Korea, Dem. People's Rep.",
    "Democratic People's Republic of Korea": "Korea, Dem. People's Rep.",
    'Republic of Korea': 'Korea, Rep.',
    'Korea': 'Korea, Rep.',
    'Sint Maarten': 'Sint Maarten (Dutch part)',
    'Puerto Rico': 'Puerto Rico (US)',
    'Turkey': 'Turkiye',
    'Egypt': 'Egypt, Arab Rep.',
    'Yemen': 'Yemen, Rep.',
    'Gambia': 'Gambia, The',
    'Bahamas': 'Bahamas, The',
    'Congo': 'Congo, Rep.',
    'Democratic Republic of the Congo': 'Congo, Dem. Rep.',
    'United Republic of Tanzania': 'Tanzania',
    'Republic of Moldova': 'Moldova',
    'Kyrgyzstan': 'Kyrgyz Republic',
    "Lao People's Democratic Republic": 'Lao PDR',
    'Slovakia': 'Slovak Republic',
    'Micronesia (Federated States of)': 'Micronesia, Fed. Sts.',
    'Netherlands (Kingdom of the)': 'Netherlands',
    # NCD-RisC-féle írásmódok
    'DR Congo': 'Congo, Dem. Rep.',
    'Federated States of Micronesia': 'Micronesia, Fed. Sts.',
    'Guinea Bissau': 'Guinea-Bissau',
    'North Korea': "Korea, Dem. People's Rep.",
    'South Korea': 'Korea, Rep.',
    'State of Palestine': 'West Bank and Gaza',
    'occupied Palestinian territory': 'West Bank and Gaza',
    'occupied Palestinian territory, including east Jerusalem': 'West Bank and Gaza',
    'Palestine, State of': 'West Bank and Gaza',
    'Palestine': 'West Bank and Gaza',
    'Czech Republic': 'Czechia',
    'Somalia': 'Somalia, Fed. Rep.',
    'Brunei': 'Brunei Darussalam',
    'Cape Verde': 'Cabo Verde',
    'Swaziland': 'Eswatini',
    'Republic of North Macedonia': 'North Macedonia',
    'Nauru': 'Naoero',
    'United States Virgin Islands': 'Virgin Islands (U.S.)',
    # WHO-féle írásmódok (20-70 évesek halálozási aránya)
    'Kosovo (under UNSC res. 1244)': 'Kosovo',
    'Micronesia (Fed. States of)': 'Micronesia, Fed. Sts.',
    'Republic of Congo': 'Congo, Rep.',
    # UNDP-féle írásmódok (HDI)
    'Hong Kong, China (SAR)': 'Hong Kong SAR, China',
    'Korea (Republic of)': 'Korea, Rep.',
    "Korea (Democratic People's Rep. of)": "Korea, Dem. People's Rep.",
    'Moldova (Republic of)': 'Moldova',
    'Eswatini (Kingdom of)': 'Eswatini',
    'Tanzania (United Republic of)': 'Tanzania',
    'Congo (Democratic Republic of the)': 'Congo, Dem. Rep.',
}

# FÜGGVÉNYEK --------------------------------------------------------------------------------------------------------------------------------------------------------------

# Világbank adatbeolvasáshoz
def load_worldbank_csv(file_path):
    with open(file_path, 'r', encoding='utf-8-sig') as f:
        lines = [line.strip('\r\n') for line in f if line.strip()]

    header = lines[0].split(',')

    data = []
    for line in lines[1:]:
        cleaned_line = line.replace('""', '"')

        if cleaned_line.startswith('"') and cleaned_line.endswith('"'):
            cleaned_line = cleaned_line[1:-1]

        row = list(csv.reader([cleaned_line]))[0]

        if len(row) > len(header):
            extra = len(row) - len(header)
            merged_country_name = ','.join(row[: extra + 1])
            row = [merged_country_name] + row[extra + 1 :]

        data.append(row)

    df = pd.DataFrame(data, columns=header)

    for col in df.columns[4:]:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    return df


# Világbank-CSV beolvasása az összes előforduló formában (a fejléc helyét és a formátumot maga ismeri fel):
def load_worldbank_rugalmas(file_path):
    with open(file_path, 'r', encoding='utf-8-sig') as f:
        sorok = f.readlines()
    tisztit = lambda sor: sor.strip().lstrip('\ufeff').lstrip('"').lstrip()
    fejlec = next((i for i, sor in enumerate(sorok) if tisztit(sor).startswith('Country Name')), None)
    if fejlec is None:
        raise ValueError(f'{file_path}: nem található "Country Name" fejléc. Az első 3 sor:\n' + ''.join(sorok[:3]))
    csomagolt_minta = re.compile(r'^\ufeff?"[^"]*,""')
    def kicsomagol(sor):
        sor = sor.rstrip('\r\n').lstrip('\ufeff')
        if csomagolt_minta.match(sor):
            sor = sor.replace('""', '"')
            if sor.startswith('"'):
                sor = sor[1:]
            if sor.endswith('"') and sor.count('"') % 2 == 1:
                sor = sor[:-1]
        return sor
    csomagolt = any(csomagolt_minta.match(sor) for sor in sorok[fejlec + 1:fejlec + 6])
    if csomagolt:
        tiszta = [kicsomagol(sor) for sor in sorok[fejlec:] if sor.strip()]
        fej = next(csv.reader([tiszta[0]]))
        adat = []
        for sor in tiszta[1:]:
            mezok = next(csv.reader([sor]))
            if len(mezok) > len(fej):   # vesszős országnév (pl. "Bahamas, The") idézőjel nélkül -> összevonás
                tobb = len(mezok) - len(fej)
                mezok = [','.join(mezok[:tobb + 1])] + mezok[tobb + 1:]
            adat.append(mezok + [''] * (len(fej) - len(mezok)))
        df = pd.DataFrame(adat, columns=fej)
        sep = ','
    else:
        fejlec_sor = sorok[fejlec].strip()
        sep = ';' if fejlec_sor.count(';') > fejlec_sor.count(',') else ','
        df = pd.read_csv(file_path, skiprows=fejlec, sep=sep, encoding='utf-8-sig', dtype=str)
    df.columns = [str(c).strip().strip('"').strip() for c in df.columns]
    df.columns = [re.match(r'^(\d{4})', c).group(1) if re.match(r'^\d{4}\b', c) else c for c in df.columns]
    df = df.loc[:, [c for c in df.columns if c and not c.startswith('Unnamed')]]
    if 'Country Name' not in df.columns:
        raise ValueError(f'{file_path}: a beolvasás után sincs "Country Name" oszlop. Oszlopok: {list(df.columns)[:6]}')
    df = df[df['Country Name'].notna() & (df['Country Name'].astype(str).str.strip() != '')]
    for col in [c for c in df.columns if c.isdigit()]:
        ertek = df[col].astype(str).str.strip().str.strip('"')
        if sep == ';':
            ertek = ertek.str.replace(',', '.', regex=False)   
        df[col] = pd.to_numeric(ertek.replace({'..': None, '': None}), errors='coerce')
    return df


# NCD-RisC (ncdrisc.org, "Country-specific data for all countries", diabetes) beolvasása.
def load_ncdrisc(file_path):
    with open(file_path, 'r', encoding='utf-8-sig') as f:
        sorok = [sor.rstrip('\r\n') for sor in f if sor.strip()]
    def kicsomagol(sor):
        if re.match(r'^"[^"]*,""', sor):
            sor = sor[1:]
            if sor.endswith('"'):
                sor = sor[:-1]
            sor = sor.replace('""', '"')
        return sor
    mezok = [next(csv.reader([kicsomagol(sor)])) for sor in sorok]
    df = pd.DataFrame(mezok[1:], columns=[c.strip() for c in mezok[0]])
    orszag_oszlop = next(c for c in df.columns if c.startswith('Country'))
    kezelt = next(c for c in df.columns if 'treated' in c and 'uncertainty' not in c)
    kezelt_also = next(c for c in df.columns if 'treated' in c and 'lower' in c)
    kezelt_felso = next(c for c in df.columns if 'treated' in c and 'upper' in c)
    gyak = next(c for c in df.columns if c.startswith('Prevalence') and 'uncertainty' not in c)
    if 'Age' in df.columns:
        df = df[df['Age'].str.contains('standard', case=False)]
    for c in [kezelt, kezelt_also, kezelt_felso, gyak]:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df['Ev'] = pd.to_numeric(df['Year'], errors='coerce').astype(int)
    if df[kezelt].max() <= 1:
        for c in [kezelt, kezelt_also, kezelt_felso]:
            df[c] = df[c] * 100
    if 'Sex' in df.columns and not df['Sex'].str.contains('Both', case=False).any():
        for uj, c in [('Value', kezelt), ('CI_also', kezelt_also), ('CI_felso', kezelt_felso)]:
            df[uj] = df[c] * df[gyak]
        ki = df.groupby([orszag_oszlop, 'Ev'])[['Value', 'CI_also', 'CI_felso', gyak]].sum().reset_index()
        for uj in ['Value', 'CI_also', 'CI_felso']:
            ki[uj] = ki[uj] / ki[gyak]
        ki = ki.drop(columns=gyak)
    else:
        df = df[df['Sex'].str.contains('Both', case=False)] if 'Sex' in df.columns else df
        ki = df.rename(columns={kezelt: 'Value', kezelt_also: 'CI_also', kezelt_felso: 'CI_felso'})[[orszag_oszlop, 'Ev', 'Value', 'CI_also', 'CI_felso']]
    ki = ki.rename(columns={orszag_oszlop: 'Orszag'}).dropna(subset=['Value'])
    print(f'NCD-RisC: {ki["Orszag"].nunique()} ország, {ki["Ev"].min()}–{ki["Ev"].max()}; mutató: {kezelt} (a nemek gyakorisággal súlyozva)')
    ki.attrs['ev_forras'] = 'fájlból (évoszlop)'
    return ki


# A soronként idézőjelbe csomagolt, pontosvesszős fejlécű WHO-exportokhoz (pl. orvosok száma 10 000 főre)
def load_hibas_who_csv(file_path):
    cleaned_lines = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            if i == 0 and line.endswith(';'):
                line = line[:-1]
            if line.startswith('"') and line.endswith('"'):
                line = line[1:-1].replace('""', '"')
            cleaned_lines.append(line)
    return pd.read_csv(io.StringIO('\n'.join(cleaned_lines)))


# Névtisztítás
def tisztitott_nev(text):
    if not isinstance(text, str):
        return text
    text = unicodedata.normalize('NFKD', text).encode('ASCII', 'ignore').decode('utf-8')
    text = text.replace('St.', 'Saint').replace('&', 'and').strip()
    return text


# Illesztési kulcs: tisztított név -> szótár szerinti egységes név -> kisbetűs
def illesztesi_kulcs(nevek):
    return nevek.astype(str).map(tisztitott_nev).replace(nev_szotar).map(str.lower)


def who_hosszu(df, ev=None, dim1=None, orszag_oszlop='Location', ev_oszlop='Period', ertek_oszlop='Value'):
    if dim1 is not None and 'Dim1' in df.columns:
        df = df[df['Dim1'] == dim1]
    if ev_oszlop in df.columns:
        evek = df[ev_oszlop]
        ev_forras = 'fájlból (évoszlop)'
    elif ev is not None:
        evek = ev
        ev_forras = 'feltételezett (nincs évoszlop, a fájlnév alapján)'
        print(f'MEGJEGYZÉS: a fájlban nincs "{ev_oszlop}" oszlop, az év a fájl forrása (Diabetes Atlas) alapján: {ev}.')
    else:
        raise ValueError('Nincs évoszlop a táblában, add meg az "ev" paramétert.')
    ertekek = szamma(df[ertek_oszlop])
    # Az újabb WHO-exportokban a 'Value' szöveges kijelzés, a szám a 'FactValueNumeric' oszlopban van.
    if ertekek.isna().all() and 'FactValueNumeric' in df.columns:
        ertekek = szamma(df['FactValueNumeric'])
    ki = pd.DataFrame({'Orszag': df[orszag_oszlop], 'Ev': evek, 'Value': ertekek}).dropna(subset=['Orszag', 'Value'])
    ki.attrs['ev_forras'] = ev_forras
    if ki.empty:
        raise ValueError(
            f'Nem maradt egyetlen számértékű sor sem. Oszlopok: {list(df.columns)}\n'
            f'A(z) "{ertek_oszlop}" oszlop első értékei: {df[ertek_oszlop].head(5).tolist()}'
        )
    return ki


# Szöveges értékek számmá alakítása
def szamma(sorozat):
    if pd.api.types.is_numeric_dtype(sorozat):
        return sorozat
    s = sorozat.astype(str).str.replace(r'\s', '', regex=True)
    return pd.to_numeric(s.str.extract(r'^(-?[0-9]+(?:[.,][0-9]+)?)')[0].str.replace(',', '.'), errors='coerce')

# Világbank-féle széles (évek oszlopokban) táblából hosszú formátum, a nem-ország sorok nélkül.
def wb_hosszu(df):
    if 'Indicator Name' in df.columns and df['Indicator Name'].notna().any():
        print('Mutató:', df['Indicator Name'].dropna().iloc[0])
    evek = [c for c in df.columns if str(c).isdigit()]
    ki = df[~df['Country Name'].isin(nem_orszagok)].melt(
        id_vars='Country Name', value_vars=evek, var_name='Ev', value_name='Value'
    )
    ki = ki.rename(columns={'Country Name': 'Orszag'})
    ki['Ev'] = ki['Ev'].astype(int)
    ki = ki.dropna(subset=['Value'])
    ki.attrs['ev_forras'] = 'fájlból (évoszlop)'
    return ki


# A felhasznált év felirata ábracímekhez (vizsgálatonként egyetlen év van, pl. '2023').
def ev_felirat(df):
    evek = sorted(df['Ev'].dropna().astype(int).unique())
    if len(evek) == 1:
        return str(evek[0])
    return f'{evek[0]}-{evek[-1]} (országonként a legfrissebb)'


# Ábrák készítéséhez
def abrak_keszitese_es_mentese(
    data_series,
    boxplot_title,
    boxplot_ylabel,
    boxplot_filename,
    hist_title,
    hist_xlabel,
    hist_ylabel='Gyakoriság (Esetszám)',
    hist_filename=None,
    figsize_box=(6, 7),
    figsize_hist=(8, 5),
    bins=25,
    color='#4C72B0'
):
    df = pd.DataFrame({'Értékek': data_series})

    sns.set_theme(style="whitegrid")

    plt.figure(figsize=figsize_box)

    # Dobozábra
    boxplot = sns.boxplot(
        data=df,
        y='Értékek',
        color=color,
        width=0.35,
        linewidth=2,
        fliersize=5,
        boxprops=dict(alpha=0.85)
    )

    boxplot.set_title(boxplot_title, fontsize=14, pad=15, fontweight='bold')
    boxplot.set_ylabel(boxplot_ylabel, fontsize=12, labelpad=10)

    plt.tight_layout()
    plt.savefig(boxplot_filename, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    # Hisztogram
    plt.figure(figsize=figsize_hist)

    hist = sns.histplot(
        data=df,
        x='Értékek',
        kde=True,
        color=color,
        bins=bins,
        alpha=0.6,
        edgecolor='white',
        linewidth=1.2
    )

    hist.set_title(hist_title, fontsize=14, pad=15, fontweight='bold')
    hist.set_xlabel(hist_xlabel, fontsize=12, labelpad=10)
    hist.set_ylabel(hist_ylabel, fontsize=12, labelpad=10)

    plt.tight_layout()
    if hist_filename:
        plt.savefig(hist_filename, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


# LEÍRÓ STATISZTIKA SZÁMÍTÁS
def leiro_statisztika_szamitasa(data_series):

    clean_series = data_series.dropna()

    q1 = clean_series.quantile(0.25)
    q3 = clean_series.quantile(0.75)
    iqt = q3 - q1

    gyakorisagok = clean_series.value_counts()
    max_gyakorisag = int(gyakorisagok.max()) if not gyakorisagok.empty else 0

    if max_gyakorisag <= 1:
        moduszok = []
        modusz_jelleg = 'nincs módusz (minden érték egyszer fordul elő)'
    else:
        moduszok = sorted(gyakorisagok[gyakorisagok == max_gyakorisag].index.tolist())
        modusz_jelleg = 'egymóduszú' if len(moduszok) == 1 else f'többmóduszú ({len(moduszok)} módusz)'

    statisztikak = {
        'Átlag': np.mean(clean_series),
        'Medián': np.median(clean_series),
        'Módusz': moduszok[0] if len(moduszok) == 1 else np.nan,
        'Módusz jellege': modusz_jelleg,
        'Móduszok száma': len(moduszok),
        'Összes módusz': '; '.join(f'{m:.6g}' for m in moduszok) if moduszok else '-',
        'Módusz gyakorisága': max_gyakorisag if moduszok else np.nan,
        'Szórás': np.std(clean_series, ddof=0),
        'Elemszám': len(clean_series),
        'Minimum': clean_series.min(),
        'Maximum': clean_series.max(),
        'Terjedelem': clean_series.max() - clean_series.min(),
        'Ferdeség (Skewness)': clean_series.skew(),
        'Csúcsosság (Kurtosis)': clean_series.kurtosis(),
        'Q1 (25%)': q1,
        'Q3 (75%)': q3,
        'IQT (Interkvartilis terjedelem)': iqt
    }

    return pd.Series(statisztikak)


# A népességi adat nélküli nevek kiírása (a névillesztés ellenőrzéséhez).
def illesztesi_jelentes(nev, hosszu):
    adat_nevek = pd.Series(hosszu['Orszag'].dropna().unique())
    nep_kulcsok = set(nepesseg_hosszu['Match_Key'])
    nincs = sorted(n for n, k in zip(adat_nevek, illesztesi_kulcs(adat_nevek)) if k not in nep_kulcsok)
    print(f'Népességi adat nélküli nevek a(z) {nev} adatban ({len(nincs)} db) - a nyers statisztikában bent maradnak, '
          f'a népességarányosításnál kimaradnak: {nincs}')


# A VIZSGÁLATOK KÖZÖS SZÁMÍTÁSI LÉPÉSSORA (ábrát nem készít, azok vizsgálatonként külön szerepelnek lent)

eredmenyek = {}   # vizsgálatonként a feldolgozott tábla
osszesito = {}    # (vizsgálat, lépés) -> leíró statisztika
bemenetek = {}    # a vizsgálatok bemenete (az adatév-érzékenységi ellenőrzéshez)

def vizsgalat(nev, hosszu, aranyositas=True, egyseg_szorzo=1, per_fo=1, also_ev=None, felso_ev=None, mentes=True):
    also_ev = ALSO_EV if also_ev is None else also_ev
    felso_ev = FELSO_EV if felso_ev is None else felso_ev
    print('\n' + '=' * 100 + f'\n{nev}  (évablak: {also_ev}–{felso_ev}{"" if mentes else ", ellenőrző futás"})\n' + '=' * 100)

    if mentes:
        bemenetek[nev] = dict(hosszu=hosszu.copy(), aranyositas=aranyositas, egyseg_szorzo=egyseg_szorzo, per_fo=per_fo,
                              ev_forras=hosszu.attrs.get('ev_forras', 'fájlból (évoszlop)'))
        illesztesi_jelentes(nev, hosszu)

    ev_forras = hosszu.attrs.get('ev_forras', 'fájlból (évoszlop)')
    hosszu = hosszu.copy()
    hosszu['Match_Key'] = illesztesi_kulcs(hosszu['Orszag'])
    hosszu['Ev'] = hosszu['Ev'].astype(int)

    duplikalt = hosszu[hosszu.duplicated(['Match_Key', 'Ev'], keep=False)]
    if not duplikalt.empty:
        print('FIGYELEM: egy országhoz egy évben több sor is tartozik, az első marad meg:')
        print(duplikalt[['Orszag', 'Ev', 'Value']].to_string(index=False))
        hosszu = hosszu.drop_duplicates(['Match_Key', 'Ev'], keep='first')

    aggregatum_kulcsok = set(illesztesi_kulcs(pd.Series(nem_orszagok + nem_orszagok_egyeb)))
    aggregatum = hosszu['Match_Key'].isin(aggregatum_kulcsok) | hosszu['Orszag'].astype(str).str.contains(r'\((?:WHO|WB)\)', regex=True)
    if aggregatum.any():
        print(f'Kiszűrt régiós / összesítő sorok: {sorted(hosszu.loc[aggregatum, "Orszag"].unique())}')
    hosszu = hosszu[~aggregatum]

    nep = nepesseg_hosszu[['Match_Key', 'Ev', 'Nepesseg', 'Orszag']].rename(columns={'Orszag': 'Country Name'})
    jeloltek = hosszu[hosszu['Ev'].between(also_ev, felso_ev)].merge(nep, on=['Match_Key', 'Ev'], how='left')
    jeloltek['Country Name'] = jeloltek['Country Name'].fillna(jeloltek['Orszag'])
    if jeloltek.empty:
        raise ValueError(f'{nev}: nincs {also_ev}-{felso_ev} közötti adat, adathiány miatt a vizsgálat nem végezhető el.')

    orszagok_szama = nep.groupby('Ev')['Match_Key'].nunique()
    lefedettseg = pd.DataFrame({'db': jeloltek[jeloltek['Nepesseg'].notna()].groupby('Ev')['Match_Key'].nunique()})
    lefedettseg = lefedettseg.reindex(sorted(jeloltek['Ev'].unique()), fill_value=0)
    lefedettseg['osszes'] = orszagok_szama.reindex(lefedettseg.index)
    lefedettseg['arany'] = lefedettseg['db'] / lefedettseg['osszes']

    megfelelo_evek = lefedettseg[lefedettseg['arany'] >= LEFEDETTSEGI_KUSZOB]
    if not megfelelo_evek.empty:
        valasztott_ev = int(megfelelo_evek.index.max())
        lefedettsegi_megjegyzes = f'legalább {LEFEDETTSEGI_KUSZOB:.0%}'
    else:
        valasztott_ev = int(lefedettseg.sort_index().iloc[::-1]['arany'].idxmax())
        lefedettsegi_megjegyzes = f'FIGYELEM: egyik évben sem éri el a {LEFEDETTSEGI_KUSZOB:.0%}-ot'
        print(f'FIGYELEM: {nev} - egyik {also_ev}-{felso_ev} közötti évben sem éri el a lefedettség a {LEFEDETTSEGI_KUSZOB:.0%}-ot, '
              f'a legjobb lefedettségű év ({valasztott_ev}) kerül felhasználásra.')

    lefedettseg_szoveg = '; '.join(
        f"{ev}: {int(r['db'])}/{int(r['osszes'])} ország ({r['arany']:.0%}){' <- választott' if ev == valasztott_ev else ''}"
        for ev, r in lefedettseg.sort_index().iterrows()
    )
    valasztott_lefedettseg = lefedettseg.loc[valasztott_ev, 'arany']

    df = jeloltek[jeloltek['Ev'] == valasztott_ev].copy()

    # Azok az országok, amelyekről van adat a fájlban, de a választott évből hiányoznak.
    kulcs_nev = dict(zip(hosszu['Match_Key'], hosszu['Orszag']))
    kulcs_nev.update(dict(zip(nep['Match_Key'], nep['Country Name'])))
    hianyzik = sorted(kulcs_nev[k] for k in set(hosszu['Match_Key']) - set(df['Match_Key']))
    if hianyzik:
        print(f'A választott évből ({valasztott_ev}) hiányzó, ezért kimaradó országok ({len(hianyzik)} db): {hianyzik}')
    print(f'Felhasznált év: {valasztott_ev}  |  lefedettség évenként: {lefedettseg_szoveg}  |  az év forrása: {ev_forras}')

    if aranyositas:
        df['Aranyos'] = df['Value'] * egyseg_szorzo / df['Nepesseg'] * per_fo
        nep_nelkul = sorted(df.loc[df['Nepesseg'].isna(), 'Country Name'])
        if nep_nelkul:
            print(f'Népességi adat hiányában a népességarányos és a log statisztikából kimaradó országok '
                  f'({len(nep_nelkul)} db, a nyersben bent vannak): {nep_nelkul}')
    else:
        df['Aranyos'] = df['Value']

    df['Log_Value'] = np.log1p(df['Aranyos'])
    df = df[['Country Name', 'Ev', 'Value', 'Nepesseg', 'Aranyos', 'Log_Value']].sort_values('Country Name')

    if not mentes:
        return df

    lepesek = {'nyers': 'Value'}
    if aranyositas:
        lepesek['nepessegaranyos'] = 'Aranyos'
    lepesek['log'] = 'Log_Value'

    for kulcs in [k for k in osszesito if k[0] == nev]:
        del osszesito[kulcs]

    ev_info = pd.Series({
        'Felhasznált év': valasztott_ev,
        'Lefedettség a választott évben': round(float(valasztott_lefedettseg), 3),
        'Lefedettség évenként': lefedettseg_szoveg,
        'Lefedettségi megjegyzés': lefedettsegi_megjegyzes,
        'Az év forrása': ev_forras,
    })
    for lepes, oszlop in lepesek.items():
        osszesito[(nev, lepes)] = pd.concat([ev_info, leiro_statisztika_szamitasa(df[oszlop])])
    print(pd.DataFrame({lepes: osszesito[(nev, lepes)] for lepes in lepesek}))
    for lepes in lepesek:
        if osszesito[(nev, lepes)]['Móduszok száma'] > 1:
            print(f'FIGYELEM: a(z) {lepes} adatsor többmóduszú: {osszesito[(nev, lepes)]["Összes módusz"]}')

    eredmenyek[nev] = df
    return df


# KIEGÉSZÍTŐ KIÍRÁSOK A DOLGOZAT SZÖVEGÉHEZ -----------------------------------------------------------------------------------------------------------------------

def _szam(v, tized=2):
    if pd.isna(v):
        return '-'
    return f'{v:,.{tized}f}'.replace(',', ' ').replace('.', ',')


def reszletes_jelentes(df, cim, tized=2):
    """A dolgozatban használt kiegészítő értékek kiírása a nyers és (ha van) a népességarányos változatra."""
    print('\n' + '-' * 100 + f'\nKIEGÉSZÍTŐ ÉRTÉKEK – {cim}  (év: {ev_felirat(df)})\n' + '-' * 100)
    oszlopok = [('nyers', 'Value')]
    if not df['Aranyos'].equals(df['Value']):
        oszlopok.append(('népességarányos', 'Aranyos'))
    for nev, oszlop in oszlopok:
        s = df.set_index('Country Name')[oszlop].dropna()
        if s.empty:
            continue
        q1, med, q3 = s.quantile(0.25), s.median(), s.quantile(0.75)
        terj = s.max() - s.min()
        print(f'[{nev}]  elemszám: {len(s)}')
        print('   3 legkisebb: ' + '; '.join(f'{o} ({_szam(v, tized)})' for o, v in s.nsmallest(3).items()))
        print('   3 legnagyobb: ' + '; '.join(f'{o} ({_szam(v, tized)})' for o, v in s.nlargest(3).items()))
        print(f'   terjedelem: {_szam(terj, tized)}', end='')
        if s.min() > 0:
            print(f'  |  maximum / minimum: {_szam(s.max() / s.min(), 1)}-szoros', end='')
        print(f'  |  IQT a terjedelem %-ában: {_szam((q3 - q1) / terj * 100, 1)}%'
              f'  |  relatív szórás (szórás / átlag): {_szam(np.std(s) / s.mean() * 100, 1)}%')
        # módusz országai
        gyak = s.round(10).value_counts()
        if not gyak.empty and gyak.max() > 1:
            for m in sorted(gyak[gyak == gyak.max()].index):
                print(f'   módusz {_szam(m, tized)} ({gyak.max()} ország): {", ".join(sorted(s[s.round(10) == m].index))}')
        # Magyarország helye
        if 'Hungary' in s.index:
            hu = s['Hungary']
            if hu < q1: hely = 'az alsó kvartilis alatt'
            elif hu < med: hely = 'az alsó kvartilis és a medián között'
            elif hu <= q3: hely = 'a medián és a felső kvartilis között'
            else: hely = 'a felső kvartilis fölött (felső negyed)'
            print(f'   Magyarország: {_szam(hu, tized)}  |  {hely}  |  az átlag {"fölött" if hu > s.mean() else "alatt"}'
                  f'  |  az országok {_szam((s < hu).mean() * 100, 1)}%-ánál nagyobb')
        else:
            print('   Magyarország: nincs adat')


def sav_arany(df, oszlop, also=None, felso=None, cimke=''):
    """Hány ország esik egy adott értéksávba (a dolgozat szövegében szereplő sávos állításokhoz)."""
    s = df[oszlop].dropna()
    feltetel = pd.Series(True, index=s.index)
    if also is not None: feltetel &= s >= also
    if felso is not None: feltetel &= s <= felso
    sav = f'{also if also is not None else "-"} … {felso if felso is not None else "-"}'
    print(f'   {cimke} [{oszlop}] {sav}: {feltetel.sum()} ország ({_szam(feltetel.mean() * 100, 1)}%)')
    return df.loc[feltetel[feltetel].index, 'Country Name'].tolist()


# NÉPESSÉG (minden vizsgálat ehhez illesztődik) ---------------------------------------------------------------------------------------------------------------------------

teljes_pop = load_worldbank_csv('teljes_populacio.csv')
nepesseg_hosszu = wb_hosszu(teljes_pop).rename(columns={'Value': 'Nepesseg'})
nepesseg_hosszu['Match_Key'] = illesztesi_kulcs(nepesseg_hosszu['Orszag'])


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ELSŐ VIZSGÁLAT - 20 és 70 éves kor közötti cukorbetegségre visszavezethető halálozási arány 2024

felnottkori_cukorbeteg_halalozas_arany = pd.read_csv("20_es_70_ev_kozottiek_cukorbetegsegre_visszavezetheto_halalozasi_aranyai_regionkent_orszagonkent_2024.csv", sep=',')

df_felnottkori_arany = vizsgalat('felnottkori_cukorbeteg_halalozasi_arany',
                                 who_hosszu(felnottkori_cukorbeteg_halalozas_arany, ev=2024),
                                 aranyositas=False)
df_felnottkori_arany.to_excel('felnottkori_cukorbeteg_halalozasi_arany_2024.xlsx', index=False)
reszletes_jelentes(df_felnottkori_arany, '1. vizsgálat – 20–79 évesek cukorbetegség miatti halálozási aránya', 1)

# Nyers értékek - dobozábra és hisztogram
abrak_keszitese_es_mentese(df_felnottkori_arany['Value'],
                           f'Felnőttkori cukorbetegségre visszavezethető halálesetek aránya - {ev_felirat(df_felnottkori_arany)}',
                           'Cukorbetegség miatti halálesetek aránya',
                           'felnottkori_cukorbeteg_halalozasi_arany_boxplot.png',
                           hist_title=f'Felnőttkori cukorbetegségre visszavezethető halálesetek eloszlása - {ev_felirat(df_felnottkori_arany)}',
                           hist_xlabel='Cukorbetegség miatti halálesetek aránya (%)',
                           hist_filename='felnottkori_cukorbeteg_halalozasi_arany_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_felnottkori_arany['Log_Value'],
                           f'Felnőttkori cukorbetegségre visszavezethető halálesetek aránya log transzformáció után - {ev_felirat(df_felnottkori_arany)}',
                           'Halálesetek aránya (log)',
                           'felnottkori_cukorbeteg_halalozasi_arany_log_boxplot.png',
                           hist_title=f'Felnőttkori cukorbetegségre visszavezethető halálesetek eloszlása log transzformáció után - {ev_felirat(df_felnottkori_arany)}',
                           hist_xlabel='Cukorbetegség miatti halálesetek aránya (log)',
                           hist_filename='felnottkori_cukorbeteg_halalozasi_arany_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# MÁSODIK VIZSGÁLAT - cukorbetegek száma országos bontásban 2024

cukorbetegek_szama_orsz_2024 = pd.read_csv("cukorbetegek_szama_ezerben_orszagonkent_2024.csv", sep=',')

df_cukorbetegszam = vizsgalat('cukorbetegek_szama',
                              who_hosszu(cukorbetegek_szama_orsz_2024, ev=2024),
                              egyseg_szorzo=1000, per_fo=100)
df_cukorbetegszam.to_excel('cukorbetegek_szama_orszagonkent_2024.xlsx', index=False)
reszletes_jelentes(df_cukorbetegszam, '2. vizsgálat – cukorbetegek száma (ezer fő; arányos: a népesség %-a)')
print('   Kína: cukorbetegek a népesség %-ában:', _szam(df_cukorbetegszam.loc[df_cukorbetegszam['Country Name'] == 'China', 'Aranyos'].iloc[0], 2))

# Nyers értékek
abrak_keszitese_es_mentese(df_cukorbetegszam['Value'],
                           f'Cukorbetegek száma ezer főben - {ev_felirat(df_cukorbetegszam)}',
                           'Cukorbetegek száma ezer főben',
                           'cukorbetegek_szama_ezer_foben_2024_boxplot.png',
                           hist_title=f'Cukorbetegek számának eloszlása ezer főben mérve - {ev_felirat(df_cukorbetegszam)}',
                           hist_xlabel='Cukorbetegek száma ezer főben',
                           hist_filename='cukorbetegek_szama_ezer_foben_2024_hist.png')

# Népességarányosítás után
abrak_keszitese_es_mentese(df_cukorbetegszam['Aranyos'],
                           f'Népességszámarányos cukorbetegeloszlás országos szinten - {ev_felirat(df_cukorbetegszam)}',
                           'Cukorbetegek aránya a népességhez képest (%)',
                           'nepessegszamaranyos_cukorbetegszam.png',
                           hist_title=f'Cukorbetegek népességarányos eloszlása - {ev_felirat(df_cukorbetegszam)}',
                           hist_xlabel='Cukorbetegek népességaránya (%)',
                           hist_filename='cukorbetegek_nepessegaranyos_eloszlasa_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_cukorbetegszam['Log_Value'],
                           f'Népességarányos cukorbetegszám log transzformáció után - {ev_felirat(df_cukorbetegszam)}',
                           'Cukorbetegek népességaránya (log)',
                           'cukorbetegek_nepessegaranyos_log_boxplot.png',
                           hist_title=f'Cukorbetegek népességarányos eloszlása log transzformáció után - {ev_felirat(df_cukorbetegszam)}',
                           hist_xlabel='Cukorbetegek népességaránya (log)',
                           hist_filename='cukorbetegek_nepessegaranyos_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# HARMADIK VIZSGÁLAT - cukorbetegséggel kapcsolatos egészségügyi kiadás egy főre, USD 2024

cuk_egeszs_kiad_egyemberre_2024 = pd.read_csv('cukorbetegseggel_kapcsolatos_egeszsegugyi_kiadasok_egy_emberre_USD_regionkent_orszagonkent_2024.csv', sep=',')

df_cuk_egeszs_kiad = vizsgalat('cuk_egeszs_kiad_egy_fore',
                               who_hosszu(cuk_egeszs_kiad_egyemberre_2024, ev=2024),
                               aranyositas=False)
df_cuk_egeszs_kiad.to_excel('cukorbeteg_egeszsegugyi_kiadasok_2024.xlsx', index=False)
reszletes_jelentes(df_cuk_egeszs_kiad, '3. vizsgálat – kiadás egy cukorbetegre (USD)', 1)
sav_arany(df_cuk_egeszs_kiad, 'Value', felso=2000, cimke='2000 USD alatti')

# Nyers értékek
abrak_keszitese_es_mentese(df_cuk_egeszs_kiad['Value'],
                           f'Cukorbetegséggel kapcsolatos egészségügyi kiadások dobozábrája - {ev_felirat(df_cuk_egeszs_kiad)}',
                           'Cukorbetegséggel kapcsolatos egészségügyi kiadások (USD)',
                           'cukorbeteg_egeszsegugyi_kiadasok_dollar_2024_boxplot.png',
                           hist_title=f'Cukorbetegséggel kapcsolatos egészségügyi kiadások USD/fő - {ev_felirat(df_cuk_egeszs_kiad)}',
                           hist_xlabel='Cukorbetegségre fordított támogatás (egy emberre, USD)',
                           hist_filename='cukorbeteg_egeszs_kiad_egyemberre_USD_2024_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_cuk_egeszs_kiad['Log_Value'],
                           f'Cukorbetegséggel kapcsolatos egészségügyi kiadások log transzformáció utáni dobozábrája - {ev_felirat(df_cuk_egeszs_kiad)}',
                           'Cukorbetegséggel kapcsolatos egészségügyi kiadások (USD) log transzformálva',
                           'cukorbeteg_egeszsegugyi_kiadasok_dollar_2024_log_transz_boxplot.png',
                           hist_title=f'Cukorbetegséggel kapcsolatos egészségügyi kiadások USD/fő - {ev_felirat(df_cuk_egeszs_kiad)} log transzformáció után',
                           hist_xlabel='Cukorbetegségre fordított támogatás (egy főre, USD, log)',
                           hist_filename='cukorbeteg_egeszs_kiad_egyemberre_USD_2024_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# NEGYEDIK VIZSGÁLAT - cukorbetegségre visszavezethető halálozások száma országonként 2024

cuk_hal_2024 = pd.read_csv('cukorbetegsegre_visszavezetheto_halalok_szama_regionkent_orszagonkent_2024.csv', sep=',')

df_cuk_hal = vizsgalat('cukorbeteg_halalesetek',
                       who_hosszu(cuk_hal_2024, ev=2024),
                       per_fo=100_000)
df_cuk_hal.to_excel('cukorbetegsegre_visszavezetheto_halalozasi_szamok_2024.xlsx', index=False)
reszletes_jelentes(df_cuk_hal, '4. vizsgálat – halálesetek (nyers: fő; arányos: 100 000 lakosra)')

# Nyers értékek
abrak_keszitese_es_mentese(df_cuk_hal['Value'],
                           f'Cukorbetegségre visszavezethető halálesetek számának dobozábrája - {ev_felirat(df_cuk_hal)}',
                           'Cukorbetegségre visszavezethető halálesetek száma',
                           'cukorbetegsegre_visszavezetheto_halalesetszamok_2024_boxplot.png',
                           hist_title=f'Cukorbetegségre visszavezethető halálesetek számának eloszlása - {ev_felirat(df_cuk_hal)}',
                           hist_xlabel='Cukorbetegségre visszavezethető halálesetek száma',
                           hist_filename='cukorbetegsegre_visszavezetheto_halalesetszamok_2024_hist.png')

# Népességarányosítás után
abrak_keszitese_es_mentese(df_cuk_hal['Aranyos'],
                           f'Cukorbetegségre visszavezethető halálesetek 100 000 lakosra - {ev_felirat(df_cuk_hal)}',
                           'Halálesetek 100 000 lakosra',
                           'cukorbetegsegre_visszavezetheto_halalesetek_nepessegaranyos_2024_boxplot.png',
                           hist_title=f'Cukorbetegségre visszavezethető halálesetek népességarányos eloszlása - {ev_felirat(df_cuk_hal)}',
                           hist_xlabel='Halálesetek 100 000 lakosra',
                           hist_filename='cukorbetegsegre_visszavezetheto_halalesetek_nepessegaranyos_2024_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_cuk_hal['Log_Value'],
                           f'Cukorbetegségre visszavezethető halálesetek log transzformáció után - {ev_felirat(df_cuk_hal)}',
                           'Halálesetek 100 000 lakosra (log)',
                           'cukorbetegsegre_visszavezetheto_halalesetek_log_2024_boxplot.png',
                           hist_title=f'Cukorbetegségre visszavezethető halálesetek eloszlása log transzformáció után - {ev_felirat(df_cuk_hal)}',
                           hist_xlabel='Halálesetek 100 000 lakosra (log)',
                           hist_filename='cukorbetegsegre_visszavezetheto_halalesetek_log_2024_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ÖTÖDIK VIZSGÁLAT - egyes típusú cukorbetegek (0-19 évesek) száma ezer főben kifejezve 2024

fiatalkori_cuk_ezer_fo = pd.read_csv('egyes_tipusu_cukorbetegek_0-19_regionkent_orszagonkent_ezer_fo_2024.csv', sep=',')

df_elemezheto_fiatal = vizsgalat('egyes_tipusu_fiatalok',
                                 who_hosszu(fiatalkori_cuk_ezer_fo, ev=2024),
                                 egyseg_szorzo=1, per_fo=100_000)
df_elemezheto_fiatal.to_excel('egyes_tipusu_fiatalok_2024.xlsx', index=False)
reszletes_jelentes(df_elemezheto_fiatal, '5. vizsgálat – 1-es típusú cukorbetegek 0–19 év (arányos: 100 000 lakosra)')
sav_arany(df_elemezheto_fiatal, 'Aranyos', felso=40, cimke='40 alatti (100 000 lakosra)')

# A teljes népességgel osztunk (0-19 éves népességadat nincs). A 100 000 lakosra vetítés azért kell, mert %-ban
# nagyon kicsik az értékek, és a log1p ott alig változtat.

# Nyers értékek
abrak_keszitese_es_mentese(df_elemezheto_fiatal['Value'],
                           f'Egyes típusú cukorbetegek számának dobozábrája - {ev_felirat(df_elemezheto_fiatal)} (0-19 évesek körében)',
                           'Egyes típusú cukorbeteg fiatalkorúak száma (ezer fő)',
                           'egyes_tipusu_fiatalok_2024_boxplot.png',
                           hist_title=f'Egyes típusú cukorbeteg fiatalkorúak eloszlása - {ev_felirat(df_elemezheto_fiatal)}',
                           hist_xlabel='Egyes típusú cukorbeteg fiatalkorúak száma (ezer fő)',
                           hist_filename='egyes_tipusu_fiatalok_2024_hist.png')

# Népességarányosítás után
abrak_keszitese_es_mentese(df_elemezheto_fiatal['Aranyos'],
                           f'Egyes típusú cukorbetegek számának népességarányosítás utáni dobozábrája - {ev_felirat(df_elemezheto_fiatal)} (0-19 évesek körében)',
                           'Egyes típusú cukorbeteg fiatalkorúak 100 000 lakosra',
                           'egyes_tipusu_fiatalok_nepessegaranyositott_2024_boxplot.png',
                           hist_title=f'Egyes típusú cukorbeteg fiatalkorúak eloszlása népességarányosítás után - {ev_felirat(df_elemezheto_fiatal)}',
                           hist_xlabel='Egyes típusú cukorbeteg fiatalkorúak 100 000 lakosra',
                           hist_filename='egyes_tipusu_fiatalok_nepessegaranyos_2024_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_elemezheto_fiatal['Log_Value'],
                           f'Egyes típusú cukorbetegek számának log transzformált népességarányosítás utáni dobozábrája - {ev_felirat(df_elemezheto_fiatal)} (0-19 évesek körében)',
                           'Egyes típusú cukorbeteg fiatalkorúak népességaránya (log)',
                           'egyes_tipusu_fiatalok_log_2024_boxplot.png',
                           hist_title=f'Egyes típusú cukorbeteg fiatalkorúak eloszlása log transzformáció után - {ev_felirat(df_elemezheto_fiatal)}',
                           hist_xlabel='Egyes típusú cukorbeteg fiatalkorúak népességaránya (log)',
                           hist_filename='egyes_tipusu_fiatalok_log_2024_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# HATODIK VIZSGÁLAT - egyes típusú cukorbetegek száma az összes korosztályban összesen 2024

egyes_tip_minden_korosztaly_2024 = pd.read_csv('egyes_tipusu_cukorbetegek_minden_korosztalyban_regionkent_orszagonkent_2024.csv', sep=',')

df_elemezheto_minden_korosztaly = vizsgalat('egyes_tipusu_minden_korosztaly',
                                            who_hosszu(egyes_tip_minden_korosztaly_2024, ev=2024),
                                            egyseg_szorzo=1, per_fo=100_000)
df_elemezheto_minden_korosztaly.to_excel('egyes_tipusu_minden_korosztaly_2024.xlsx', index=False)
reszletes_jelentes(df_elemezheto_minden_korosztaly, '6. vizsgálat – 1-es típusú cukorbetegek, minden korosztály (arányos: 100 000 lakosra)')

# Nyers értékek
abrak_keszitese_es_mentese(df_elemezheto_minden_korosztaly['Value'],
                           f'Egyes típusú cukorbetegek számának dobozábrája - {ev_felirat(df_elemezheto_minden_korosztaly)}',
                           'Egyes típusú cukorbetegek száma (ezer fő)',
                           'egyes_tipusu_2024_boxplot.png',
                           hist_title=f'Egyes típusú cukorbetegek eloszlása - {ev_felirat(df_elemezheto_minden_korosztaly)}',
                           hist_xlabel='Egyes típusú cukorbetegek száma (ezer fő)',
                           hist_filename='egyes_tipusu_2024_hist.png')

# Népességarányosítás után
abrak_keszitese_es_mentese(df_elemezheto_minden_korosztaly['Aranyos'],
                           f'Egyes típusú cukorbetegek számának népességarányosítás utáni dobozábrája - {ev_felirat(df_elemezheto_minden_korosztaly)}',
                           'Egyes típusú cukorbetegek 100 000 lakosra',
                           'egyes_tipusu_nepessegaranyositott_2024_boxplot.png',
                           hist_title=f'Egyes típusú cukorbetegek eloszlása népességarányosítás után - {ev_felirat(df_elemezheto_minden_korosztaly)}',
                           hist_xlabel='Egyes típusú cukorbetegek 100 000 lakosra',
                           hist_filename='egyes_tipusu_nepessegaranyos_2024_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_elemezheto_minden_korosztaly['Log_Value'],
                           f'Egyes típusú cukorbetegek számának log transzformált népességarányosítás utáni dobozábrája - {ev_felirat(df_elemezheto_minden_korosztaly)}',
                           'Egyes típusú cukorbetegek népességaránya (log)',
                           'egyes_tipusu_log_2024_boxplot.png',
                           hist_title=f'Egyes típusú cukorbetegek eloszlása log transzformáció után - {ev_felirat(df_elemezheto_minden_korosztaly)}',
                           hist_xlabel='Egyes típusú cukorbetegek népességaránya (log)',
                           hist_filename='egyes_tipusu_log_2024_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# HETEDIK VIZSGÁLAT - felnőttek elhízási arányai (BMI > 30) 2024

felnottek_elhizas_szazalekban = pd.read_csv('felnottek_kozotti_elhizas_szazalekban_bmi_nagyobb_30.csv', sep=',')

df_felnottek_elhizas = vizsgalat('felnottkori_elhizas',
                                 who_hosszu(felnottek_elhizas_szazalekban, dim1='Both sexes'),
                                 aranyositas=False)
df_felnottek_elhizas.to_excel('felnottkori_elhizas_2024.xlsx', index=False)
reszletes_jelentes(df_felnottek_elhizas, '7. vizsgálat – felnőttkori elhízás (%)', 1)
sav_arany(df_felnottek_elhizas, 'Value', 10, 40, cimke='10–40% közötti')
sav_arany(df_felnottek_elhizas, 'Value', also=60, cimke='60% fölötti (csendes-óceáni csoport)')

# Nyers értékek
abrak_keszitese_es_mentese(df_felnottek_elhizas['Value'],
                           f'Felnőttkori elhízás népességarányának dobozábrája - {ev_felirat(df_felnottek_elhizas)}',
                           'Felnőttkori elhízás népességaránya (%)',
                           'felnottkori_elhizas_boxplot.png',
                           hist_title=f'Felnőttkori elhízás népességarányának eloszlása országok szintjén - {ev_felirat(df_felnottek_elhizas)}',
                           hist_xlabel='Felnőttkori elhízás népességaránya (%)',
                           hist_filename='felnottkori_elhizas_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_felnottek_elhizas['Log_Value'],
                           f'Felnőttkori elhízás népességarányának dobozábrája log transzformáció után - {ev_felirat(df_felnottek_elhizas)}',
                           'Felnőttkori elhízás népességaránya (log)',
                           'felnottkori_elhizas_log_boxplot.png',
                           hist_title=f'Felnőttkori elhízás népességarányának eloszlása log transzformáció után - {ev_felirat(df_felnottek_elhizas)}',
                           hist_xlabel='Felnőttkori elhízás népességaránya (log)',
                           hist_filename='felnottkori_elhizas_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# NYOLCADIK VIZSGÁLAT - GDP-arányos egészségügyi kiadások

gdp_aranyos_egeszs_kiad = load_worldbank_csv('gdp_aranyos_egeszsegugyi_kiadasok.csv')

df_gdp_aranyos_egeszs_kiad = vizsgalat('gdp_aranyos_egeszs_kiad',
                                       wb_hosszu(gdp_aranyos_egeszs_kiad),
                                       aranyositas=False)
df_gdp_aranyos_egeszs_kiad.to_excel('gdp_aranyos_egeszsegugyi_kiadasok.xlsx', index=False)
reszletes_jelentes(df_gdp_aranyos_egeszs_kiad, '8. vizsgálat – GDP-arányos egészségügyi kiadás (%)')
sav_arany(df_gdp_aranyos_egeszs_kiad, 'Value', 3, 10, cimke='3–10% közötti')

# Nyers értékek
abrak_keszitese_es_mentese(df_gdp_aranyos_egeszs_kiad['Value'],
                           f'GDP-arányos egészségügyi kiadások dobozábrája - {ev_felirat(df_gdp_aranyos_egeszs_kiad)}',
                           'Egészségügyi kiadás (a GDP %-a)',
                           'gdp_aranyos_egeszsegugyi_kiadasok_boxplot.png',
                           hist_title=f'GDP-arányos egészségügyi kiadások eloszlása - {ev_felirat(df_gdp_aranyos_egeszs_kiad)}',
                           hist_xlabel='Egészségügyi kiadás (a GDP %-a)',
                           hist_filename='gdp_aranyos_egeszsegugyi_kiadasok_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_gdp_aranyos_egeszs_kiad['Log_Value'],
                           f'GDP-arányos egészségügyi kiadások dobozábrája log transzformáció után - {ev_felirat(df_gdp_aranyos_egeszs_kiad)}',
                           'Egészségügyi kiadás (a GDP %-a, log)',
                           'gdp_aranyos_egeszsegugyi_kiadasok_log_boxplot.png',
                           hist_title=f'GDP-arányos egészségügyi kiadások eloszlása log transzformáció után - {ev_felirat(df_gdp_aranyos_egeszs_kiad)}',
                           hist_xlabel='Egészségügyi kiadás (a GDP %-a, log)',
                           hist_filename='gdp_aranyos_egeszsegugyi_kiadasok_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# KILENCEDIK VIZSGÁLAT - egy főre jutó GDP

gdp_per_capita = load_worldbank_csv('gdp_per_capita.csv')

df_gdp_per_capita = vizsgalat('gdp_per_capita',
                              wb_hosszu(gdp_per_capita),
                              aranyositas=False)
df_gdp_per_capita.to_excel('gdp_per_capita.xlsx', index=False)
reszletes_jelentes(df_gdp_per_capita, '9. vizsgálat – egy főre jutó GDP (PPP)', 0)

# Nyers értékek
abrak_keszitese_es_mentese(df_gdp_per_capita['Value'],
                           f'Egy főre jutó GDP dobozábrája - {ev_felirat(df_gdp_per_capita)}',
                           'Egy főre jutó GDP (USD)',
                           'gdp_per_capita_boxplot.png',
                           hist_title=f'Egy főre jutó GDP eloszlása - {ev_felirat(df_gdp_per_capita)}',
                           hist_xlabel='Egy főre jutó GDP (USD)',
                           hist_filename='gdp_per_capita_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_gdp_per_capita['Log_Value'],
                           f'Egy főre jutó GDP dobozábrája log transzformáció után - {ev_felirat(df_gdp_per_capita)}',
                           'Egy főre jutó GDP (USD, log)',
                           'gdp_per_capita_log_boxplot.png',
                           hist_title=f'Egy főre jutó GDP eloszlása log transzformáció után - {ev_felirat(df_gdp_per_capita)}',
                           hist_xlabel='Egy főre jutó GDP (USD, log)',
                           hist_filename='gdp_per_capita_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIZEDIK VIZSGÁLAT - gondozók száma ezer emberre

gondozok_szama = load_worldbank_csv('gondozok_szama.csv')

df_gondozok_szama = vizsgalat('gondozok_szama',
                              wb_hosszu(gondozok_szama),
                              aranyositas=False)
df_gondozok_szama.to_excel('gondozok_szama.xlsx', index=False)
reszletes_jelentes(df_gondozok_szama, '10. vizsgálat – gondozók 1000 lakosra')

# Nyers értékek
abrak_keszitese_es_mentese(df_gondozok_szama['Value'],
                           f'Gondozók száma ezer lakosra - dobozábra - {ev_felirat(df_gondozok_szama)}',
                           'Gondozók 1000 lakosra',
                           'gondozok_szama_boxplot.png',
                           hist_title=f'Gondozók számának eloszlása ezer lakosra - {ev_felirat(df_gondozok_szama)}',
                           hist_xlabel='Gondozók 1000 lakosra',
                           hist_filename='gondozok_szama_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_gondozok_szama['Log_Value'],
                           f'Gondozók száma ezer lakosra log transzformáció után - dobozábra - {ev_felirat(df_gondozok_szama)}',
                           'Gondozók 1000 lakosra (log)',
                           'gondozok_szama_log_boxplot.png',
                           hist_title=f'Gondozók számának eloszlása log transzformáció után - {ev_felirat(df_gondozok_szama)}',
                           hist_xlabel='Gondozók 1000 lakosra (log)',
                           hist_filename='gondozok_szama_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIZENEGYEDIK VIZSGÁLAT - kórházi ágyak száma

korhazi_agyak = load_worldbank_csv('korhazi_agyak_szama.csv')

df_korhazi_agyak = vizsgalat('korhazi_agyak',
                             wb_hosszu(korhazi_agyak),
                             aranyositas=False)
df_korhazi_agyak.to_excel('korhazi_agyak_szama.xlsx', index=False)
reszletes_jelentes(df_korhazi_agyak, '11. vizsgálat – kórházi ágyak 1000 lakosra')

# Nyers értékek
abrak_keszitese_es_mentese(df_korhazi_agyak['Value'],
                           f'Kórházi ágyak száma ezer lakosra - dobozábra - {ev_felirat(df_korhazi_agyak)}',
                           'Kórházi ágyak 1000 lakosra',
                           'korhazi_agyak_boxplot.png',
                           hist_title=f'Kórházi ágyak számának eloszlása ezer lakosra - {ev_felirat(df_korhazi_agyak)}',
                           hist_xlabel='Kórházi ágyak 1000 lakosra',
                           hist_filename='korhazi_agyak_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_korhazi_agyak['Log_Value'],
                           f'Kórházi ágyak száma log transzformáció után - dobozábra - {ev_felirat(df_korhazi_agyak)}',
                           'Kórházi ágyak 1000 lakosra (log)',
                           'korhazi_agyak_log_boxplot.png',
                           hist_title=f'Kórházi ágyak számának eloszlása log transzformáció után - {ev_felirat(df_korhazi_agyak)}',
                           hist_xlabel='Kórházi ágyak 1000 lakosra (log)',
                           hist_filename='korhazi_agyak_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIZENKETTEDIK VIZSGÁLAT - korstandardizált cukorbetegek aránya 2024

korstand_cuk_ar = pd.read_csv("korstandardizalt_cukorbetegek_aranya_orszagonkent_2024.csv", sep=',')

df_korstand_cuk_ar = vizsgalat('korstand_cukorbeteg_arany',
                               who_hosszu(korstand_cuk_ar, ev=2024),
                               aranyositas=False)
df_korstand_cuk_ar.to_excel('korstandardizalt_cukorbetegek_aranya_2024.xlsx', index=False)
reszletes_jelentes(df_korstand_cuk_ar, '12. vizsgálat – korstandardizált cukorbeteg-arány (%)', 1)

# Nyers értékek
abrak_keszitese_es_mentese(df_korstand_cuk_ar['Value'],
                           f'Korstandardizált cukorbetegek aránya - dobozábra - {ev_felirat(df_korstand_cuk_ar)}',
                           'Korstandardizált cukorbetegek aránya (%)',
                           'korstandardizalt_cukorbetegek_aranya_2024_boxplot.png',
                           hist_title=f'Korstandardizált cukorbetegek arányának eloszlása - {ev_felirat(df_korstand_cuk_ar)}',
                           hist_xlabel='Korstandardizált cukorbetegek aránya (%)',
                           hist_filename='korstandardizalt_cukorbetegek_aranya_2024_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_korstand_cuk_ar['Log_Value'],
                           f'Korstandardizált cukorbetegek aránya log transzformáció után - dobozábra - {ev_felirat(df_korstand_cuk_ar)}',
                           'Korstandardizált cukorbetegek aránya (log)',
                           'korstandardizalt_cukorbetegek_aranya_log_2024_boxplot.png',
                           hist_title=f'Korstandardizált cukorbetegek arányának eloszlása log transzformáció után - {ev_felirat(df_korstand_cuk_ar)}',
                           hist_xlabel='Korstandardizált cukorbetegek aránya (log)',
                           hist_filename='korstandardizalt_cukorbetegek_aranya_log_2024_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIZENHARMADIK VIZSGÁLAT - orvosok száma

medical_doc_num = load_hibas_who_csv('medical_doctors_number.csv')

df_orvosok_szama = vizsgalat('orvosok_szama',
                             who_hosszu(medical_doc_num),
                             per_fo=10_000)
df_orvosok_szama.to_excel('orvosok_szama.xlsx', index=False)
reszletes_jelentes(df_orvosok_szama, '13. vizsgálat – orvosok (nyers: fő; arányos: 10 000 lakosra)')
print('   5 alatti országok:', sav_arany(df_orvosok_szama, 'Aranyos', felso=5, cimke='5 alatti (10 000 lakosra)'))
sav_arany(df_orvosok_szama, 'Aranyos', 25, 45, cimke='25–45 közötti (10 000 lakosra)')

# 10 000 lakosra vetítve, így összevethető a 14. vizsgálat WHO-féle mutatójával.

# Nyers értékek
abrak_keszitese_es_mentese(df_orvosok_szama['Value'],
                           f'Orvosok száma - dobozábra - {ev_felirat(df_orvosok_szama)}',
                           'Orvosok száma (fő)',
                           'orvosok_szama_boxplot.png',
                           hist_title=f'Orvosok számának eloszlása - {ev_felirat(df_orvosok_szama)}',
                           hist_xlabel='Orvosok száma (fő)',
                           hist_filename='orvosok_szama_hist.png')

# Népességarányosítás után
abrak_keszitese_es_mentese(df_orvosok_szama['Aranyos'],
                           f'Orvosok száma 10 000 lakosra - dobozábra - {ev_felirat(df_orvosok_szama)}',
                           'Orvosok 10 000 lakosra',
                           'orvosok_szama_nepessegaranyos_boxplot.png',
                           hist_title=f'Orvosok számának népességarányos eloszlása - {ev_felirat(df_orvosok_szama)}',
                           hist_xlabel='Orvosok 10 000 lakosra',
                           hist_filename='orvosok_szama_nepessegaranyos_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_orvosok_szama['Log_Value'],
                           f'Orvosok száma 10 000 lakosra log transzformáció után - dobozábra - {ev_felirat(df_orvosok_szama)}',
                           'Orvosok 10 000 lakosra (log)',
                           'orvosok_szama_log_boxplot.png',
                           hist_title=f'Orvosok számának népességarányos eloszlása log transzformáció után - {ev_felirat(df_orvosok_szama)}',
                           hist_xlabel='Orvosok 10 000 lakosra (log)',
                           hist_filename='orvosok_szama_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIZENNEGYEDIK VIZSGÁLAT - orvosok száma 10 000 főre

orvosok_tizezer = load_hibas_who_csv('orvosok_szama_tizezer_fore.csv')

df_orvosok_tizezer = vizsgalat('orvosok_tizezer_fore',
                               who_hosszu(orvosok_tizezer),
                               aranyositas=False)
df_orvosok_tizezer.to_excel('orvosok_szama_tizezer_fore.xlsx', index=False)
reszletes_jelentes(df_orvosok_tizezer, '14. vizsgálat – orvosok 10 000 lakosra')
sav_arany(df_orvosok_tizezer, 'Value', felso=5, cimke='5 alatti')
sav_arany(df_orvosok_tizezer, 'Value', 25, 45, cimke='25–45 közötti')

# Nyers értékek
abrak_keszitese_es_mentese(df_orvosok_tizezer['Value'],
                           f'Orvosok száma 10 000 főre - dobozábra - {ev_felirat(df_orvosok_tizezer)}',
                           'Orvosok 10 000 lakosra',
                           'orvosok_tizezer_fore_boxplot.png',
                           hist_title=f'Orvosok 10 000 főre jutó számának eloszlása - {ev_felirat(df_orvosok_tizezer)}',
                           hist_xlabel='Orvosok 10 000 lakosra',
                           hist_filename='orvosok_tizezer_fore_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_orvosok_tizezer['Log_Value'],
                           f'Orvosok száma 10 000 főre log transzformáció után - dobozábra - {ev_felirat(df_orvosok_tizezer)}',
                           'Orvosok 10 000 lakosra (log)',
                           'orvosok_tizezer_fore_log_boxplot.png',
                           hist_title=f'Orvosok 10 000 főre jutó számának eloszlása log transzformáció után - {ev_felirat(df_orvosok_tizezer)}',
                           hist_xlabel='Orvosok 10 000 lakosra (log)',
                           hist_filename='orvosok_tizezer_fore_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIZENÖTÖDIK VIZSGÁLAT - teljes cukorbetegséggel kapcsolatos kiadások USD 2024

cuk_kiad_2024_usd = pd.read_csv('teljes_cukorbetegseggel_kapcsolatos_egeszsegugyi_kiadasok_USD_regionkent_orszagonkent_2024.csv', sep=',')

df_cuk_kiad_teljes = vizsgalat('cuk_kiadas_teljes',
                               who_hosszu(cuk_kiad_2024_usd, ev=2024),
                               egyseg_szorzo=1_000_000)
df_cuk_kiad_teljes.to_excel('teljes_cukorbeteg_kiadasok_2024.xlsx', index=False)
reszletes_jelentes(df_cuk_kiad_teljes, '15. vizsgálat – teljes kiadás (nyers: millió USD; arányos: USD/lakos)')

# Nyers értékek
abrak_keszitese_es_mentese(df_cuk_kiad_teljes['Value'],
                           f'Cukorbetegséggel kapcsolatos teljes kiadások - dobozábra - {ev_felirat(df_cuk_kiad_teljes)}',
                           'Teljes kiadás (USD)',
                           'teljes_cukorbeteg_kiadasok_2024_boxplot.png',
                           hist_title=f'Cukorbetegséggel kapcsolatos teljes kiadások eloszlása - {ev_felirat(df_cuk_kiad_teljes)}',
                           hist_xlabel='Teljes kiadás (USD)',
                           hist_filename='teljes_cukorbeteg_kiadasok_2024_hist.png')

# Népességarányosítás után
abrak_keszitese_es_mentese(df_cuk_kiad_teljes['Aranyos'],
                           f'Cukorbetegséggel kapcsolatos kiadások egy lakosra - dobozábra - {ev_felirat(df_cuk_kiad_teljes)}',
                           'Kiadás egy lakosra (USD/fő)',
                           'teljes_cukorbeteg_kiadasok_nepessegaranyos_2024_boxplot.png',
                           hist_title=f'Cukorbetegséggel kapcsolatos kiadások egy lakosra jutó eloszlása - {ev_felirat(df_cuk_kiad_teljes)}',
                           hist_xlabel='Kiadás egy lakosra (USD/fő)',
                           hist_filename='teljes_cukorbeteg_kiadasok_nepessegaranyos_2024_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_cuk_kiad_teljes['Log_Value'],
                           f'Cukorbetegséggel kapcsolatos kiadások egy lakosra log transzformáció után - dobozábra - {ev_felirat(df_cuk_kiad_teljes)}',
                           'Kiadás egy lakosra (USD/fő, log)',
                           'teljes_cukorbeteg_kiadasok_log_2024_boxplot.png',
                           hist_title=f'Cukorbetegséggel kapcsolatos kiadások egy lakosra log transzformáció után - {ev_felirat(df_cuk_kiad_teljes)}',
                           hist_xlabel='Kiadás egy lakosra (USD/fő, log)',
                           hist_filename='teljes_cukorbeteg_kiadasok_log_2024_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIZENHATODIK VIZSGÁLAT - teljes populáció

df_teljes_pop = vizsgalat('teljes_populacio',
                          nepesseg_hosszu[['Orszag', 'Ev', 'Nepesseg']].rename(columns={'Nepesseg': 'Value'}),
                          aranyositas=False)
df_teljes_pop.to_excel('teljes_populacio.xlsx', index=False)
reszletes_jelentes(df_teljes_pop, '16. vizsgálat – teljes népesség (fő)', 0)
sav_arany(df_teljes_pop, 'Value', also=100_000_000, cimke='100 millió fölötti')
sav_arany(df_teljes_pop, 'Value', felso=1_000_000, cimke='1 millió alatti')

# Nyers értékek
abrak_keszitese_es_mentese(df_teljes_pop['Value'],
                           f'Teljes népesség - dobozábra - {ev_felirat(df_teljes_pop)}',
                           'Népesség (fő)',
                           'teljes_populacio_boxplot.png',
                           hist_title=f'Országok népességének eloszlása - {ev_felirat(df_teljes_pop)}',
                           hist_xlabel='Népesség (fő)',
                           hist_filename='teljes_populacio_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_teljes_pop['Log_Value'],
                           f'Teljes népesség log transzformáció után - dobozábra - {ev_felirat(df_teljes_pop)}',
                           'Népesség (log)',
                           'teljes_populacio_log_boxplot.png',
                           hist_title=f'Országok népességének eloszlása log transzformáció után - {ev_felirat(df_teljes_pop)}',
                           hist_xlabel='Népesség (log)',
                           hist_filename='teljes_populacio_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIZENHETEDIK VIZSGÁLAT - universal health coverage index

health_cov_idx = pd.read_csv('universal-health-coverage-index.csv', sep=',')

df_health_cov_idx = vizsgalat('health_coverage_index',
                              who_hosszu(health_cov_idx, orszag_oszlop='Entity', ev_oszlop='Year',
                                         ertek_oszlop='UHC Service Coverage Index (SDG 3.8.1)'),
                              aranyositas=False)
df_health_cov_idx.to_excel('universal_health_coverage_index.xlsx', index=False)
reszletes_jelentes(df_health_cov_idx, '17. vizsgálat – UHC index (pont)', 0)
print('   40 pont alatti országok:', sav_arany(df_health_cov_idx, 'Value', felso=40, cimke='40 pont alatti'))

# Pontskála (0-100), nem százalék.

# Nyers értékek
abrak_keszitese_es_mentese(df_health_cov_idx['Value'],
                           f'Egyetemes egészségügyi lefedettség (UHC) indexe - dobozábra - {ev_felirat(df_health_cov_idx)}',
                           'UHC index (pont)',
                           'health_coverage_index_boxplot.png',
                           hist_title=f'Egyetemes egészségügyi lefedettség (UHC) indexének eloszlása - {ev_felirat(df_health_cov_idx)}',
                           hist_xlabel='UHC index (pont)',
                           hist_filename='health_coverage_index_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_health_cov_idx['Log_Value'],
                           f'UHC index log transzformáció után - dobozábra - {ev_felirat(df_health_cov_idx)}',
                           'UHC index (log)',
                           'health_coverage_index_log_boxplot.png',
                           hist_title=f'UHC index eloszlása log transzformáció után - {ev_felirat(df_health_cov_idx)}',
                           hist_xlabel='UHC index (log)',
                           hist_filename='health_coverage_index_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIZENNYOLCADIK VIZSGÁLAT - városiasodás

varosiasodas = load_worldbank_csv('varosiasodas.csv')

df_varosiasodas = vizsgalat('varosiasodas',
                            wb_hosszu(varosiasodas),
                            aranyositas=False)
df_varosiasodas.to_excel('varosiasodas.xlsx', index=False)
reszletes_jelentes(df_varosiasodas, '18. vizsgálat – városi népesség aránya (%)')
print('   100%-os városi népességű országok:', sav_arany(df_varosiasodas, 'Value', also=100, cimke='100%'))

# Nyers értékek
abrak_keszitese_es_mentese(df_varosiasodas['Value'],
                           f'Városi népesség aránya - dobozábra - {ev_felirat(df_varosiasodas)}',
                           'Városi népesség aránya (%)',
                           'varosiasodas_boxplot.png',
                           hist_title=f'Városi népesség arányának eloszlása - {ev_felirat(df_varosiasodas)}',
                           hist_xlabel='Városi népesség aránya (%)',
                           hist_filename='varosiasodas_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_varosiasodas['Log_Value'],
                           f'Városi népesség aránya log transzformáció után - dobozábra - {ev_felirat(df_varosiasodas)}',
                           'Városi népesség aránya (log)',
                           'varosiasodas_log_boxplot.png',
                           hist_title=f'Városi népesség arányának eloszlása log transzformáció után - {ev_felirat(df_varosiasodas)}',
                           hist_xlabel='Városi népesség aránya (log)',
                           hist_filename='varosiasodas_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TIZENKILENCEDIK VIZSGÁLAT - Humán fejlettségi index (HDI) 2023

hdi_index = pd.read_excel('iskolazottsag_human_kepzettsegi_index.xlsx', sheet_name='HDI_olvashato')

df_hdi = vizsgalat('hdi_index',
                   who_hosszu(hdi_index, ev=2023, orszag_oszlop='Country',
                              ertek_oszlop='Human_Development_Index_2023'),
                   aranyositas=False)
df_hdi.to_excel('human_fejlettsegi_index_2023.xlsx', index=False)
reszletes_jelentes(df_hdi, '19. vizsgálat – HDI', 3)

# Nyers értékek
abrak_keszitese_es_mentese(df_hdi['Value'],
                           f'Humán fejlettségi index (HDI) dobozábrája - {ev_felirat(df_hdi)}',
                           'HDI (0-1)',
                           'hdi_index_boxplot.png',
                           hist_title=f'Humán fejlettségi index (HDI) eloszlása - {ev_felirat(df_hdi)}',
                           hist_xlabel='HDI (0-1)',
                           hist_filename='hdi_index_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_hdi['Log_Value'],
                           f'Humán fejlettségi index (HDI) dobozábrája log transzformáció után - {ev_felirat(df_hdi)}',
                           'HDI (log)',
                           'hdi_index_log_boxplot.png',
                           hist_title=f'Humán fejlettségi index (HDI) eloszlása log transzformáció után - {ev_felirat(df_hdi)}',
                           hist_xlabel='HDI (log)',
                           hist_filename='hdi_index_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# HUSZADIK VIZSGÁLAT - nem diagnosztizált cukorbetegek aránya a 20–79 éves cukorbetegek között 2024 (IDF)
# A K modell ellenőrzéséhez és a halálozási ellenőrzéshez. A fájlban nincs évoszlop, ez az IDF 2024-es becslése.

ALULDIAG_FAJL = 'nem_diagnosztizalt_cukorbetegek_eloszlasa.csv'
WHO_CEL_ALULDIAG = 20.0
# WHO-cél: a cukorbetegek legalább 80%-a diagnosztizált, vagyis legfeljebb 20% nem diagnosztizált.

nem_diagnosztizalt_arany = pd.read_csv(ALULDIAG_FAJL, sep=',')

df_aluldiag = vizsgalat('aluldiag_arany',
                        who_hosszu(nem_diagnosztizalt_arany, ev=2024),
                        aranyositas=False)
df_aluldiag.to_excel('nem_diagnosztizalt_cukorbetegek_aranya_2024.xlsx', index=False)
reszletes_jelentes(df_aluldiag, '20. vizsgálat – nem diagnosztizált cukorbetegek aránya (a cukorbetegek %-a)', 1)
print('   WHO-célt teljesítő országok:',
      sav_arany(df_aluldiag, 'Value', felso=WHO_CEL_ALULDIAG, cimke=f'a WHO-célt teljesítő ({WHO_CEL_ALULDIAG:.0f}% vagy kevesebb)'))
sav_arany(df_aluldiag, 'Value', also=50, cimke='50% fölötti (a cukorbetegek többsége nem tud a betegségéről)')
sav_arany(df_aluldiag, 'Value', also=75, cimke='75% fölötti (négyből legalább három cukorbeteg nem tud róla)')

_ism = df_aluldiag['Value'].value_counts()
_ism = _ism[_ism >= 3]
print(f'   Legalább 3 országnál azonos értékek ({len(_ism)} db, összesen {int(_ism.sum())} ország):')
for _ertek, _db in _ism.items():
    _orsz = df_aluldiag.loc[df_aluldiag['Value'] == _ertek, 'Country Name'].tolist()
    print(f'      {_szam(_ertek, 1)}% – {_db} ország: {", ".join(_orsz[:8])}{" …" if len(_orsz) > 8 else ""}')

# Nyers értékek (arányszám, népességarányosítás nincs)
abrak_keszitese_es_mentese(df_aluldiag['Value'],
                           f'Nem diagnosztizált cukorbetegek aránya - dobozábra - {ev_felirat(df_aluldiag)}',
                           'Nem diagnosztizált arány (a cukorbetegek %-a)',
                           'aluldiag_arany_boxplot.png',
                           hist_title=f'Nem diagnosztizált cukorbetegek arányának eloszlása - {ev_felirat(df_aluldiag)}',
                           hist_xlabel='Nem diagnosztizált arány (a cukorbetegek %-a)',
                           hist_filename='aluldiag_arany_hist.png')

# Log transzformáció után
abrak_keszitese_es_mentese(df_aluldiag['Log_Value'],
                           f'Nem diagnosztizált cukorbetegek aránya log transzformáció után - dobozábra - {ev_felirat(df_aluldiag)}',
                           'Nem diagnosztizált arány (log)',
                           'aluldiag_arany_log_boxplot.png',
                           hist_title=f'Nem diagnosztizált cukorbetegek arányának eloszlása log transzformáció után - {ev_felirat(df_aluldiag)}',
                           hist_xlabel='Nem diagnosztizált arány (log)',
                           hist_filename='aluldiag_arany_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# HUSZONEGYEDIK ÉS HUSZONKETTEDIK VIZSGÁLAT - demográfiai kontrollváltozók (Világbank)
# A halálozási ellenőrzésnél a korösszetétel és az általános halálozási szint kiszűrésére szolgálnak.

def demografiai_vizsgalat(nev, fajl, cim, tengely):
    df = vizsgalat(nev, wb_hosszu(load_worldbank_rugalmas(fajl)), aranyositas=False)
    df.to_excel(f'{nev}.xlsx', index=False)
    reszletes_jelentes(df, cim, 1)
    abrak_keszitese_es_mentese(df['Value'], f'{tengely} - dobozábra - {ev_felirat(df)}', tengely, f'{nev}_boxplot.png',
                               hist_title=f'{tengely} eloszlása - {ev_felirat(df)}', hist_xlabel=tengely,
                               hist_filename=f'{nev}_hist.png')
    abrak_keszitese_es_mentese(df['Log_Value'], f'{tengely} log transzformáció után - {ev_felirat(df)}', f'{tengely} (log)',
                               f'{nev}_log_boxplot.png',
                               hist_title=f'{tengely} eloszlása log transzformáció után - {ev_felirat(df)}',
                               hist_xlabel=f'{tengely} (log)', hist_filename=f'{nev}_log_hist.png')
    return df


# 21. vizsgálat: nyers halálozási arány (összes haláleset 1000 lakosra) – az általános halálozási szint
df_nyers_halalozas = demografiai_vizsgalat('nyers_halalozasi_arany', 'nyers_halalozasi_arany_1000_lakosra_orszagonkent.csv',
                                           '21. vizsgálat – nyers halálozási arány (1000 lakosra)',
                                           'Nyers halálozási arány (1000 lakosra)')

# 22. vizsgálat: 65 év felettiek aránya a népességben (%) – a korösszetétel
df_idosek = demografiai_vizsgalat('idosek_aranya', '65_ev_felettiek_aranya_a_nepessegben_orszagonkent.csv',
                                  '22. vizsgálat – 65 év felettiek aránya a népességben (%)',
                                  '65 év felettiek aránya (%)')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# HUSZONHARMADIK VIZSGÁLAT - kezeletlen cukorbetegek aránya (NCD-RisC, korstandardizált, 30+ évesek) - a FŐ CÉLVÁLTOZÓ (K modell)
# A kezeletlen arány = 100 − a gyógyszeresen kezelt cukorbetegek aránya. Aki nincs diagnosztizálva, az kezelést sem kap,
# így ez a mutató a felismerés hiányát és a diagnózis utáni ellátási hiányt együtt méri ("cascade of care" első két lépése).
# Előnye az IDF nem diagnosztizált arányával szemben: független forrásból (NCD-RisC, vérvizsgálatos felmérések, saját modell)
# származik, és kb. 200 országra ad becslést. A legfrissebb év 2022.

ncdrisc = load_ncdrisc('kezelt_cukorbetegek_aranya_ncdrisc.csv')

_kez = ncdrisc[['Orszag', 'Ev']].copy()
_kez['Value'] = 100 - ncdrisc['Value']
_kez.attrs['ev_forras'] = ncdrisc.attrs['ev_forras']
df_kezeletlen = vizsgalat('kezeletlen_arany', _kez, aranyositas=False)
df_kezeletlen.to_excel('kezeletlen_cukorbetegek_aranya.xlsx', index=False)
reszletes_jelentes(df_kezeletlen, '23. vizsgálat – kezeletlen cukorbetegek aránya (a cukorbetegek %-a, NCD-RisC)', 1)
sav_arany(df_kezeletlen, 'Value', also=50, cimke='50% fölötti (a cukorbetegek többsége nem kap kezelést)')
sav_arany(df_kezeletlen, 'Value', also=75, cimke='75% fölötti')

_ev_K = int(df_kezeletlen['Ev'].iloc[0])
kezeles_ci = ncdrisc[ncdrisc['Ev'] == _ev_K].assign(Match_Key=lambda d: illesztesi_kulcs(d['Orszag']))
kezeles_ci['CI_szelesseg'] = kezeles_ci['CI_felso'] - kezeles_ci['CI_also']
print(f'   A kezelési arány 95%-os bizonytalansági intervallumának átlagos szélessége: {_szam(kezeles_ci["CI_szelesseg"].mean(), 1)} százalékpont')

abrak_keszitese_es_mentese(df_kezeletlen['Value'],
                           f'Kezeletlen cukorbetegek aránya - dobozábra - {ev_felirat(df_kezeletlen)}',
                           'Kezeletlen arány (a cukorbetegek %-a)', 'kezeletlen_arany_boxplot.png',
                           hist_title=f'Kezeletlen cukorbetegek arányának eloszlása - {ev_felirat(df_kezeletlen)}',
                           hist_xlabel='Kezeletlen arány (a cukorbetegek %-a)', hist_filename='kezeletlen_arany_hist.png')
abrak_keszitese_es_mentese(df_kezeletlen['Log_Value'],
                           f'Kezeletlen cukorbetegek aránya log transzformáció után - {ev_felirat(df_kezeletlen)}',
                           'Kezeletlen arány (log)', 'kezeletlen_arany_log_boxplot.png',
                           hist_title=f'Kezeletlen cukorbetegek arányának eloszlása log transzformáció után - {ev_felirat(df_kezeletlen)}',
                           hist_xlabel='Kezeletlen arány (log)', hist_filename='kezeletlen_arany_log_hist.png')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ÖSSZESÍTÉS

# Ellenőrzés
print('Az összesítőben szereplő vizsgálatok és lépések:')
for nev in dict.fromkeys(k[0] for k in osszesito):
    print(f'  {nev}: {[k[1] for k in osszesito if k[0] == nev]}')

# Minden leíró statisztika egy táblában:
osszesito_tabla = pd.DataFrame(osszesito).T
osszesito_tabla.index.names = ['Vizsgálat', 'Lépés']
osszesito_tabla = osszesito_tabla.reset_index().infer_objects()  # infer_objects: a számoszlopok újra számként
print(osszesito_tabla[['Vizsgálat', 'Lépés', 'Felhasznált év', 'Elemszám', 'Ferdeség (Skewness)', 'Csúcsosság (Kurtosis)']].round(3))


# JAVASOLT MODELLVÁLTOZAT VIZSGÁLATONKÉNT

LEPES_OSZLOP = {'nyers': 'Value', 'nepessegaranyos': 'Aranyos', 'log': 'Log_Value'}
LEPES_ELOTAG = {'nyers': 'nyers', 'nepessegaranyos': 'arany', 'log': 'log'}

def alak_eltérés(stat):
    return abs(stat['Ferdeség (Skewness)']) + abs(stat['Csúcsosság (Kurtosis)'])

def elfogadhato(stat):
    return abs(stat['Ferdeség (Skewness)']) < NORMALITASI_KUSZOB and abs(stat['Csúcsosság (Kurtosis)']) < NORMALITASI_KUSZOB

def javasolt_valtozat(nev):
    lepesek = {k[1]: osszesito[k] for k in osszesito if k[0] == nev}
    alap = 'nepessegaranyos' if 'nepessegaranyos' in lepesek else 'nyers'
    if elfogadhato(lepesek[alap]):
        return alap, 'az alapváltozat eloszlása már elfogadható, log transzformáció nem szükséges'
    if alak_eltérés(lepesek['log']) < alak_eltérés(lepesek[alap]):
        megj = 'log transzformáció szükséges és javít'
        if not elfogadhato(lepesek['log']):
            megj += ' (de a küszöböt log után sem éri el teljesen)'
        return 'log', megj
    return alap, 'a log transzformáció nem javít az eloszláson, ezért az alapváltozat marad'

javaslatok = {nev: javasolt_valtozat(nev) for nev in dict.fromkeys(k[0] for k in osszesito)}

# Második munkalap
lepes_sorrend = ['nyers', 'nepessegaranyos', 'log']
alak_osszevetes = osszesito_tabla.pivot(index='Vizsgálat', columns='Lépés',
                                        values=['Ferdeség (Skewness)', 'Csúcsosság (Kurtosis)'])
alak_osszevetes = alak_osszevetes.reindex(columns=lepes_sorrend, level=1)
alak_osszevetes.columns = [f'{mutato} - {lepes}' for mutato, lepes in alak_osszevetes.columns]
alak_osszevetes = alak_osszevetes.reindex(list(dict.fromkeys(osszesito_tabla['Vizsgálat']))).reset_index()
ev_oszlopok = osszesito_tabla.drop_duplicates('Vizsgálat')[
    ['Vizsgálat', 'Felhasznált év', 'Lefedettség a választott évben', 'Lefedettségi megjegyzés', 'Az év forrása']]
alak_osszevetes = ev_oszlopok.merge(alak_osszevetes, on='Vizsgálat', how='right')
alak_osszevetes['Népességarányosítás'] = alak_osszevetes['Vizsgálat'].map(
    lambda n: 'igen' if (n, 'nepessegaranyos') in osszesito else 'nem (eleve arány / egy főre jutó / index)')
alak_osszevetes['Javasolt modellváltozat'] = alak_osszevetes['Vizsgálat'].map(lambda n: javaslatok[n][0])
alak_osszevetes['Indoklás'] = alak_osszevetes['Vizsgálat'].map(lambda n: javaslatok[n][1])

with pd.ExcelWriter('leiro_statisztikak_osszesito.xlsx') as writer:
    osszesito_tabla.to_excel(writer, sheet_name='Leiro_statisztikak', index=False)
    alak_osszevetes.to_excel(writer, sheet_name='Ferdeseg_csucsossag', index=False)

print(alak_osszevetes[['Vizsgálat', 'Felhasznált év', 'Népességarányosítás', 'Javasolt modellváltozat']].to_string(index=False))


# MODELLBEMENET

modell_df = None
modell_log_df = None
modell_evek = []
for nev, df in eredmenyek.items():
    lepes = javaslatok[nev][0]
    oszlopnev = f'{LEPES_ELOTAG[lepes]}_{nev}'
    resz = df[['Country Name', LEPES_OSZLOP[lepes]]].rename(columns={LEPES_OSZLOP[lepes]: oszlopnev})
    resz_log = df[['Country Name', 'Log_Value']].rename(columns={'Log_Value': f'log_{nev}'})
    modell_df = resz if modell_df is None else modell_df.merge(resz, on='Country Name', how='outer')
    modell_log_df = resz_log if modell_log_df is None else modell_log_df.merge(resz_log, on='Country Name', how='outer')
    modell_evek.append({'Változó': oszlopnev, 'Vizsgálat': nev, 'Felhasznált év': int(df['Ev'].iloc[0]),
                        'Elemszám': len(df), 'Javasolt változat': lepes, 'Indoklás': javaslatok[nev][1]})

modell_df = modell_df.sort_values('Country Name')
modell_log_df = modell_log_df.sort_values('Country Name')
modell_evek = pd.DataFrame(modell_evek)

print(modell_df.shape)
print(modell_df.notna().sum())
print(modell_evek[['Változó', 'Felhasznált év', 'Elemszám']].to_string(index=False))
with pd.ExcelWriter('modell_adatok.xlsx') as writer:
    modell_df.to_excel(writer, sheet_name='Modell_bemenet', index=False)
    modell_evek.to_excel(writer, sheet_name='Valtozok_es_evek', index=False)
    modell_log_df.to_excel(writer, sheet_name='Log_ertekek', index=False)

# =========================================================================================================================================================================
# REGRESSZIÓS MODELLEK
# =========================================================================================================================================================================
#
# Minden modell országszintű, keresztmetszeti regresszió, az alapváltozat OLS heteroszkedaszticitás-robusztus (HC3) standard hibákkal.
#
#   K modell (fő modell) – a kezeletlen cukorbetegek aránya (23. vizsgálat, NCD-RisC); a kód végén.
#   B modell – az aluldetektálás és a halálozás kapcsolata (14. táblázat első sora):
#       célváltozó: a cukorbetegség miatti halálesetek 1000 cukorbetegre (20–79 évesek; 4. és 2. vizsgálat hányadosa), log;
#       fő magyarázó változó: az IDF nem diagnosztizált aránya (20. vizsgálat), kontrollváltozók: fejlettség és ellátórendszer.
#   H1 / H2 – a cukorbetegek és a halálesetek népességarányos számának különbsége, parciális korreláció az IDF-aránnyal
#       (14. táblázat második és harmadik sora).

import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.stattools import jarque_bera

# BEÁLLÍTÁSOK
VIF_KUSZOB = 10.0
MIN_LEFEDETTSEG = 0.60
SZIGNIF = 0.05

# A magyarázó változók tematikus blokkjai (a vizsgalat() nevei szerint). A sorrend a hierarchikus belépés sorrendje.
VALTOZO_BLOKKOK = {
    'Gazdasági': ['gdp_per_capita', 'gdp_aranyos_egeszs_kiad', 'cuk_egeszs_kiad_egy_fore'],
    'Szociológiai-demográfiai': ['hdi_index', 'varosiasodas', 'teljes_populacio'],
    'Ellátórendszeri': ['health_coverage_index', 'orvosok_tizezer_fore', 'korhazi_agyak', 'gondozok_szama'],
    'Epidemiológiai': ['korstand_cukorbeteg_arany', 'felnottkori_elhizas'],
}

# Ezek a vizsgálatok szándékosan maradnak ki a magyarázó változók közül:
# A leíró részben feldolgozott, de a magyarázó változók közül kimaradó vizsgálatok és a kimaradás oka:
KIZART_VALTOZOK = {
    'orvosok_szama': 'népességarányosítva ugyanazt méri, mint az orvosok 10 000 lakosra jutó száma (WHO)',
    'cuk_kiadas_teljes': 'egy lakosra vetítve ugyanazt az információt hordozza, mint az egy cukorbetegre jutó kiadás',
    'cukorbetegek_szama': 'nyers, korösszetételtől függő gyakoriság; helyette a korstandardizált arány szerepel',
    'egyes_tipusu_fiatalok': '1-es típusú cukorbetegség: a nyilvántartott betegek száma a felismerési képességet tükrözi',
    'egyes_tipusu_minden_korosztaly': '1-es típusú cukorbetegség: a nyilvántartott betegek száma a felismerési képességet tükrözi',
    'felnottkori_cukorbeteg_halalozasi_arany': 'halálozási kimenet, nem ok',
    'cukorbeteg_halalesetek': 'halálozási kimenet (a B modell célváltozójának része)',
    'aluldiag_arany': 'az IDF nem diagnosztizált aránya: ellenőrzésre és a B modell fő magyarázó változójaként szerepel',
    'nyers_halalozasi_arany': 'demográfiai kontroll – csak a halálozási ellenőrzésben szerepel',
    'idosek_aranya': 'demográfiai kontroll – csak a halálozási ellenőrzésben szerepel',
    'kezeletlen_arany': 'a fő célváltozó (K modell)',
}

# Olvasható nevek a táblázatokhoz és ábrákhoz
VALTOZO_CIMKE = {
    'aluldiag_arany': 'Nem diagnosztizált arány (%)',
    'gdp_per_capita': 'Egy főre jutó GDP (PPP)',
    'gdp_aranyos_egeszs_kiad': 'Egészségügyi kiadás a GDP %-ában',
    'cuk_egeszs_kiad_egy_fore': 'Cukorbetegséggel kapcsolatos kiadás / fő',
    'hdi_index': 'Humán fejlettségi index',
    'varosiasodas': 'Városi népesség aránya',
    'teljes_populacio': 'Népesség',
    'health_coverage_index': 'UHC-index',
    'orvosok_tizezer_fore': 'Orvosok / 10 000 lakos',
    'korhazi_agyak': 'Kórházi ágyak',
    'gondozok_szama': 'Ápolók és szülésznők',
    'korstand_cukorbeteg_arany': 'Korstandardizált cukorbeteg-arány',
    'felnottkori_elhizas': 'Felnőttkori elhízás',
    'egyes_tipusu_fiatalok': '1-es típusú cukorbetegek (0–19 évesek) / 100 000 lakos',
    'nyers_halalozasi_arany': 'Nyers halálozási arány (1000 lakosra)',
    'idosek_aranya': '65 év felettiek aránya (%)',
    'kezeletlen_arany': 'Kezeletlen cukorbetegek aránya (%)',
}

# IDF-régiók (a Diabetes Atlas 7 régiója) a regionális összesítéshez. Az IDF-fájlok írásmódját használja,
# az illesztesi_kulcs() a többi forrás névváltozataira is ráilleszti. A be nem sorolt országokat a kód kiírja.
_IDF_REGIOK_LISTA = {
    'Afrika (AFR)': [
        'Angola', 'Benin', 'Botswana', 'Burkina Faso', 'Burundi', 'Cabo Verde', 'Cameroon', 'Central African Republic', 'Chad',
        'Comoros', 'Republic of Congo', 'Congo', "Côte d'Ivoire", 'Democratic Republic of the Congo', 'Djibouti', 'Equatorial Guinea',
        'Eritrea', 'Eswatini', 'Ethiopia', 'Gabon', 'Gambia', 'Ghana', 'Guinea', 'Guinea-Bissau', 'Kenya', 'Lesotho', 'Liberia',
        'Madagascar', 'Malawi', 'Mali', 'Mauritania', 'Mayotte', 'Mozambique', 'Namibia', 'Niger', 'Nigeria', 'Reunion', 'Rwanda',
        'Sao Tome and Principe', 'Senegal', 'Seychelles', 'Sierra Leone', 'Somalia', 'South Africa', 'South Sudan', 'Togo', 'Uganda',
        'United Republic of Tanzania', 'Western Sahara', 'Zambia', 'Zimbabwe'],
    'Európa (EUR)': [
        'Albania', 'Andorra', 'Armenia', 'Austria', 'Azerbaijan', 'Belarus', 'Belgium', 'Bosnia and Herzegovina', 'Bulgaria',
        'Channel Islands', 'Croatia', 'Cyprus', 'Czechia', 'Denmark', 'Estonia', 'Faroe Islands', 'Finland', 'France', 'Georgia',
        'Germany', 'Gibraltar', 'Greece', 'Greenland', 'Guernsey', 'Hungary', 'Iceland', 'Ireland', 'Isle of Man', 'Israel', 'Italy',
        'Jersey', 'Kazakhstan', 'Kosovo (under UNSC res. 1244)', 'Kosovo', 'Kyrgyzstan', 'Latvia', 'Liechtenstein', 'Lithuania',
        'Luxembourg', 'Malta', 'Monaco', 'Montenegro', 'Netherlands', 'North Macedonia', 'Norway', 'Poland', 'Portugal',
        'Republic of Moldova', 'Romania', 'Russian Federation', 'San Marino', 'Serbia', 'Slovakia', 'Slovenia', 'Spain', 'Sweden',
        'Switzerland', 'Tajikistan', 'Türkiye', 'Turkmenistan', 'Ukraine', 'United Kingdom', 'Uzbekistan'],
    'Közel-Kelet és Észak-Afrika (MENA)': [
        'Afghanistan', 'Algeria', 'Bahrain', 'Egypt', 'Iran (Islamic Republic of)', 'Iraq', 'Jordan', 'Kuwait', 'Lebanon', 'Libya',
        'Morocco', 'Oman', 'Pakistan', 'State of Palestine', 'West Bank and Gaza', 'Qatar', 'Saudi Arabia', 'Sudan',
        'Syrian Arab Republic', 'Tunisia', 'United Arab Emirates', 'Yemen'],
    'Észak-Amerika és Karib-térség (NAC)': [
        'Anguilla', 'Antigua and Barbuda', 'Aruba', 'Bahamas', 'Barbados', 'Belize', 'Bermuda', 'British Virgin Islands', 'Canada',
        'Cayman Islands', 'Curaçao', 'Dominica', 'Grenada', 'Guadeloupe', 'Guyana', 'Haiti', 'Jamaica', 'Martinique', 'Mexico',
        'Montserrat', 'Saint Kitts and Nevis', 'St. Martin (French part)', 'Saint Lucia', 'Saint Vincent and the Grenadines', 'Sint Maarten', 'Saint Martin',
        'Suriname', 'Trinidad and Tobago', 'Turks and Caicos Islands', 'United States of America', 'United States Virgin Islands'],
    'Dél- és Közép-Amerika (SACA)': [
        'Argentina', 'Bolivia (Plurinational State of)', 'Brazil', 'Chile', 'Colombia', 'Costa Rica', 'Cuba', 'Dominican Republic',
        'Ecuador', 'El Salvador', 'French Guiana', 'Guatemala', 'Honduras', 'Nicaragua', 'Panama', 'Paraguay', 'Peru',
        'Puerto Rico', 'Uruguay', 'Venezuela (Bolivarian Republic of)'],
    'Délkelet-Ázsia (SEA)': ['Bangladesh', 'Bhutan', 'India', 'Maldives', 'Mauritius', 'Nepal', 'Sri Lanka'],
    'Nyugat-Csendes-óceán (WP)': [
        'American Samoa', 'Australia', 'Brunei Darussalam', 'Cambodia', 'China', 'China, Hong Kong SAR', 'China, Macao SAR',
        'Cook Islands', "Dem. People's Republic of Korea", 'Fiji', 'French Polynesia', 'Guam', 'Indonesia', 'Japan', 'Kiribati',
        "Lao People's Democratic Republic", 'Malaysia', 'Marshall Islands', 'Micronesia (Fed. States of)', 'Mongolia', 'Myanmar',
        'Nauru', 'New Caledonia', 'New Zealand', 'Niue', 'Northern Mariana Islands', 'Palau', 'Papua New Guinea', 'Philippines',
        'Republic of Korea', 'Samoa', 'Singapore', 'Solomon Islands', 'Taiwan, Province of China', 'Thailand', 'Timor-Leste',
        'Tokelau', 'Tonga', 'Tuvalu', 'Vanuatu', 'Viet Nam', 'Wallis and Futuna'],
}
IDF_REGIO = {}
for _regio, _orszagok in _IDF_REGIOK_LISTA.items():
    for _kulcs in illesztesi_kulcs(pd.Series(_orszagok)):
        IDF_REGIO[_kulcs] = _regio


def regio_hozzarendeles(orszagnevek):
    regiok = illesztesi_kulcs(pd.Series(list(orszagnevek), index=list(orszagnevek))).map(IDF_REGIO)
    hianyzo = sorted(regiok[regiok.isna()].index)
    if hianyzo:
        print(f'FIGYELEM: IDF-régióba nem sorolt országok ({len(hianyzo)} db), a regionális összesítésből kimaradnak: {hianyzo}')
    return regiok.fillna('besorolatlan')


# SEGÉDFÜGGVÉNYEK
def modell_valtozo(nev, valtozat=None):
    """Egy vizsgálat modellbe kerülő változata országonként (alapból a javasolt: nyers / arányos / log)."""
    lepes = valtozat or javasolt_valtozat(nev)[0]
    return eredmenyek[nev].set_index('Country Name')[LEPES_OSZLOP[lepes]].rename(nev), lepes


def vif_tabla(X):
    Xc = sm.add_constant(X, has_constant='add')
    return pd.Series([variance_inflation_factor(Xc.values, i) for i in range(1, Xc.shape[1])], index=X.columns)


def vif_szures(adat, valtozok, cel, vedett=()):
    """Iteratívan elhagyja a legnagyobb VIF-ű változót, amíg minden VIF < VIF_KUSZOB. A védett változó (pl. a B modell
    fő magyarázó változója) nem esik ki. Minden lépés után újra listaszerű törlés, így a minta nőhet."""
    valtozok = list(valtozok)
    while len(valtozok) > 1:
        minta = adat[[cel] + valtozok].dropna()
        v = vif_tabla(minta[valtozok]).drop(list(vedett), errors='ignore')
        if v.empty or v.max() < VIF_KUSZOB:
            break
        ki = v.idxmax()
        print(f'   multikollinearitás: {ki} kimarad (VIF = {v.max():.1f})')
        valtozok.remove(ki)
    return valtozok


def ols(y, X, robusztus=True):
    modell = sm.OLS(y, sm.add_constant(X, has_constant='add'))
    return modell.fit(cov_type='HC3') if robusztus else modell.fit()


def egyutthato_tabla(eredmeny, y, X, lepesek):
    """Együtthatók HC3-hibákkal, 95%-os intervallummal, standardizált bétával és VIF-fel."""
    std_beta = eredmeny.params.drop('const') * X.std() / y.std()
    ci = eredmeny.conf_int()
    t = pd.DataFrame({
        'Változó': [VALTOZO_CIMKE.get(v, v) for v in X.columns],
        'Kód': X.columns,
        'Modellváltozat': [lepesek.get(v, 'nyers') for v in X.columns],
        'Béta (B)': eredmeny.params.drop('const').values,
        'Robusztus SE (HC3)': eredmeny.bse.drop('const').values,
        'p-érték': eredmeny.pvalues.drop('const').values,
        'CI 95% alsó': ci.loc[X.columns, 0].values,
        'CI 95% felső': ci.loc[X.columns, 1].values,
        'Standardizált béta': std_beta.values,
        'VIF': vif_tabla(X).values,
    })
    t['Szignifikáns (5%)'] = np.where(t['p-érték'] < SZIGNIF, 'igen', 'nem')
    return t.sort_values('Standardizált béta', key=abs, ascending=False)


def ertelmezes(sor, cel_egyseg):
    """Egy szignifikáns együttható szöveges értelmezése a változat (log / nyers) szerint."""
    b, nev = sor['Béta (B)'], sor['Változó']
    irany = 'nő' if b > 0 else 'csökken'
    if sor['Modellváltozat'] == 'log':
        return f'{nev}: 10%-os növekedése mellett a célváltozó {abs(b) * np.log(1.1):.2f} {cel_egyseg} {irany}'
    return f'{nev}: egységnyi növekedése mellett a célváltozó {abs(b):.3f} {cel_egyseg} {irany}'


def diagnosztika(eredmeny, nev):
    bp = het_breuschpagan(eredmeny.resid, eredmeny.model.exog)
    jb = jarque_bera(eredmeny.resid)
    print(f'   Diagnosztika ({nev}): Breusch–Pagan p = {bp[1]:.3f} '
          f'({"heteroszkedasztikus, ezért indokolt a HC3" if bp[1] < SZIGNIF else "nincs jele heteroszkedaszticitásnak"}); '
          f'Jarque–Bera p = {jb[1]:.3f} ({"a reziduumok nem normálisak" if jb[1] < SZIGNIF else "a reziduumok közel normálisak"}); '
          f'reziduum ferdeség = {jb[2]:.2f}')
    cook = pd.Series(eredmeny.get_influence().cooks_distance[0], index=eredmeny.model.data.row_labels)
    nagy = cook[cook > 4 / len(cook)].sort_values(ascending=False)
    print(f'   Befolyásos megfigyelések (Cook-távolság > 4/n): {len(nagy)} db; legnagyobbak: '
          + ', '.join(f'{o} ({c:.2f})' for o, c in nagy.head(5).items()))
    return {'Breusch–Pagan p': bp[1], 'Jarque–Bera p': jb[1], 'Befolyásos megfigyelések': len(nagy)}


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# KÖZÖS ORSZÁGSZINTŰ ADATTÁBLA: az IDF nem diagnosztizált aránya + minden jelölt magyarázó változó a javasolt változatban

df_aluldiag = eredmenyek['aluldiag_arany']
reg = df_aluldiag.set_index('Country Name')[['Value']].rename(columns={'Value': 'aluldiag_arany'})
valtozo_lepes = {}
for blokk, nevek in VALTOZO_BLOKKOK.items():
    for nev in nevek:
        s, lepes = modell_valtozo(nev)
        reg = reg.join(s, how='left')
        valtozo_lepes[nev] = lepes
reg['Regio'] = regio_hozzarendeles(reg.index).values

print('\nA magyarázó változók közül kimaradó vizsgálatok:')
for nev, ok in KIZART_VALTOZOK.items():
    print(f'   {nev}: {ok}')

# =========================================================================================================================================================================
# "B" MODELL – AZ ALULDETEKTÁLÁS ÉS A CUKORBETEGEK HALÁLOZÁSA
# =========================================================================================================================================================================

print('\n' + '=' * 100 + '\n"B" MODELL – az aluldetektálás és a cukorbetegek halálozása\n' + '=' * 100)

# Halálesetek (fő, 20–79 évesek) és cukorbetegek (a fájlban ezer főben, 20–79 évesek), mindkettő IDF 2024.
halal = eredmenyek['cukorbeteg_halalesetek'].set_index('Country Name')['Value'].rename('halalesetek')
betegek = (eredmenyek['cukorbetegek_szama'].set_index('Country Name')['Value'] * 1000).rename('cukorbetegek')
reg = reg.drop(columns=['halalesetek', 'cukorbetegek'], errors='ignore').join(halal, how='left').join(betegek, how='left')
reg['halal_1000_betegre'] = reg['halalesetek'] / reg['cukorbetegek'] * 1000
reg['log_halal_1000_betegre'] = np.log(reg['halal_1000_betegre'].where(reg['halal_1000_betegre'] > 0))

# Kontrollváltozók: a halálozási kockázatot az aluldetektálástól függetlenül is befolyásoló fejlettség és ellátórendszer.
# (A korstandardizált gyakoriság itt nem kell, mert a célváltozó már cukorbetegekre vetített.)
kontroll_B = [v for v in ['gdp_per_capita', 'health_coverage_index', 'orvosok_tizezer_fore', 'cuk_egeszs_kiad_egy_fore',
                          'varosiasodas', 'felnottkori_elhizas']
              if v in reg.columns and reg.loc[reg['aluldiag_arany'].notna(), v].notna().mean() >= MIN_LEFEDETTSEG]
print('Multikollinearitás-szűrés (VIF), az aluldetektálás védett változó:')
valtozok_B = vif_szures(reg, ['aluldiag_arany'] + kontroll_B, 'log_halal_1000_betegre', vedett=('aluldiag_arany',))
minta_B = reg[['log_halal_1000_betegre'] + valtozok_B].dropna()
print(f'"B" modell: {len(valtozok_B)} magyarázó változó, {len(minta_B)} ország')

# Összehasonlításképp: csak az aluldetektálással (kontroll nélkül) – megmutatja, mennyit "visz el" a fejlettség.
modell_B0 = ols(minta_B['log_halal_1000_betegre'], minta_B[['aluldiag_arany']])
modell_B = ols(minta_B['log_halal_1000_betegre'], minta_B[valtozok_B])
lepes_B = dict(valtozo_lepes, aluldiag_arany='nyers')
egyutthatok_B = egyutthato_tabla(modell_B, minta_B['log_halal_1000_betegre'], minta_B[valtozok_B], lepes_B)
print(egyutthatok_B.drop(columns='Kód').round(4).to_string(index=False))
print(f'R² = {modell_B.rsquared:.3f}, korrigált R² = {modell_B.rsquared_adj:.3f}, N = {int(modell_B.nobs)}')
diag_B = diagnosztika(modell_B, 'B modell')

beta = modell_B.params['aluldiag_arany']
b_also, b_felso = modell_B.conf_int().loc['aluldiag_arany']
p_B = modell_B.pvalues['aluldiag_arany']
print(f'\nAz aluldetektálás együtthatója: kontroll nélkül {modell_B0.params["aluldiag_arany"]:.4f} '
      f'(p = {modell_B0.pvalues["aluldiag_arany"]:.4f}), kontrollokkal {beta:.4f} (p = {p_B:.4f}, 95% CI: {b_also:.4f} – {b_felso:.4f})')
print(f'Értelmezés: 1 százalékponttal magasabb nem diagnosztizált arány mellett a cukorbetegek halálozása '
      f'{(np.exp(beta) - 1) * 100:+.2f}%-kal tér el (10 százalékpontnál {(np.exp(10 * beta) - 1) * 100:+.1f}%).')
if p_B >= SZIGNIF:
    print('FIGYELEM: az aluldetektálás hatása 5%-os szinten NEM szignifikáns.')

modell_osszefoglalo = pd.DataFrame([
    {'Modell': 'B – halálozás és aluldetektálás', 'Célváltozó': 'log(halálesetek / 1000 cukorbeteg)', 'N': int(modell_B.nobs),
     'R²': modell_B.rsquared, 'Korrigált R²': modell_B.rsquared_adj, 'F-próba p (HC3)': modell_B.f_pvalue,
     'Aluldetektálás béta': beta, 'CI 95% alsó': b_also, 'CI 95% felső': b_felso, 'p': p_B,
     'Kontroll nélküli béta': modell_B0.params['aluldiag_arany'], **diag_B},
])
with pd.ExcelWriter('regresszio_eredmenyek.xlsx') as writer:
    modell_osszefoglalo.to_excel(writer, sheet_name='Modellek', index=False)
    egyutthatok_B.to_excel(writer, sheet_name='B_egyutthatok', index=False)
print('\nMentve: regresszio_eredmenyek.xlsx')


# =========================================================================================================================================================================
# A DIABETES ATLAS ÉRTÉKEINEK BESOROLÁSA (MÉRT / PÓTOLT) ÉS A K MODELL KÖZÖS BEÁLLÍTÁSAI
# =========================================================================================================================================================================

from sklearn.linear_model import Lasso, LassoCV
from sklearn.model_selection import RepeatedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from statsmodels.regression.quantile_regression import QuantReg
from scipy import stats

VELETLEN_MAG = 42
BOOTSTRAP_N = 200
STABILITAS_KUSZOB = 0.70
CV_REDOK, CV_ISMETLES = 5, 20
SZIMULACIO_N = 2000
FEJLETTSEGI_VALTOZOK = ['gdp_per_capita', 'hdi_index', 'health_coverage_index', 'orvosok_tizezer_fore',
                        'cuk_egeszs_kiad_egy_fore', 'gdp_aranyos_egeszs_kiad']
ELMELETI_MODELL = ['gdp_per_capita', 'health_coverage_index', 'orvosok_tizezer_fore', 'korstand_cukorbeteg_arany']
VALTOZO_CIMKE['fejlettsegi_fokomponens'] = 'Fejlettségi-ellátási főkomponens'

print('\n' + '=' * 100 + '\nA DIABETES ATLAS NEM DIAGNOSZTIZÁLT ÉRTÉKEINEK BESOROLÁSA\n' + '=' * 100)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MÉRT ÉS PÓTOLT ÉRTÉKEK
# Az IDF a felmérés nélküli országokba az azonos IDF-régió és jövedelmi csoport átlagát írja, ezért ezek az értékek a
# régión belül ismétlődnek. Szabály: ha az érték az adott IDF-régióban csak egy országnál fordul elő -> saját felmérés (mért),
# pontosan két országnál -> bizonytalan, háromnál vagy többnél -> pótolt.


reg['ismetlodes_regioban'] = reg.groupby(['Regio', 'aluldiag_arany'])['aluldiag_arany'].transform('size')
reg['potlas_statusz'] = np.select([reg['ismetlodes_regioban'] == 1, reg['ismetlodes_regioban'] == 2],
                                  ['mért', 'bizonytalan'], 'pótolt')
reg.loc[reg['aluldiag_arany'].isna(), 'potlas_statusz'] = 'nincs adat'
reg['valoszinuleg_potolt'] = reg['potlas_statusz'] == 'pótolt'
_mert = reg.index[reg['potlas_statusz'] == 'mért']
print('1. A célváltozó értékeinek besorolása: ' + ', '.join(f'{k}: {v}' for k, v in reg['potlas_statusz'].value_counts().items()))
print('   Bizonytalan (páros egyezés): ' + ', '.join(
    f'{"–".join(sorted(cs.index))} ({e:.1f}%)' for (_, e), cs in reg[reg['potlas_statusz'] == 'bizonytalan'].groupby(['Regio', 'aluldiag_arany'])))
print(reg[reg['valoszinuleg_potolt']].groupby('Regio')['aluldiag_arany'].agg(['size', 'nunique'])
      .rename(columns={'size': 'pótolt országok', 'nunique': 'különböző értékek'}).to_string())
if 'Hungary' in reg.index:
    print(f'   Magyarország: {reg.loc["Hungary", "potlas_statusz"]} '
          f'– társai: {sorted(reg.index[(reg["Regio"] == reg.loc["Hungary", "Regio"]) & (reg["aluldiag_arany"] == reg.loc["Hungary", "aluldiag_arany"])].drop("Hungary"))}')


# Keresztvalidált előrejelzési hiba (RMSE) és R² ismételt k-részes keresztvalidációval. A terv(tanito, teszt) függvény adja a
# két rész magyarázóváltozó-mátrixát; a főkomponenst így minden tanítórészen újraillesztjük, a tesztrész információja nem szivárog be.
def cv_rmse(terv, adat, cel):
    hibak, r2 = [], []
    for tanito_i, teszt_i in RepeatedKFold(n_splits=CV_REDOK, n_repeats=CV_ISMETLES, random_state=VELETLEN_MAG).split(adat):
        tanito, teszt = adat.iloc[tanito_i], adat.iloc[teszt_i]
        Xtr, Xte = terv(tanito, teszt)
        b = sm.OLS(tanito[cel], sm.add_constant(Xtr, has_constant='add')).fit().params
        becsult = sm.add_constant(Xte, has_constant='add') @ b
        hibak.append(np.sqrt(np.mean((teszt[cel] - becsult) ** 2)))
        r2.append(1 - np.sum((teszt[cel] - becsult) ** 2) / np.sum((teszt[cel] - tanito[cel].mean()) ** 2))
    return np.mean(hibak), np.std(hibak) / np.sqrt(len(hibak)), np.mean(r2)


# =========================================================================================================================================================================
# MODELLTÍPUS-VÁLASZTÁS – MELYIK REGRESSZIÓ ILLIK AZ ADOTT CÉLVÁLTOZÓHOZ?
# =========================================================================================================================================================================
#
# A célváltozó típusa határozza meg a jelölteket:
#   'arany'   – 0 és 100% közötti arány (K modell):        lineáris (OLS) | logit-transzformált OLS | frakcionális logit (GLM)
#   'pozitiv' – pozitív, jobbra ferde ráta (B modell):     lineáris (OLS) | log-OLS (Duan-korrekcióval) | Gamma-GLM log linkkel
#
# Lépések:
#   1. RESET-teszt a lineáris modellen: ha p < 0,05, a lineáris függvényforma nem megfelelő.
#   2. Minden jelölt keresztvalidált hibája (RMSE) az eredeti skálán, így a jelöltek összevethetők.
#   3. Választás: a legkisebb hibájú jelölthöz képest legfeljebb egy standard hibával rosszabbak közül az adattípushoz
#      elméletileg legjobban illőt választjuk.
#   4. Viszonyítási alap: véletlen erdő, amely magától megtalálja a nemlinearitásokat és kölcsönhatásokat. Ha a választott
#      modellnél több mint 10%-kal pontosabb, kimaradt nemlinearitás valószínű. A véletlen erdő nem választható.

from statsmodels.stats.diagnostic import linear_reset
from sklearn.ensemble import RandomForestRegressor

MODELLTIPUS_CV = dict(n_splits=5, n_repeats=10, random_state=VELETLEN_MAG)
RF_ELONY_KUSZOB = 0.10
ELMELETI_SORREND = {
    'arany': ['Frakcionális logit (GLM)', 'Logit-transzformált OLS', 'Lineáris (OLS)'],
    'pozitiv': ['Gamma-GLM (log link)', 'Log-OLS (Duan-korrekcióval)', 'Lineáris (OLS)'],
}
modelltipus_osszevetesek, modelltipus_valasztasok = [], []


def _konst(X):
    return sm.add_constant(X, has_constant='add')


def jelolt_illesztok(tipus):
    """Jelöltenként egy (illesztő, becslő) pár. Az illesztő: (X, y, robusztus) -> modell; a becslő: (modell, X) -> becslés
    az eredeti skálán. robusztus=True esetén heteroszkedaszticitás-robusztus standard hibák (OLS: HC3, GLM: HC0)."""
    ols_cov = lambda r: dict(cov_type='HC3') if r else {}
    glm_cov = lambda r: dict(cov_type='HC0') if r else {}
    if tipus == 'arany':
        return {
            'Lineáris (OLS)': (lambda X, y, r=False: sm.OLS(y, _konst(X)).fit(**ols_cov(r)),
                               lambda m, X: m.predict(_konst(X))),
            'Logit-transzformált OLS': (lambda X, y, r=False: sm.OLS(np.log(y / (100 - y)), _konst(X)).fit(**ols_cov(r)),
                                        lambda m, X: 100 / (1 + np.exp(-m.predict(_konst(X))))),
            'Frakcionális logit (GLM)': (lambda X, y, r=False: sm.GLM(y / 100, _konst(X), family=sm.families.Binomial()).fit(**glm_cov(r)),
                                         lambda m, X: 100 * m.predict(_konst(X))),
        }
    def log_ols(X, y, r=False):
        m = sm.OLS(np.log(y), _konst(X)).fit(**ols_cov(r))
        m.duan = np.mean(np.exp(m.resid))   # Duan-féle korrekció: a log-skáláról visszaszámolt átlag torzításának javítása
        return m
    return {
        'Lineáris (OLS)': (lambda X, y, r=False: sm.OLS(y, _konst(X)).fit(**ols_cov(r)),
                           lambda m, X: m.predict(_konst(X))),
        'Log-OLS (Duan-korrekcióval)': (log_ols, lambda m, X: np.exp(m.predict(_konst(X))) * m.duan),
        'Gamma-GLM (log link)': (lambda X, y, r=False: sm.GLM(y, _konst(X), family=sm.families.Gamma(sm.families.links.Log())).fit(**glm_cov(r)),
                                 lambda m, X: m.predict(_konst(X))),
    }


def modelltipus_valasztas(kod, cim, tipus, X, y, kiemelt=None):
    """Egy modell regressziótípusának kiválasztása. X: magyarázó változók, y: célváltozó az eredeti skálán.
    kiemelt: a változó, amelynek hatását külön kiírjuk (pl. aluldetektálás)."""
    print('\n' + '-' * 100 + f'\nMODELLTÍPUS-VÁLASZTÁS – {kod}: {cim}  (célváltozó típusa: {tipus}, N = {len(y)})\n' + '-' * 100)
    X = X.astype(float)
    y = y.astype(float)
    jeloltek = jelolt_illesztok(tipus)
    megjegyzesek = []

    # 1. RESET-teszt (függvényforma) a lineáris / log-lineáris alapmodellen
    alap_y = y if tipus == 'arany' else np.log(y)
    reset_p = linear_reset(sm.OLS(alap_y, _konst(X)).fit(), power=2, use_f=True).pvalue
    megjegyzesek.append(f'RESET-teszt p = {reset_p:.3f} ({"a lineáris függvényforma NEM megfelelő" if reset_p < SZIGNIF else "a függvényforma megfelelő"})')
    sorrend = list(ELMELETI_SORREND[tipus])

    # 2. Keresztvalidált hiba az eredeti skálán + véletlen erdő viszonyítási alapként
    hibak = {n: [] for n in list(jeloltek) + ['Véletlen erdő (csak viszonyítás)']}
    for tr, te in RepeatedKFold(**MODELLTIPUS_CV).split(X):
        Xtr, Xte = X.iloc[tr], X.iloc[te]
        warnings.filterwarnings('ignore')
        warnings.filterwarnings('default', category=DeprecationWarning)
        for n, (illeszt, becsul) in jeloltek.items():
            try:
                becsles = becsul(illeszt(Xtr, y.iloc[tr]), Xte)
                hibak[n].append(np.sqrt(np.mean((y.iloc[te].values - np.asarray(becsles)) ** 2)))
            except Exception:
                pass
        rf = RandomForestRegressor(n_estimators=300, min_samples_leaf=3, random_state=VELETLEN_MAG).fit(Xtr, y.iloc[tr])
        hibak['Véletlen erdő (csak viszonyítás)'].append(np.sqrt(np.mean((y.iloc[te].values - rf.predict(Xte)) ** 2)))
    t = pd.DataFrame({'Modelltípus': list(hibak),
                      'CV-RMSE (eredeti skála)': [np.mean(h) if h else np.nan for h in hibak.values()],
                      'SE': [np.std(h) / np.sqrt(len(h)) if h else np.nan for h in hibak.values()],
                      'Sikeres CV-futások': [len(h) for h in hibak.values()]})

    # 3. Választás: a legjobb + 1 SE sávon belül az elméletileg legjobban illő
    param = t[t['Modelltípus'].isin(jeloltek)].dropna(subset=['CV-RMSE (eredeti skála)'])
    legjobb = param.loc[param['CV-RMSE (eredeti skála)'].idxmin()]
    sav = param[param['CV-RMSE (eredeti skála)'] <= legjobb['CV-RMSE (eredeti skála)'] + legjobb['SE']]['Modelltípus'].tolist()
    valasztott = next(n for n in sorrend if n in sav)
    t['Választott'] = np.where(t['Modelltípus'] == valasztott, 'igen', '')
    t.insert(0, 'Modell', kod)

    # 4. Véletlen erdő összevetés
    rf_hiba = t.loc[t['Modelltípus'].str.startswith('Véletlen'), 'CV-RMSE (eredeti skála)'].iloc[0]
    v_hiba = t.loc[t['Modelltípus'] == valasztott, 'CV-RMSE (eredeti skála)'].iloc[0]
    rf_elony = (v_hiba - rf_hiba) / v_hiba
    megjegyzesek.append(f'véletlen erdő {rf_elony:+.0%}-kal {"pontosabb" if rf_elony > 0 else "pontatlanabb"} '
                        f'({"jelentős – kimaradt nemlinearitás / kölcsönhatás valószínű" if rf_elony > RF_ELONY_KUSZOB else "nem jelentős – a választott modell elegendő"})')
    print(t.drop(columns='Modell').round(3).to_string(index=False))
    for m_ in megjegyzesek:
        print('   ' + m_)
    indok = (f'{valasztott}: ' + ('a legkisebb hibájú' if valasztott == legjobb['Modelltípus'] else
                                  f'a legkisebb hibájúhoz ({legjobb["Modelltípus"]}) képest egy standard hibán belül, és az adattípushoz elméletileg jobban illik'))
    print(f'   VÁLASZTOTT MODELLTÍPUS: {indok}')

    # A választott modell a teljes mintán, robusztus standard hibákkal, értelmezhető hatásokkal
    illeszt, _ = jeloltek[valasztott]
    m = illeszt(X, y, True)
    k = X.shape[1] + 1
    params = pd.Series(np.asarray(m.params)[:k], index=['const'] + list(X.columns))
    pvals = pd.Series(np.asarray(m.pvalues)[:k], index=params.index)
    ci = pd.DataFrame(np.asarray(m.conf_int())[:k], index=params.index)

    # Értelmezés: log-változatú magyarázó változónál 10%-os növekedésre, egyébként egységnyi növekedésre (a főkomponens értéke
    # nincs újra standardizálva, ezért ott is egységnyi); arány-célváltozónál százalékpontban, log linkes modellnél %-ban.

    def lepeskoz(v):
        if v == 'fejlettsegi_fokomponens':
            return 1.0, 'egységnyi növekedésére'
        lepes = valtozo_lepes.get(v) or (javasolt_valtozat(v)[0] if v in eredmenyek else 'nyers')
        return (np.log(1.1), '10%-os növekedésére') if lepes == 'log' else (1.0, 'egységnyi növekedésére')
    p_ = 1 / (1 + np.exp(-(_konst(X).values @ params.values))) if (tipus == 'arany' and valasztott != 'Lineáris (OLS)') else None
    hatas = {}
    for v in X.columns:
        d, szoveg = lepeskoz(v)
        if p_ is not None:
            hatas[v] = f'{np.mean(p_ * (1 - p_)) * params[v] * d * 100:+.3f} %pont a változó {szoveg} (átlagos marginális hatás)'
        elif tipus == 'arany':
            hatas[v] = f'{params[v] * d:+.3f} %pont a változó {szoveg}'
        elif valasztott == 'Lineáris (OLS)':
            hatas[v] = f'{params[v] * d:+.3f} egység a változó {szoveg}'
        else:
            hatas[v] = f'{(np.exp(params[v] * d) - 1) * 100:+.2f}% a változó {szoveg}'
    egy = pd.DataFrame({'Modell': kod, 'Modelltípus': valasztott, 'Változó': [VALTOZO_CIMKE.get(v, v) for v in X.columns],
                        'Együttható': params[X.columns].values, 'p-érték': pvals[X.columns].values,
                        'Értelmezett hatás': [hatas[v] for v in X.columns]})
    print(egy.drop(columns=['Modell', 'Modelltípus']).round(4).to_string(index=False))
    kiemelt_szoveg = ''
    if kiemelt is not None and kiemelt in X.columns and tipus != 'arany' and valasztott != 'Lineáris (OLS)':
        b, (lo, hi) = params[kiemelt], ci.loc[kiemelt]
        kiemelt_szoveg = (f'{VALTOZO_CIMKE.get(kiemelt, kiemelt)} +10 egység: {(np.exp(10 * b) - 1) * 100:+.1f}% '
                          f'(95% CI: {(np.exp(10 * lo) - 1) * 100:+.1f}% – {(np.exp(10 * hi) - 1) * 100:+.1f}%), p = {pvals[kiemelt]:.3f}')
        print('   A választott modelltípussal: ' + kiemelt_szoveg)
    modelltipus_osszevetesek.append(t)
    modelltipus_valasztasok.append({'Modell': kod, 'Célváltozó': cim, 'Típus': tipus, 'N': len(y), 'Választott modelltípus': valasztott,
                                    'Indoklás': indok, 'Tesztek': '; '.join(megjegyzesek), 'Kiemelt hatás': kiemelt_szoveg,
                                    'Együtthatók': egy})
    return valasztott, m


# Alkalmazás a B modellre (a K modellre a K szakaszban fut)
modelltipus_valasztas('B', 'halálesetek 1000 cukorbetegre', 'pozitiv',
                      minta_B[valtozok_B], np.exp(minta_B['log_halal_1000_betegre']), kiemelt='aluldiag_arany')


# =========================================================================================================================================================================
# HALÁLOZÁSI ELLENŐRZÉS – A CUKORBETEGEK ÉS A CUKORBETEGSÉG MIATTI HALÁLESETEK NÉPESSÉGARÁNYOS SZÁMÁNAK KÜLÖNBSÉGE
# =========================================================================================================================================================================
#
#   H1 – nyers különbség:         cukorbetegek / 100 000 lakos − halálesetek / 100 000 lakos
#   H2 – logaritmikus különbség:  ln(cukorbetegek / 100 000 lakos) − ln(halálesetek / 100 000 lakos) = ln(cukorbetegek / halálesetek)
#
# Mindkettőt parciális korrelációval vetjük össze az IDF nem diagnosztizált arányával: a célváltozóból és az IDF-arányból is
# kiszűrjük a demográfiai kontrollt (65 év felettiek aránya, nyers halálozási arány), és a maradékok korrelációját nézzük,
# minden országon és csak a saját felméréssel rendelkező (mért) országokon. Elvárt előjel: negatív (több felismeretlen beteg
# mellett egy cukorbetegre több haláleset jut).

print('\n' + '=' * 100 + '\nHALÁLOZÁSI ELLENŐRZÉS – a cukorbetegek és a halálesetek népességarányos számának különbsége\n' + '=' * 100)

FO_EGYSEG = 100_000
DEMOGRAFIAI_KONTROLL = ['idosek_aranya', 'nyers_halalozasi_arany']
VALTOZO_CIMKE.update({'H1': 'Nyers különbség (cukorbetegek − halálesetek, 100 000 lakosra)',
                      'H2': 'Logaritmikus különbség (ln cukorbetegek − ln halálesetek)'})

hreg = pd.DataFrame({
    # 2. vizsgálat: az arányos érték a népesség %-a -> 100 000 lakosra
    'cukorbetegek_100e': eredmenyek['cukorbetegek_szama'].set_index('Country Name')['Aranyos'] * FO_EGYSEG / 100,
    # 4. vizsgálat: az arányos érték már 100 000 lakosra jutó haláleset
    'halalesetek_100e': eredmenyek['cukorbeteg_halalesetek'].set_index('Country Name')['Aranyos'],
}).dropna()
_nem_pozitiv = hreg.index[(hreg <= 0).any(axis=1)]
if len(_nem_pozitiv):
    print(f'Kimarad (nulla haláleset vagy beteg, a log nem értelmezhető): {sorted(_nem_pozitiv)}')
hreg = hreg.drop(_nem_pozitiv)
hreg['H1'] = hreg['cukorbetegek_100e'] - hreg['halalesetek_100e']
hreg['H2'] = np.log(hreg['cukorbetegek_100e']) - np.log(hreg['halalesetek_100e'])
for _v in DEMOGRAFIAI_KONTROLL:
    hreg = hreg.join(modell_valtozo(_v)[0], how='left')
hreg = hreg.join(reg[['aluldiag_arany', 'potlas_statusz']], how='left')
hreg['potlas_statusz'] = hreg['potlas_statusz'].fillna('nincs adat')

_r_prev = hreg[['H1', 'cukorbetegek_100e']].corr().iloc[0, 1]
print(f'A H1 korrelációja a népességarányos cukorbetegszámmal: r = {_r_prev:.4f}; a halálesetek a cukorbetegek számának '
      f'átlagosan {(hreg["halalesetek_100e"] / hreg["cukorbetegek_100e"]).mean() * 100:.2f}%-át teszik ki.')


def parcialis_korrelacio(adat, a, b, kontroll):
    d = adat[[a, b] + kontroll].dropna()
    if len(d) < len(kontroll) + 10:
        return dict(r=np.nan, p=np.nan, rho=np.nan, N=len(d))
    if kontroll:
        Xk = sm.add_constant(d[kontroll], has_constant='add')
        ra, rb = sm.OLS(d[a], Xk).fit().resid, sm.OLS(d[b], Xk).fit().resid
    else:
        ra, rb = d[a] - d[a].mean(), d[b] - d[b].mean()
    r = np.corrcoef(ra, rb)[0, 1]
    df_ = len(d) - 2 - len(kontroll)
    p = 2 * stats.t.sf(abs(r * np.sqrt(df_ / (1 - r ** 2))), df_)
    return dict(r=r, p=p, rho=stats.spearmanr(ra, rb)[0], N=len(d))


valid_sorok = []
for kod in ['H1', 'H2']:
    for minta_nev, maszk in [('minden ország', hreg['aluldiag_arany'].notna()), ('mért országok', hreg['potlas_statusz'] == 'mért')]:
        for kontroll_nev, kontroll in [('kontroll nélkül', []), ('demográfiai kontrollal', DEMOGRAFIAI_KONTROLL)]:
            k = parcialis_korrelacio(hreg[maszk], kod, 'aluldiag_arany', kontroll)
            valid_sorok.append({'Célváltozó': kod, 'Minta': minta_nev, 'Kontroll': kontroll_nev, 'N': k['N'],
                                'Pearson r': k['r'], 'p': k['p'], 'Spearman ρ': k['rho'],
                                'Várt irány (negatív)': 'igen' if k['r'] < 0 else 'nem'})
validacio_H = pd.DataFrame(valid_sorok)
print('\nA célváltozók kapcsolata az IDF nem diagnosztizált arányával (parciális korreláció):')
print(validacio_H.round(3).to_string(index=False))

with pd.ExcelWriter('regresszio_halalozas_kulonbseg.xlsx') as writer:
    hreg.reset_index(names='Country Name').to_excel(writer, sheet_name='Celvaltozok_orszagonkent', index=False)
    validacio_H.to_excel(writer, sheet_name='Validacio_IDF', index=False)
print('\nMentve: regresszio_halalozas_kulonbseg.xlsx')

# =========================================================================================================================================================================
# FŐ MODELL (K) – A KEZELETLEN CUKORBETEGEK ARÁNYA: OKOK, ELLENŐRZÉS AZ IDF-ADATOKKAL, MEGOLDÁSI FORGATÓKÖNYVEK
# =========================================================================================================================================================================
#
# A kezeletlen arány (23. vizsgálat, NCD-RisC) független forrásból, vérvizsgálatos felmérésekre épülve, kb. 200 országra méri,
# hogy a cukorbetegek hány százaléka nem kap kezelést. Aki nincs diagnosztizálva, az kezelést sem kap, ezért ez a felismerés
# és az ellátás együttes hiányát méri. Az IDF nem diagnosztizált aránya kb. 85 országban régiós átlag, ezért itt csak
# ellenőrzésre szolgál: ha a két független forrás ugyanazokat az okokat mutatja, az mindkettőt megerősíti.
#
# Lépések:
#   1. Adattábla: célváltozó + a magyarázó változók javasolt változata + IDF-arány és pótlási státusz + bizonytalanság
#   2. Okcsoportok magyarázó ereje (hierarchikus blokkok, VIF-szűrt változókkal)
#   3. Változószelekció: LASSO-stabilitás (bootstrap) és 5 versengő modell keresztvalidációval, "legjobb + 1 SE" szabály
#   4. A választott modell együtthatói és diagnosztikája
#   5. Robusztusság 7 specifikációban; 5b. adatév-érzékenység; 5c. a kihagyott változók ellenőrzése
#   6. Modelltípus-választás (arány-célváltozó: OLS / logit-OLS / frakcionális logit)
#   7. Ellenőrzés: korreláció az IDF nem diagnosztizált arányával, és ugyanaz a modell az IDF-arányra a mért országokon
#   8. Országok: a vártnál rosszabb / jobb helyzetűek (reziduumok), Magyarország
#   9. Megoldási forgatókönyvek: ha az ellátórendszeri tényezők elérnék a saját IDF-régió felső negyedét, mennyivel
#      csökkenne a kezeletlen arány és hány cukorbeteg kerülne kezelésbe (95%-os szimulációs intervallummal)
#  10. Mentés: regresszio_kezeletlen_arany.xlsx, regresszio_K_*.png

K_CEL = 'kezeletlen_arany'
FORGATOKONYV_BEAVATKOZAS = ['health_coverage_index', 'orvosok_tizezer_fore']   # az egészségpolitika által alakítható tényezők
FORGATOKONYV_KVANTILIS = 0.75     # a cél: a saját IDF-régió felső negyedének határa
VALTOZO_CIMKE[K_CEL] = 'Kezeletlen cukorbetegek aránya (%)'

print('\n' + '=' * 100 + '\nFŐ MODELL (K) – a kezeletlen cukorbetegek aránya (NCD-RisC)\n' + '=' * 100)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. ADATTÁBLA
kreg = eredmenyek[K_CEL].set_index('Country Name')[['Value']].rename(columns={'Value': K_CEL})
for _v in valtozo_lepes:
    kreg = kreg.join(modell_valtozo(_v)[0], how='left')
kreg = kreg.join(reg[['aluldiag_arany', 'potlas_statusz']], how='left')
kreg['potlas_statusz'] = kreg['potlas_statusz'].fillna('nincs adat')
kreg['Regio'] = regio_hozzarendeles(kreg.index).values
_ci = kezeles_ci.set_index('Match_Key')['CI_szelesseg']
kreg['CI_szelesseg'] = illesztesi_kulcs(pd.Series(kreg.index, index=kreg.index)).map(_ci).values
kreg['cukorbetegek_fo'] = eredmenyek['cukorbetegek_szama'].set_index('Country Name')['Value'].reindex(kreg.index) * 1000
print(f'1. {len(kreg)} ország; a kezeletlen arány átlaga {kreg[K_CEL].mean():.1f}%, mediánja {kreg[K_CEL].median():.1f}%')

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. OKCSOPORTOK MAGYARÁZÓ EREJE (mint az A modellben)
van_K = kreg[K_CEL].notna()
jelolt_K = [v for v in valtozo_lepes if kreg.loc[van_K, v].notna().mean() >= MIN_LEFEDETTSEG]
print('2. Multikollinearitás-szűrés (VIF):')
fo_K = vif_szures(kreg, jelolt_K, K_CEL)
minta_K = kreg[[K_CEL] + fo_K].dropna()
blokk_K, elozo, elozo_nr, bent = [], None, None, []
for blokk, nevek in VALTOZO_BLOKKOK.items():
    uj = [v for v in nevek if v in fo_K]
    if not uj:
        continue
    bent = bent + uj
    m = ols(minta_K[K_CEL], minta_K[bent])
    m_nr = ols(minta_K[K_CEL], minta_K[bent], robusztus=False)
    sor = {'Lépés': f'+ {blokk}', 'Belépő változók': ', '.join(VALTOZO_CIMKE.get(v, v) for v in uj), 'Változók száma': len(bent),
           'N': int(m.nobs), 'R²': m.rsquared, 'Korrigált R²': m.rsquared_adj,
           'ΔR²': m.rsquared - (elozo.rsquared if elozo is not None else 0)}
    if elozo_nr is not None:
        sor['F (ΔR²)'], sor['p (ΔR²)'], _ = m_nr.compare_f_test(elozo_nr)
    else:
        sor['F (ΔR²)'], sor['p (ΔR²)'] = m_nr.fvalue, m_nr.f_pvalue
    blokk_K.append(sor)
    elozo, elozo_nr = m, m_nr
blokk_K = pd.DataFrame(blokk_K)
print(blokk_K.round(3).to_string(index=False))

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. VÁLTOZÓSZELEKCIÓ: LASSO-stabilitás (bootstrap) és versengő modellek keresztvalidációval
szel_K = [v for v in valtozo_lepes if kreg.loc[van_K, v].notna().mean() >= MIN_LEFEDETTSEG]
minta_Ksz = kreg[[K_CEL] + szel_K].dropna()
X_ = StandardScaler().fit_transform(minta_Ksz[szel_K])
y_ = minta_Ksz[K_CEL].values
lasso_K = LassoCV(cv=CV_REDOK, random_state=VELETLEN_MAG, max_iter=50_000).fit(X_, y_)
rng_K = np.random.default_rng(VELETLEN_MAG)
kiv = np.zeros(len(szel_K))
for _ in range(BOOTSTRAP_N):
    idx = rng_K.integers(0, len(y_), len(y_))
    kiv += Lasso(alpha=lasso_K.alpha_, max_iter=50_000).fit(X_[idx], y_[idx]).coef_ != 0
gyak_K = pd.Series(kiv / BOOTSTRAP_N, index=szel_K)
stabil_K = gyak_K[gyak_K >= STABILITAS_KUSZOB].sort_values(ascending=False).index.tolist() or \
    gyak_K.sort_values(ascending=False).head(3).index.tolist()
szures_K = pd.DataFrame({
    'Változó': [VALTOZO_CIMKE.get(v, v) for v in szel_K],
    'Korreláció a célváltozóval (r)': [minta_Ksz[[K_CEL, v]].corr().iloc[0, 1] for v in szel_K],
    'VIF': vif_tabla(minta_Ksz[szel_K]).values,
    'LASSO-együttható (std.)': lasso_K.coef_,
    'LASSO kiválasztási gyakoriság': gyak_K.values,
    'A VIF-szűrt modellben': ['igen' if v in fo_K else 'nem' for v in szel_K]}, index=szel_K)
print(f'\n3. Változószelekció (N = {len(minta_Ksz)}); LASSO-stabil (≥ {STABILITAS_KUSZOB:.0%}): {[VALTOZO_CIMKE.get(v, v) for v in stabil_K]}')
print(szures_K.round(2).to_string(index=False))

fejl_K = [v for v in FEJLETTSEGI_VALTOZOK if v in szel_K]
egyeb_K = [v for v in szel_K if v not in fejl_K]

def pca_terv_K(tanito, teszt, tovabbi):
    sk = StandardScaler().fit(tanito[fejl_K])
    pca = PCA(n_components=1).fit(sk.transform(tanito[fejl_K]))
    elojel = np.sign(pca.components_[0][fejl_K.index('gdp_per_capita')] if 'gdp_per_capita' in fejl_K else 1)
    def tr(d):
        ki = d[tovabbi].copy()
        ki['fejlettsegi_fokomponens'] = elojel * pca.transform(sk.transform(d[fejl_K]))[:, 0]
        return ki
    return tr(tanito), tr(teszt)

_valaszt = lambda lista: (lambda tr, te: (tr[[v for v in lista if v in tr]], te[[v for v in lista if v in te]]))
versengo_K = {
    'M1 – VIF-szűrt teljes modell': _valaszt(fo_K),
    'M2 – LASSO-stabil változók': _valaszt(stabil_K),
    'M3 – fejlettségi főkomponens + egyéb változók': lambda tr, te: pca_terv_K(tr, te, egyeb_K),
    'M4 – főkomponens + LASSO-stabil nem fejlettségi változók': lambda tr, te: pca_terv_K(tr, te, [v for v in stabil_K if v not in fejl_K]),
    'M5 – elméleti szűk modell': _valaszt(ELMELETI_MODELL),
}
sorok_K, illesztett_K = [], {}
for nev, terv in versengo_K.items():
    Xt, _ = terv(minta_Ksz, minta_Ksz)
    m = ols(minta_Ksz[K_CEL], Xt)
    rmse, se, cvr2 = cv_rmse(terv, minta_Ksz, K_CEL)
    illesztett_K[nev] = (m, Xt)
    sorok_K.append({'Modell': nev, 'Változók': ', '.join(VALTOZO_CIMKE.get(c, c) for c in Xt.columns), 'k': Xt.shape[1],
                    'N': int(m.nobs), 'R²': m.rsquared, 'Korrigált R²': m.rsquared_adj, 'AIC': m.aic, 'BIC': m.bic,
                    'CV-RMSE (%pont)': rmse, 'CV-RMSE SE': se, 'CV-R²': cvr2,
                    'Max. VIF': vif_tabla(Xt).max() if Xt.shape[1] > 1 else 1.0})
versengo_tabla_K = pd.DataFrame(sorok_K)
_legjobb = versengo_tabla_K.loc[versengo_tabla_K['CV-RMSE (%pont)'].idxmin()]
_hatar = _legjobb['CV-RMSE (%pont)'] + _legjobb['CV-RMSE SE']
valasztott_K = versengo_tabla_K[versengo_tabla_K['CV-RMSE (%pont)'] <= _hatar].sort_values(['k', 'CV-RMSE (%pont)']).iloc[0]['Modell']
versengo_tabla_K['Választott'] = np.where(versengo_tabla_K['Modell'] == valasztott_K, 'igen', '')
print(f'\n   Versengő modellek (CV: {CV_REDOK}-részes, {CV_ISMETLES}× ismételve):')
print(versengo_tabla_K.drop(columns='Változók').round(3).to_string(index=False))
print(f'   Választott ("legjobb + 1 SE"): {valasztott_K}')

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4. A VÁLASZTOTT MODELL
modell_K, X_K = illesztett_K[valasztott_K]
lepes_K = dict(valtozo_lepes, fejlettsegi_fokomponens='nyers')
egyutthatok_K = egyutthato_tabla(modell_K, minta_Ksz[K_CEL], X_K, lepes_K)
print(f'\n4. A választott modell (célváltozó: kezeletlen arány, %): R² = {modell_K.rsquared:.3f}, '
      f'korrigált R² = {modell_K.rsquared_adj:.3f}, N = {int(modell_K.nobs)}')
if 'fejlettsegi_fokomponens' in X_K.columns:
    _sk = StandardScaler().fit(minta_Ksz[fejl_K])
    _pca = PCA(n_components=1).fit(_sk.transform(minta_Ksz[fejl_K]))
    _tolt = pd.Series(_pca.components_[0], index=[VALTOZO_CIMKE.get(v, v) for v in fejl_K])
    _tolt *= np.sign(_tolt.get(VALTOZO_CIMKE['gdp_per_capita'], 1))
    print(f'   A fejlettségi főkomponens a {len(fejl_K)} változó szórásának {_pca.explained_variance_ratio_[0]:.0%}-át sűríti; súlyok: '
          + ', '.join(f'{k} {v:+.2f}' for k, v in _tolt.items()))
print(egyutthatok_K.drop(columns='Kód').round(4).to_string(index=False))
diag_K = diagnosztika(modell_K, 'K modell')
print('   Szignifikáns tényezők értelmezése (ceteris paribus):')
for _, sor in egyutthatok_K[egyutthatok_K['p-érték'] < SZIGNIF].iterrows():
    if sor['Kód'] == 'fejlettsegi_fokomponens':
        print(f'      {sor["Változó"]}: egységnyi növekedése mellett a kezeletlen arány {abs(sor["Béta (B)"]):.2f} százalékponttal '
              f'{"nő" if sor["Béta (B)"] > 0 else "csökken"}')
    else:
        print('      ' + ertelmezes(sor, 'százalékponttal'))

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 5. ROBUSZTUSSÁG – 7 specifikáció; robusztus az a tényező, amely mindenhol azonos előjelű és legalább 5-ben szignifikáns
Xv_K = list(X_K.columns)
rob_K = kreg.loc[minta_Ksz.index].copy()
rob_K[Xv_K] = X_K
rob_K['logit_cel'] = np.log(rob_K[K_CEL] / (100 - rob_K[K_CEL]))
rob_K['precizios_suly'] = 1 / rob_K['CI_szelesseg'].clip(lower=1) ** 2   # szűkebb intervallum = pontosabb becslés = nagyobb súly
cook_K = pd.Series(modell_K.get_influence().cooks_distance[0], index=minta_Ksz.index)

def rob_K_becsles(nev, adat, cel=K_CEL, regio_fix=False, wls=False, median=False):
    X = adat[Xv_K].copy()
    if regio_fix:
        X = X.join(pd.get_dummies(adat['Regio'], prefix='regio', drop_first=True, dtype=float))
    X = sm.add_constant(X, has_constant='add')
    if median:
        m = QuantReg(adat[cel], X).fit(q=0.5)
    elif wls:
        ok = adat['precizios_suly'].notna()
        m = sm.WLS(adat.loc[ok, cel], X[ok], weights=adat.loc[ok, 'precizios_suly']).fit(cov_type='HC3')
    else:
        m = sm.OLS(adat[cel], X).fit(cov_type='HC3')
    return nev, int(m.nobs), m.params[Xv_K], m.pvalues[Xv_K]

spec_K = [
    rob_K_becsles('Alap (OLS, HC3)', rob_K),
    rob_K_becsles('IDF-régió dummy változókkal', rob_K, regio_fix=True),
    rob_K_becsles('Logit-transzformált célváltozó', rob_K, cel='logit_cel'),
    rob_K_becsles('A becslés pontosságával súlyozva (WLS)', rob_K, wls=True),
    rob_K_becsles('Adatgazdag országok (IDF saját felmérés)', rob_K[rob_K['potlas_statusz'] == 'mért']),
    rob_K_becsles('Befolyásos országok nélkül (Cook > 4/n)', rob_K[cook_K <= 4 / len(cook_K)]),
    rob_K_becsles('Medián-regresszió', rob_K, median=True),
]
rob_tabla_K = pd.DataFrame({f'{n} (N = {nn})': [f'{b:+.3f} ({p:.3f}){"*" if p < SZIGNIF else ""}' for b, p in zip(bb, pp)]
                            for n, nn, bb, pp in spec_K}, index=[VALTOZO_CIMKE.get(v, v) for v in Xv_K])
_elojel = pd.DataFrame({n: np.sign(bb.values) for n, _, bb, _ in spec_K}, index=Xv_K)
_szign = pd.DataFrame({n: (pp.values < SZIGNIF) for n, _, _, pp in spec_K}, index=Xv_K).sum(axis=1)
rob_tabla_K['Azonos előjel mindenhol'] = np.where(_elojel.nunique(axis=1).values == 1, 'igen', 'nem')
rob_tabla_K[f'Szignifikáns (db / {len(spec_K)})'] = _szign.values
rob_tabla_K['Robusztus tényező'] = np.where((rob_tabla_K['Azonos előjel mindenhol'] == 'igen') & (_szign.values >= 5), 'IGEN', 'nem')
print('\n5. Robusztusság (béta (p), * = 5%-on szignifikáns; a logit-oszlop skálája eltér, ott csak az előjel és a p számít):')
print(rob_tabla_K.T.to_string())
print(f'   Robusztus tényezők: {rob_tabla_K.index[rob_tabla_K["Robusztus tényező"] == "IGEN"].tolist() or "nincs"}')

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 5b. ADATÉV-ÉRZÉKENYSÉG: a választott modell újrabecslése úgy, hogy a Világbank- és WHO-változók legfeljebb 2022-esek
# (a célváltozó éve). Csak az évválasztás változik; a változók köre, a transzformáció (nyers / arányos / log), a
# főkomponens felépítése és a becslési eljárás ugyanaz. Az IDF-változók (csak 2024) és a HDI (csak 2023) nem cserélhetők,
# ezek változatlanok maradnak – ezt a kiírt táblázat jelzi.
ERZ_ALSO_EV, ERZ_FELSO_EV = 2019, 2022     # 2019-től, hogy a csak páratlan években közölt UHC-indexnek is legyen értéke (2021)
print(f'\n5b. Adatév-érzékenység: a magyarázó változók {ERZ_ALSO_EV}–{ERZ_FELSO_EV} közötti legfrissebb értékével')

_alapvaltozok = list(dict.fromkeys((fejl_K if 'fejlettsegi_fokomponens' in Xv_K else [])
                                   + [v for v in Xv_K if v != 'fejlettsegi_fokomponens']))
kreg_erz = kreg.copy()
ev_sorok_erz = []
for _v in _alapvaltozok:
    _fo_ev = int(eredmenyek[_v]['Ev'].iloc[0])
    _sor = {'Változó': VALTOZO_CIMKE.get(_v, _v), 'Fő modell éve': _fo_ev}
    if _fo_ev <= ERZ_FELSO_EV:
        _sor.update({'Ellenőrzés éve': _fo_ev, 'Megjegyzés': 'már a fő modellben is legfeljebb 2022-es – nem változik'})
    elif _v not in bemenetek:
        _sor.update({'Ellenőrzés éve': _fo_ev, 'Megjegyzés': 'nincs eltárolt bemenet – nem változik'})
    else:
        _b = bemenetek[_v]
        try:
            _h = _b['hosszu'].copy()
            _h.attrs['ev_forras'] = _b['ev_forras']
            _d = vizsgalat(_v, _h, _b['aranyositas'], _b['egyseg_szorzo'], _b['per_fo'],
                           also_ev=ERZ_ALSO_EV, felso_ev=ERZ_FELSO_EV, mentes=False)
            kreg_erz[_v] = _d.set_index('Country Name')[LEPES_OSZLOP[valtozo_lepes[_v]]].reindex(kreg_erz.index)
            _sor.update({'Ellenőrzés éve': int(_d['Ev'].iloc[0]), 'Megjegyzés': 'korábbi évvel újraszámolva'})
        except ValueError:
            _sor.update({'Ellenőrzés éve': _fo_ev,
                         'Megjegyzés': f'nincs {ERZ_ALSO_EV}–{ERZ_FELSO_EV} közötti adat (csak {_fo_ev}) – nem változik'})
    ev_sorok_erz.append(_sor)
ev_tabla_erz = pd.DataFrame(ev_sorok_erz)
print(ev_tabla_erz.to_string(index=False))

# Ugyanaz a modell a módosított adatokon. A mintát ugyanúgy képezzük, mint a fő modellnél (szel_K változók, listaszerű törlés).
minta_erz = kreg_erz[[K_CEL] + szel_K].dropna()
if 'fejlettsegi_fokomponens' in Xv_K:
    X_erz, _ = pca_terv_K(minta_erz, minta_erz, [v for v in Xv_K if v != 'fejlettsegi_fokomponens'])
    X_erz = X_erz[Xv_K]
else:
    X_erz = minta_erz[Xv_K]
modell_erz = ols(minta_erz[K_CEL], X_erz)

# Összevetés: a fő modell, a fő modell a közös (mindkét adatváltozatban szereplő) országokon, és az ellenőrző modell
_kozos = minta_Ksz.index.intersection(minta_erz.index)
modell_kozos = ols(minta_Ksz.loc[_kozos, K_CEL], X_K.loc[_kozos])
erz_tabla = pd.DataFrame({
    'Változó': [VALTOZO_CIMKE.get(v, v) for v in Xv_K],
    'Fő modell: B': modell_K.params[Xv_K].values, 'Fő modell: p': modell_K.pvalues[Xv_K].values,
    'Fő modell, közös minta: B': modell_kozos.params[Xv_K].values, 'Fő modell, közös minta: p': modell_kozos.pvalues[Xv_K].values,
    'Korábbi adatévek: B': modell_erz.params[Xv_K].values, 'Korábbi adatévek: p': modell_erz.pvalues[Xv_K].values,
})
erz_tabla['Azonos előjel'] = np.where(np.sign(erz_tabla['Fő modell: B']) == np.sign(erz_tabla['Korábbi adatévek: B']), 'igen', 'nem')
erz_tabla['Szignifikáns mindkettőben'] = np.where((erz_tabla['Fő modell: p'] < SZIGNIF) & (erz_tabla['Korábbi adatévek: p'] < SZIGNIF), 'igen', 'nem')
print(f'   Fő modell: N = {int(modell_K.nobs)}, R² = {modell_K.rsquared:.3f} | közös minta: N = {int(modell_kozos.nobs)}, '
      f'R² = {modell_kozos.rsquared:.3f} | korábbi adatévek: N = {int(modell_erz.nobs)}, R² = {modell_erz.rsquared:.3f}')
print(erz_tabla.round(3).to_string(index=False))
if 'fejlettsegi_fokomponens' in Xv_K:
    _ci = modell_erz.conf_int().loc['fejlettsegi_fokomponens']
    print(f'   Összefoglalva: a fejlettségi-ellátási főkomponens együtthatója korábbi adatévekkel {modell_erz.params["fejlettsegi_fokomponens"]:.2f} '
          f'(95% CI: {_ci[0]:.2f} – {_ci[1]:.2f}; p = {modell_erz.pvalues["fejlettsegi_fokomponens"]:.4f}; N = {int(modell_erz.nobs)}), '
          f'a fő modellben {modell_K.params["fejlettsegi_fokomponens"]:.2f}.')

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 5c. A MODELLBŐL KIHAGYOTT VÁLTOZÓK ELLENŐRZÉSE: mindegyiket egyenként hozzáadjuk a választott modellhez, és ugyanazon a
# mintán (az adott változóval is rendelkező országokon) összevetjük a modellt a változóval és nélküle.
KIHAGYOTT_VALTOZOK_K = {
    'cuk_kiadas_teljes': 'átfed az egy cukorbetegre jutó kiadással',
    'orvosok_szama': 'azonos az orvosok 10 000 lakosra jutó számával',
    'korhazi_agyak': f'lefedettség < {MIN_LEFEDETTSEG:.0%}',
}
VALTOZO_CIMKE.update({'cuk_kiadas_teljes': 'Cukorbetegséggel kapcsolatos teljes kiadás / lakos',
                      'orvosok_szama': 'Orvosok száma (népességarányosítva)'})
terv_K = versengo_K[valasztott_K]

def terv_kiegeszitve(v):
    def terv(tanito, teszt):
        a, b = terv_K(tanito, teszt)
        return a.join(tanito[[v]]), b.join(teszt[[v]])
    return terv

kihagyott_sorok = []
for _v, _ok in KIHAGYOTT_VALTOZOK_K.items():
    _s = modell_valtozo(_v)[0]
    _d = kreg[[K_CEL] + szel_K].join(_s, how='left').dropna()
    _X0, _ = terv_K(_d, _d)
    _X1, _ = terv_kiegeszitve(_v)(_d, _d)
    _m0, _m1 = ols(_d[K_CEL], _X0), ols(_d[K_CEL], _X1)
    _sor = {'Változó': VALTOZO_CIMKE.get(_v, _v), 'Kihagyás oka': _ok,
            'Lefedettség (célváltozós országok)': _s.reindex(kreg.index[van_K]).notna().mean(), 'N': len(_d),
            'p (HC3)': _m1.pvalues[_v], 'ΔR²': _m1.rsquared - _m0.rsquared,
            'CV-RMSE nélküle (%pont)': cv_rmse(terv_K, _d, K_CEL)[0],
            'CV-RMSE vele (%pont)': cv_rmse(terv_kiegeszitve(_v), _d, K_CEL)[0]}
    if 'fejlettsegi_fokomponens' in _X0.columns:
        _sor['Főkomponens B nélküle'] = _m0.params['fejlettsegi_fokomponens']
        _sor['Főkomponens B vele'] = _m1.params['fejlettsegi_fokomponens']
    kihagyott_sorok.append(_sor)
kihagyott_K = pd.DataFrame(kihagyott_sorok)
_d_kiadas = kreg[[K_CEL] + szel_K].join(modell_valtozo('cuk_kiadas_teljes')[0], how='left').dropna()
_r_kiadas = _d_kiadas[['cuk_kiadas_teljes', 'cuk_egeszs_kiad_egy_fore']].corr().iloc[0, 1]
print('\n5c. A modellből kihagyott változók ellenőrzése (egyenként hozzáadva a választott modellhez):')
print(kihagyott_K.round(3).to_string(index=False))
print(f'   A teljes kiadás / lakos és az egy cukorbetegre jutó kiadás korrelációja a modellmintán: r = {_r_kiadas:.2f}')

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 6. MODELLTÍPUS
modelltipus_valasztas('K', 'kezeletlen cukorbetegek aránya (NCD-RisC)', 'arany', X_K, minta_Ksz[K_CEL])
modelltipus_K = modelltipus_osszevetesek[-1]

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 7. ELLENŐRZÉS AZ IDF NEM DIAGNOSZTIZÁLT ARÁNYÁVAL
# Elvárás: pozitív kapcsolat (több felismeretlen beteg -> több kezeletlen), és ugyanazok az okok, azonos előjellel.
print('\n7. Ellenőrzés az IDF nem diagnosztizált arányával')
valid_K = []
for nev, maszk in [('minden ország', kreg['aluldiag_arany'].notna()), ('mért országok', kreg['potlas_statusz'] == 'mért'),
                   ('pótolt országok', kreg['potlas_statusz'] == 'pótolt')]:
    d = kreg.loc[maszk, [K_CEL, 'aluldiag_arany']].dropna()
    if len(d) > 10:
        r, p = stats.pearsonr(d[K_CEL], d['aluldiag_arany'])
        valid_K.append({'Minta': nev, 'N': len(d), 'Pearson r': r, 'p': p, 'Spearman ρ': stats.spearmanr(d[K_CEL], d['aluldiag_arany'])[0]})
valid_K = pd.DataFrame(valid_K)
print(valid_K.round(3).to_string(index=False))
# Ugyanaz a modell (ugyanazok a változók) az IDF-arányra, csak a mért országokon
_mert_K = X_K.index[kreg.loc[X_K.index, 'potlas_statusz'] == 'mért']
_y_idf = kreg.loc[_mert_K, 'aluldiag_arany']
m_idf = ols(_y_idf, X_K.loc[_mert_K])
okok_osszevetes = pd.DataFrame({
    'Változó': [VALTOZO_CIMKE.get(v, v) for v in Xv_K],
    'K: std. béta': (modell_K.params[Xv_K] * X_K.std() / minta_Ksz[K_CEL].std()).values,
    'K: p': modell_K.pvalues[Xv_K].values,
    'IDF (mért orsz.): std. béta': (m_idf.params[Xv_K] * X_K.loc[_mert_K].std() / _y_idf.std()).values,
    'IDF (mért orsz.): p': m_idf.pvalues[Xv_K].values})
okok_osszevetes['Azonos irány'] = np.where(np.sign(okok_osszevetes['K: std. béta']) == np.sign(okok_osszevetes['IDF (mért orsz.): std. béta']), 'igen', 'nem')
print(f'   Ugyanaz a modell az IDF nem diagnosztizált arányára, a mért országokon (N = {int(m_idf.nobs)}, R² = {m_idf.rsquared:.3f}):')
print(okok_osszevetes.round(3).to_string(index=False))

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 8. ORSZÁGOK: a vártnál rosszabb / jobb helyzet (a modell szerinti várható értékhez képest)
orszag_K = minta_Ksz[[K_CEL]].copy()
orszag_K['Várható (modell)'] = modell_K.fittedvalues
orszag_K['Eltérés (%pont)'] = modell_K.resid
orszag_K['Regio'] = kreg.loc[orszag_K.index, 'Regio']
orszag_K['IDF nem diagnosztizált (%)'] = kreg.loc[orszag_K.index, 'aluldiag_arany']
print('\n8. A vártnál LÉNYEGESEN TÖBB kezeletlen cukorbeteg:')
print(orszag_K.sort_values('Eltérés (%pont)', ascending=False).head(10).round(1).to_string())
print('   A vártnál KEVESEBB (jó gyakorlat):')
print(orszag_K.sort_values('Eltérés (%pont)').head(10).round(1).to_string())
if 'Hungary' in kreg.index:
    _hu = kreg.loc['Hungary']
    print(f'   Magyarország: kezeletlen arány {_hu[K_CEL]:.1f}% (95% bizonytalansági sáv szélessége {_hu["CI_szelesseg"]:.1f} %pont); '
          f'az országok {(kreg[K_CEL] < _hu[K_CEL]).mean() * 100:.0f}%-ánál magasabb'
          + (f'; a modell szerint várható {orszag_K.loc["Hungary", "Várható (modell)"]:.1f}%, eltérés {orszag_K.loc["Hungary", "Eltérés (%pont)"]:+.1f} %pont'
             if 'Hungary' in orszag_K.index else '; a modell mintájában nem szerepel (hiányzó magyarázó változó)'))

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 9. MEGOLDÁSI FORGATÓKÖNYVEK – a VÁLASZTOTT (K) modellel
# Nem külön forgatókönyv-modell készül: a választott modell (4. lépés) frakcionális logit változatát
# használjuk (a 6. lépés ezt a modelltípust jelölte ki), ugyanazon a mintán (minta_Ksz). A beavatkozható tényezők (UHC-index,
# orvosellátottság) a fejlettségi főkomponens részei: a forgatókönyvben ezeket emeljük meg a saját IDF-régió felső negyedéig
# (aki fölötte van, annál nincs változás), a főkomponenst az eredeti standardizálással és súlyokkal újraszámoljuk, a modell
# többi változója változatlan. Ha a választott modellben nincs főkomponens, a beavatkozási változók közvetlenül módosulnak.
# Eredmény: a kezeletlen arány várható csökkenése és a kezelésbe kerülő cukorbetegek száma (20–79 éves cukorbetegek, IDF –
# közelítés, mert a kezelési arány 30+ évesekre vonatkozik). Bizonytalanság: paraméteres szimuláció (az együtthatók becsült,
# többdimenziós normális eloszlásából SZIMULACIO_N húzás).
beav = [v for v in FORGATOKONYV_BEAVATKOZAS if v in szel_K]
fg_adat = kreg.loc[minta_Ksz.index, [K_CEL, 'Regio', 'cukorbetegek_fo']].join(minta_Ksz[szel_K])
Xf = sm.add_constant(X_K, has_constant='add')
fg_modell = sm.GLM(fg_adat[K_CEL] / 100, Xf, family=sm.families.Binomial()).fit(cov_type='HC0')
print(f'\n9. Megoldási forgatókönyvek – a választott modell frakcionális logit változatával (N = {int(fg_modell.nobs)}):')
print(pd.DataFrame({'Változó': [VALTOZO_CIMKE.get(v, v) for v in X_K.columns], 'Együttható': fg_modell.params[list(X_K.columns)].values,
                    'p': fg_modell.pvalues[list(X_K.columns)].values}).round(4).to_string(index=False))
regio_cel = kreg.groupby('Regio')[beav].quantile(FORGATOKONYV_KVANTILIS)
_fg_pc = 'fejlettsegi_fokomponens' in X_K.columns
if _fg_pc:
    _fg_sk = StandardScaler().fit(minta_Ksz[fejl_K])
    _fg_pca = PCA(n_components=1).fit(_fg_sk.transform(minta_Ksz[fejl_K]))
    _fg_el = np.sign(_fg_pca.components_[0][fejl_K.index('gdp_per_capita')]) if 'gdp_per_capita' in fejl_K else 1
    _ell = _fg_el * _fg_pca.transform(_fg_sk.transform(minta_Ksz[fejl_K]))[:, 0]
    print(f'   Ellenőrzés: az újraszámolt főkomponens eltérése a modellbelitől: {np.abs(_ell - X_K["fejlettsegi_fokomponens"].values).max():.2e}')

def forgatokonyv_X(valtozok):
    nyers = fg_adat[szel_K].copy()
    for v in valtozok:
        cel = fg_adat['Regio'].map(regio_cel[v])
        nyers[v] = np.maximum(nyers[v], cel.fillna(nyers[v]))
    X = Xf.copy()
    for v in valtozok:
        if v in X.columns:
            X[v] = nyers[v]
    if _fg_pc:
        X['fejlettsegi_fokomponens'] = _fg_el * _fg_pca.transform(_fg_sk.transform(nyers[fejl_K]))[:, 0]
    erintett = (nyers[valtozok] != fg_adat[valtozok]).any(axis=1)
    return X, erintett

sim = np.random.default_rng(VELETLEN_MAG).multivariate_normal(fg_modell.params.values, fg_modell.cov_params().values, SZIMULACIO_N)
p_alap = 1 / (1 + np.exp(-(Xf.values @ sim.T)))   # országok × szimulációk
forgatokonyvek = {f'{VALTOZO_CIMKE.get(v, v)} a régió felső negyedéig': [v] for v in beav}
if len(beav) > 1:
    forgatokonyvek['Mindkettő együtt'] = beav
fg_sorok, fg_orszag = [], pd.DataFrame(index=fg_adat.index)
for nev, valt in forgatokonyvek.items():
    X_uj, erintett = forgatokonyv_X(valt)
    p_uj = 1 / (1 + np.exp(-(X_uj.values @ sim.T)))
    delta_pp = (p_alap - p_uj) * 100                            # a kezeletlen arány csökkenése (%pont), országonként
    fo = delta_pp / 100 * fg_adat['cukorbetegek_fo'].fillna(0).values[:, None]
    fg_orszag[f'{nev}: csökkenés (%pont)'] = delta_pp.mean(axis=1)
    fg_orszag[f'{nev}: kezelésbe kerülők (fő)'] = fo.mean(axis=1)
    sor = {'Forgatókönyv': nev, 'Érintett országok': int(erintett.sum()),
           'Átlagos csökkenés az érintett országokban (%pont)': delta_pp[erintett.values].mean() if erintett.any() else 0,
           'Kezelésbe kerülők összesen (millió fő)': fo.sum(axis=0).mean() / 1e6,
           'CI alsó (millió fő)': np.percentile(fo.sum(axis=0), 2.5) / 1e6,
           'CI felső (millió fő)': np.percentile(fo.sum(axis=0), 97.5) / 1e6}
    if 'Hungary' in fg_adat.index:
        _i = fg_adat.index.get_loc('Hungary')
        sor['Magyarország: csökkenés (%pont)'] = delta_pp[_i].mean()
        sor['Magyarország: kezelésbe kerülők (ezer fő)'] = fo[_i].mean() / 1000
    fg_sorok.append(sor)
forgatokonyv_tabla = pd.DataFrame(fg_sorok)
print(forgatokonyv_tabla.round(2).to_string(index=False))
fg_orszag['Regio'] = fg_adat['Regio']
forgatokonyv_regiok = fg_orszag.groupby('Regio')[[c for c in fg_orszag.columns if c.endswith('(fő)')]].sum() / 1e6
forgatokonyv_regiok.columns = [c.replace('(fő)', '(millió fő)') for c in forgatokonyv_regiok.columns]
print('   Kezelésbe kerülők IDF-régiónként (millió fő):')
print(forgatokonyv_regiok.round(2).to_string())
_fg_hatas = 'fejlettsegi_fokomponens' if _fg_pc else beav[0]
if fg_modell.pvalues[_fg_hatas] >= SZIGNIF:
    print(f'   FIGYELEM: a forgatókönyvet hordozó tényező ({VALTOZO_CIMKE.get(_fg_hatas, _fg_hatas)}) nem szignifikáns; a forgatókönyv csak tájékoztató.')

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ÁBRÁK
# a) A K modell és ugyanaz a modell az IDF-arányra: standardizált béták egymás mellett (egyezik-e az okok iránya?)
_b = okok_osszevetes.set_index('Változó')[['K: std. béta', 'IDF (mért orsz.): std. béta']]
_b = _b.reindex(_b['K: std. béta'].abs().sort_values().index)
_b.columns = ['Kezeletlen arány (K modell, NCD-RisC)', 'Nem diagnosztizált arány (IDF, mért országok)']
fig, ax = plt.subplots(figsize=(9, 1.5 + 0.6 * len(_b)))
_b.plot.barh(ax=ax, color=['#c0392b', '#7f8c8d'])
ax.axvline(0, color='black', linewidth=0.8)
ax.set_xlabel('Standardizált béta')
ax.set_ylabel('')
ax.set_title('Az okok két különböző forrás alapján')
_h, _l = ax.get_legend_handles_labels()
ax.legend(_h[::-1], _l[::-1], fontsize=8)
plt.tight_layout()
plt.savefig('regresszio_K_okok_ket_forrasban.png', dpi=300)
plt.show()

# b) Kezeletlen és nem diagnosztizált arány országonként
fig, ax = plt.subplots(figsize=(8, 6))
for _st, _szin in [('mért', '#1f77b4'), ('bizonytalan', '#ff7f0e'), ('pótolt', '#7f7f7f')]:
    _r = kreg[kreg['potlas_statusz'] == _st]
    ax.scatter(_r['aluldiag_arany'], _r[K_CEL], s=20, alpha=0.75, color=_szin, label=f'{_st} ({len(_r)})')
if 'Hungary' in kreg.index:
    ax.annotate('Magyarország', (kreg.loc['Hungary', 'aluldiag_arany'], kreg.loc['Hungary', K_CEL]), fontsize=8)
ax.set_xlabel('IDF: nem diagnosztizált arány (%)')
ax.set_ylabel('NCD-RisC: kezeletlen arány (%)')
ax.set_title('Kezeletlen és nem diagnosztizált cukorbetegek aránya országonként')
ax.legend(title='IDF-adat státusza', fontsize=8)
plt.tight_layout()
plt.savefig('regresszio_K_idf_osszevetes.png', dpi=300)
plt.show()

# c) Forgatókönyvek régiónként
_osz = [c for c in forgatokonyv_regiok.columns if c.startswith('Mindkettő')] or list(forgatokonyv_regiok.columns[:1])
_r = forgatokonyv_regiok[_osz[0]].sort_values()
fig, ax = plt.subplots(figsize=(9, 5))
ax.barh(_r.index, _r.values, color='#2a7ab0')
ax.set_xlabel('Kezelésbe kerülő cukorbetegek (millió fő)')
ax.set_title(f'Forgatókönyv: {_osz[0].split(":")[0]} – IDF-régiónként')
plt.tight_layout()
plt.savefig('regresszio_K_forgatokonyv_regiok.png', dpi=300)
plt.show()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 10. MENTÉS
with pd.ExcelWriter('regresszio_kezeletlen_arany.xlsx') as writer:
    pd.DataFrame([{'Célváltozó': 'kezeletlen cukorbetegek aránya (%, NCD-RisC, korstandardizált, 30+)', 'Választott modell': valasztott_K,
                   'N': int(modell_K.nobs), 'R²': modell_K.rsquared, 'Korrigált R²': modell_K.rsquared_adj,
                   'CV-R²': versengo_tabla_K.loc[versengo_tabla_K['Választott'] == 'igen', 'CV-R²'].iloc[0],
                   'Modelltípus': modelltipus_valasztasok[-1]['Választott modelltípus'], **diag_K}]).to_excel(writer, sheet_name='Osszefoglalo', index=False)
    blokk_K.to_excel(writer, sheet_name='Okcsoportok', index=False)
    szures_K.reset_index(names='Kód').to_excel(writer, sheet_name='Valtozoszures', index=False)
    versengo_tabla_K.to_excel(writer, sheet_name='Versengo_modellek', index=False)
    egyutthatok_K.to_excel(writer, sheet_name='Valasztott_modell', index=False)
    rob_tabla_K.to_excel(writer, sheet_name='Robusztussag')
    modelltipus_K.to_excel(writer, sheet_name='Modelltipus', index=False)
    valid_K.to_excel(writer, sheet_name='Ellenorzes_IDF', index=False)
    okok_osszevetes.to_excel(writer, sheet_name='Okok_K_vs_IDF', index=False)
    orszag_K.sort_values('Eltérés (%pont)', ascending=False).to_excel(writer, sheet_name='Orszagok')
    forgatokonyv_tabla.to_excel(writer, sheet_name='Forgatokonyvek', index=False)
    forgatokonyv_regiok.to_excel(writer, sheet_name='Forgatokonyv_regiok')
    fg_orszag.to_excel(writer, sheet_name='Forgatokonyv_orszagok')
    ev_tabla_erz.to_excel(writer, sheet_name='Adatev_valtozok', index=False)
    erz_tabla.to_excel(writer, sheet_name='Adatev_erzekenyseg', index=False)
    kihagyott_K.to_excel(writer, sheet_name='Kihagyott_valtozok', index=False)
print('\nMentve: regresszio_kezeletlen_arany.xlsx, regresszio_K_okok_ket_forrasban.png, regresszio_K_idf_osszevetes.png, '
      'regresszio_K_forgatokonyv_regiok.png')


# =========================================================================================================================================================================
# MODELLTÍPUS-VÁLASZTÁS – ÖSSZEFOGLALÓ ÉS MENTÉS (B és K modell)
# =========================================================================================================================================================================

modelltipus_tabla = pd.DataFrame([{k: v for k, v in r.items() if k != 'Együtthatók'} for r in modelltipus_valasztasok])
print('\n' + '=' * 100 + '\nMODELLTÍPUS-VÁLASZTÁS ÖSSZEFOGLALÓJA\n' + '=' * 100)
print(modelltipus_tabla[['Modell', 'Típus', 'N', 'Választott modelltípus', 'Kiemelt hatás']].to_string(index=False))
with pd.ExcelWriter('regresszio_modelltipus.xlsx') as writer:
    modelltipus_tabla.to_excel(writer, sheet_name='Valasztasok', index=False)
    pd.concat(modelltipus_osszevetesek).to_excel(writer, sheet_name='Jeloltek_CV_hibaja', index=False)
    pd.concat([r['Együtthatók'] for r in modelltipus_valasztasok]).to_excel(writer, sheet_name='Valasztott_egyutthatok', index=False)
print('Mentve: regresszio_modelltipus.xlsx')