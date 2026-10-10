# Plan: cold agents (uzivatelsky test API a dokumentace agenty bez kontextu)

Stav: navrh (2026-10-10), neschvaleny. Vychozi verze: eis_analysis v0.58.0.

## Kontext

Inspirace: driftlet (`doc/DRIFTLET_COMPARISON1.md`, bod 3). Verze 0.15.0,
0.15.1 a 0.16.1 driftletu vzesly hlavne z kol "cold agents": cerstvy agent
bez kontextu vyvoje dostane realnou ulohu se znamou teoretickou odpovedi
a k dispozici ma jen pruvodce a dokumentaci. Kazde kolo naslo tiche spatne
vysledky (0.15.0: pet, napr. o 10-25 % vic fotoproudu) a mista, kde
dokumentace nebo API sveda k chybe.

U nas to nic nepokryva. Unit testy a stress test (`doc/STRESS_TEST.md`)
overuji algoritmy volane spravne. Nikdo neoveruje, jestli se z README
a `doc/PYTHON_API.md` da spravne volani vubec sestavit: jaky obvod zvolit,
kterou funkci zavolat, jak cist vysledek, co znamena varovani.

**Cil:** zjistit, kde nas uzivatel (clovek nebo LLM agent) s dokumentaci
a API dojde ke spatnemu vysledku nebo ztrati cas, a kazdy nalez overit
a roztridit. Vystupem jsou nalezy a navrhy oprav, ne opravy samotne.

**Neni cilem:** automatizace, CI, testovaci framework. Kolo se spousti
rucne, nalezy se zapisuji do dokumentu. Kod v repozitari nevznika (krome
pripadnych oprav, kazda samostatne a po schvaleni).

## Principy

1. **Studeny start.** Agent nedostane nic z historie vyvoje. Smi cist jen
   `README.md`, `doc/PYTHON_API.md` a `python3 eis.py --help`. Zdrojaky
   `eis_analysis/`, `tests/`, ostatni `doc/` a `example/` jsou zakazane.
   Izolace je mekka (instrukce v zadani): agent ma pristup k celemu disku.
   Proto v debriefu vypise vsechny soubory, ktere precetl. Porusi-li zakaz,
   vysledek ulohy se nepocita, nalezy z jeho debriefu jen jako podnet.
2. **Realna uloha se znamou pravdou.** Zadani je formulovano jako od
   uzivatele ("urci tloustku oxidu"), ne jako test API ("zavolej
   `analyze_oxide_layer`"). Pravda je znama (synteticka data z daneho
   obvodu, nebo realny soubor s drive overenym vysledkem) a ma toleranci.
   Bez pravdy agent chybu nepozna a vyrobi jen nejake cislo.
3. **Pravda mimo dosah agenta.** Data dostane jako kopii v pracovnim
   adresari mimo repozitar (scratchpad), s neutralnim nazvem a **bez
   hlavicky s parametry** (CSV v `example/` maji pravdu v komentari na
   prvnim radku). Pravda a tolerance jsou jen v tomto planu, ktery agent
   cist nesmi.
4. **Debrief je dulezitejsi nez vysledek.** Agent, ktery narazi, si to
   obvykle obejde a pokracuje. Pevna sablona debriefu (nize) ho nuti
   priznat, kde hadal a co obchazel.
5. **Kazdy nalez overit.** Agent se muze splest sam. Pred zapisem nalezu
   ho reprodukuji (vlastnim volanim nebo ctenim dokumentace); co se
   nereprodukuje, je chyba agenta.
6. **Metodu nepredepisovat.** Volba metody (DE vs. LM, typ obvodu, lambda)
   je soucast pozorovani: vede dokumentace k DE u oxidu, jak ho fitujeme
   my? Jen kde by volba metody ulohu znehodnotila, je v zadani.

## Mechanismus

- Agent = subagent Claude Code (`general-purpose`), kazdy spusten zvlast,
  bez sdileneho kontextu. Zadani je samostatny text: role, uloha, data,
  pravidla cteni, sablona debriefu.
- Model: vychozi. Volitelne jeden agent na slabsim modelu (`sonnet` nebo
  `haiku`) na stejnou ulohu: simuluje mene zkuseneho uzivatele, a co
  zvladne slabsi model, je dobre zdokumentovane.
- Pracovni adresar agenta: vlastni podadresar scratchpadu s daty; skripty
  pise tam. Repozitar nemeni (`isolation: worktree` neni potreba, agent nema
  duvod repozitar editovat; pokud ho zmeni, `git status` to ukaze).
- Paralelne nejvyse 4 agenti (DE fity jsou CPU-heavy, pravidlo max. ctyr
  procesu).
- Jazyk zadani: cesky, jako by psal uzivatel. Debrief cesky.

### Sablona debriefu (soucast zadani)

```
1. Vysledek: hodnota(y), nejistota, jak jsi ji urcil.
2. Postup: kroky a volani (funkce, parametry, CLI prikazy), v poradi.
3. Prectene soubory: uplny seznam (cesty).
4. Kde jsi tapal nebo hadal: misto v dokumentaci (soubor, sekce), co jsi
   cekal, co ses dozvedel, co jsi nakonec predpokladal.
5. Obchvaty: co nefungovalo napoprve a jak jsi to obesel (vcetne vyjimek,
   chybovych hlasek a varovani, ktera jsi videl, doslovne).
6. Nesoulad: tvrzeni v dokumentaci, ktera neodpovidala chovani.
7. Duvera: veris vysledku? Proc ano/ne?
```

## Ulohy kola 1

Ctyri ulohy pokryvajici hlavni cesty (fit + oxid, DRT, validace, CLI).
Pravdu a tolerance doplnit v kroku 1 (zde jen zdroj a smysl).

| Id | Zadani (zkracene) | Data | Pravda | Co testuje |
|---|---|---|---|---|
| T1 | "Mereni ZrO2 vrstvy, plocha A cm^2, eps_r = 22. Urci tloustku vrstvy a jeji nejistotu." | stress `anomalous`, vetev ZrO2 (`R-(Wa/Wat|Q|C)`), vybrany pripad se sumem 1 % | d = eps0 eps_r A / C z pravdiveho C | volba obvodu (Wa, C vedle Q), DE, `analyze_oxide_layer` a vyber dominantniho elementu |
| T2 | "Kolik procesu je ve spektru, jake maji casove konstanty a odpory?" | stress `rc` nebo `cpe`, 2-3 oblouky, sum 0.3-1 % | tau a R oblouku, R_pol | DRT (auto lambda, piky, GMM), cteni `DRTResult`, R_inf |
| T3 | "Jsou data duveryhodna? Pokud ano, vyhodnot DRT." | `example/EISPOT-test1.DTA` (kopie, prejmenovana) | kapacitni NF konec (Im = -1e6 Ohm pri 3 mHz), DRT bez seriove C nevhodne (`doc/DRIFTLET_COMPARISON.md`, `doc/DRT_RINF_L_ANALYSIS_2026-09-25.md`) | Lin-KK/Z-HIT, `--kk-series-c`, varovani DRT o pile-upu; pozna agent, ze vysledek DRT je nafouknuty? |
| T4 | "Jen pres CLI: nafituj Randles s difuzi, exportuj vysledky." | `example/diffusion_RpWs_noise1.csv` bez hlavicky (pravda Rs = 10, Rp = 100, R_D = 100 Ohm, tau = 1 s) | parametry z hlavicky | `eis.py --help`, syntaxe obvodu, vyber Ws/Wo/W, export JSON/CSV |

Pripadne T5 (po pilotu, pokud zbyde): Warburgovo spektrum Ag | AgNO3 | Ag
z driftletu s analytickou pravdou (`doc/DRIFTLET_COMPARISON.md`, mereni A),
jednotky Ohm*m^2 -> narazi agent na absolutni meze?

## Kroky

### 1. Priprava uloh a pravdy

- Vybrat konkretni pripady (`family/index`) pro T1 a T2: v okne viditelne
  vsechny elementy, sum 1 %, bez L. Zapsat obvod, parametry, frekvencni
  rozsah do tohoto planu.
- Pro T1 urcit A a pravdive d; pro T3 shrnout drive overeny stav.
- Pro kazdou ulohu tolerance: co je "spravne" (napr. d v ramci 2x stderr
  fitu C, tau oblouku v ramci 0.15 dekady jako Epeak).
- Vygenerovat data do scratchpadu (CSV `freq, Z_re, Z_im`, bez komentaru).
- **Overeni:** ulohu vyresim sam pres verejne API podle dokumentace
  a dostanu pravdu v toleranci. Neda-li se to, uloha je spatne polozena
  (nebo je to uz nalez) a upravi se pred spustenim agentu.

### 2. Pilot

- Jeden agent na T4 (nejkratsi). Ucel: odladit zadani a sablonu (rozumi
  agent pravidlum cteni, vyplni debrief, nedela neco jineho).
- **Overeni:** debrief vyplneny ve vsech bodech, seznam prectenych souboru
  bez zakazanych. Jinak upravit zadani a pilot opakovat.

### 3. Kolo 1

- T1-T4, kazda jednim agentem, nejvyse 4 soucasne. Volitelne druhy agent
  na T1 (nejslozitejsi) na slabsim modelu.
- Ulozit kazdy debrief doslovne (scratchpad, pak shrnuti do vysledku).

### 4. Triage

Kazdy nalez overit (princip 5) a zaradit:

| Kategorie | Priklad | Priorita |
|---|---|---|
| A: tichy spatny vysledek | knihovna vrati spatne cislo bez varovani pri rozumnem volani | vysoka (CLAUDE.md: chyby ovlivnujici vysledky) |
| B: nespravna dokumentace | priklad nefunguje, popis neodpovida chovani | vysoka (CLAUDE.md: nespravna dokumentace) |
| C: tren v API / chybejici informace | agent hadal, obchazel, varovani nepojmenovalo pricinu | stredni |
| D: chyba agenta | dokumentace je jasna, agent ji prehledl | bez akce; opakuje-li se u vice agentu, je to C |

Vystup: `doc/COLD_AGENT.md` (jako `doc/STRESS_TEST.md`): kolo, ulohy,
vysledek proti pravde, nalezy s kategorii a dukazem (citace debriefu,
reprodukce), navrzena akce.

### 5. Opravy

Mimo tento plan: ke kazdemu nalezu A-C navrh opravy, schvaleni,
samostatny commit (`fix(...)` nebo `docs(...)`). Nalezy B jsou vetsinou
jen dokumentace (bez bumpu verze).

### 6. Kolo 2 (overeni oprav)

Nove agenty (opet bez kontextu) na stejne ulohy, ktere mely nalezy A-C,
plus pripadne nove ulohy. Nalez je uzavren, kdyz ho druhe kolo nezopakuje.
Dale jen pri vetsich zmenach API nebo dokumentace (napr. pred minor
releasem), ne pravidelne.

## Rizika a omezeni

- **Mekka izolace:** agent muze precist zdrojaky nebo pravdu. Ochrana:
  data mimo repozitar bez hlavicek, pravda jen v tomto planu, povinny
  seznam prectenych souboru. Neda se vynutit; spolehnout se na debrief.
- **Falesne nalezy:** agent nahlasi chybu, ktera neni. Ochrana: princip 5.
- **Pokryti podle uloh:** agent najde jen to, kam ho uloha zavede. Ctyri
  ulohy pokryvaji hlavni cesty, ne vse (Voigt chain, THD, export grafu
  nepokryto zamerne).
- **Cena:** kazdy agent desitky minut a odpovidajici spotrebu tokenu; DE
  fit oxidu na RK3588 minuty. Kolo 1 = pilot + 4 agenti.
- **Agent jako uzivatel neni clovek:** LLM cte dokumentaci celou a doslova,
  clovek prehlizi. Co agent zvladne, nemusi zvladnout clovek; co agent
  nezvladne, clovek skoro jiste ne.

## Otevrene otazky (rozhodnuti uzivatele)

1. Sada uloh T1-T4 (a T5), nebo jine? Hlavne T1: chceme overit spis cestu
   k tloustce, nebo i k permitivite (eps z C pri znamem d)?
2. Ma agent smet cist cely `doc/`, nebo jen `README.md` a `PYTHON_API.md`?
   Realny uzivatel knihovny vidi cely adresar; uzsi vyber testuje hlavni
   dokumentaci prisneji.
3. Slabsi model pro jednu ulohu ano/ne.
