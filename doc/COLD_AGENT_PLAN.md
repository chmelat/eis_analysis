# Plan: cold agents (uzivatelsky test API a dokumentace agenty bez kontextu)

Stav: navrh v2 (2026-10-13), neschvaleny. Vychozi verze: eis_analysis v0.59.1.

## Kontext

Inspirace: driftlet (CHANGELOG 0.15.0 az 0.16.3). Vetsina jeho oprav
vzesla z kol "cold agents": cerstvy agent bez kontextu vyvoje dostane
ucebnicovou ulohu s analytickou odpovedi (pn prechod, GITT, smiseny
potencial) a k dispozici ma jen pruvodce pro agenty. Vysledky vetsinou
sedi s teorii; kazde kolo ale naslo tiche spatne vysledky (0.15.0: pet,
napr. o 10-25 % vic fotoproudu) a mista, kde dokumentace nebo API sveda
k chybe.

U nas to nic nepokryva. Unit testy a stress test (`doc/STRESS_TEST.md`)
overuji algoritmy volane spravne. Nikdo neoveruje, jestli se z README
a navazujici dokumentace da spravne volani vubec sestavit: kterou funkci
zavolat, s jakymi parametry, jak cist vysledek, co znamena varovani.

**Cil:** zjistit, kde nas uzivatel (clovek nebo LLM agent) s dokumentaci
a API dojde ke spatnemu vysledku nebo ztrati cas, a kazdy nalez overit
a roztridit. Vystupem jsou nalezy a navrhy oprav, ne opravy samotne.

**Neni cilem:** hledat model k experimentalnim datum (to je vyzkum, ne test
dokumentace; bez znamene pravdy nalez nejde overit), automatizace, CI,
testovaci framework. Kolo se spousti rucne, nalezy se zapisuji do
dokumentu. Kod v repozitari nevznika (krome pripadnych oprav, kazda
samostatne a po schvaleni); generator dat a pravda ziji v piskovisti.

## Principy

1. **Studeny start v piskovisti.** Agent bezi v adresari mimo repozitar,
   s vlastni instalaci balicku a kopii uzivatelske dokumentace (viz
   Piskoviste). Nedostane CLAUDE.md projektu ani pamet projektu: ta obsahuje
   vysledky vyvoje (napr. volbu obvodu pro ZrO2), ktere by ulohy
   kontaminovaly. Izolace od zbytku disku je mekka (instrukce v zadani),
   proto debrief zacina kontrolou (bod 0) a obsahuje seznam prectenych
   souboru.
2. **Umela uloha s analytickou pravdou.** Synteticka data ze zname
   nahradni soustavy, pravda plyne ze vzorce (d = eps0 eps_r A / C,
   tau = R C). Zadani zni jako od uzivatele ("urci tloustku vrstvy"), ne
   jako test API ("zavolej `analyze_oxide_layer`"). Bez pravdy agent chybu
   nepozna a vyrobi jen nejake cislo.
3. **Pravda mimo dosah agenta.** Data dostane jako CSV `freq, Z_re, Z_im`
   bez komentare (CSV v `example/` maji pravdu na prvnim radku). Parametry
   jsou cerstve vygenerovane, ne z `example/`. Pravda, tolerance a generator
   jsou v `~/coldbox-private/`, ne v tomto planu a ne v repozitari.
4. **Parametry jine nez vychozi.** Kde ma funkce vychozi hodnotu (eps_r = 22,
   plocha 1 cm^2), uloha pouzije jinou; jinak test nepozna agenta, ktery
   parametr vubec nepreda.
5. **Debrief je dulezitejsi nez vysledek.** Agent, ktery narazi, si to
   obvykle obejde a pokracuje. Pevna sablona debriefu ho nuti priznat, kde
   hadal a co obchazel.
6. **Kazdy nalez overit.** Agent se muze splest sam. Pred zapisem nalezu
   ho reprodukuji (vlastnim volanim nebo ctenim dokumentace); co se
   nereprodukuje, je chyba agenta.
7. **Metodu nepredepisovat.** Volba metody (DRT lambda, KK test, fit) je
   soucast pozorovani. Jen kde by volba metody ulohu znehodnotila, je
   v zadani (T4: jen CLI).

## Piskoviste

```
~/coldbox/                 svet agenta (jeho pracovni adresar)
  .venv/                   python3 -m venv; pip install <wheel z repozitare>
  docs/                    kopie README.md a vseho, kam z nej vedou relativni
                           odkazy (tranzitivne), ve stejne strukture
  runs/r1/T1/ ...          data a zadani.md; agent sem pise skripty
~/coldbox-private/         mimo svet agenta: setup.sh, generator, pravda,
                           tolerance, debriefy, poznamky k triage
```

`~/coldbox-private/setup.sh` piskoviste (pre)stavi z HEAD repozitare
(odmitne necommitnute zmeny) a zapise hash do `~/coldbox-private/COMMIT`.

- **Instalace:** wheel z repozitare (`pip wheel . --no-deps -w ...`), do venv
  nainstalovany ne-editable i se zavislostmi. Venv nevidi uzivatelsky
  site-packages, takze ani editable instalaci odkazujici do repozitare;
  `help()` a docstringy ma agent k dispozici jako kazdy uzivatel balicku
  (docstringy jsou podle CLAUDE.md reference API). Vedlejsi efekt: odpadne
  sum prostredi vyvoje (varovani matplotlibu o Axes3D, zastarala metadata
  editable instalace, ktera `pip show` hlasi jako 0.9.2).
- **Dokumentace:** pravidlo je mechanicke, bez rucniho vyberu: README
  a vse, kam z nej vedou odkazy, i pres dalsi dokumenty (v0.59.1:
  31 souboru vcetne CHANGELOG, obrazku a `example/example_eis_data.csv`).
  Tranzitivne, aby v kopii nebyly mrtve odkazy. Plany, srovnani, stress
  test a `tests/` agent nema (nikdo na ne neodkazuje). Pravidlo odpovida
  realnemu uzivateli: dokumentace, kam ho README posle.
- **Spusteni agenta:** bezhlavy `claude -p` z `runs/<kolo>/<uloha>/`, jen
  s uzivatelskymi nastavenimi (`--setting-sources user`), vystup do
  `~/coldbox-private/<kolo>/<uloha>.md`. Volitelne `--model sonnet`.
  Presne prepinace (povolene nastroje, rezim opravneni) overit v pilotu.
  Globalni `~/.claude/CLAUDE.md` se nacte i tak; neobsahuje nic o projektu.
- **Prostredi v zadani:** cesta k `~/coldbox/.venv/bin/python` a `eis`,
  `MPLBACKEND=Agg` (zadna okna s grafy).
- **Obnova pred kolem:** novy wheel a nova kopie docs z daneho commitu;
  hash commitu do vysledku kola.
- **Paralelne nejvyse 2 agenti** (stroj ma 4 vlakna; numpy/scipy a fit
  vlakna sdileji).

### Zadani agenta (kostra)

Role (vyzkumnik, ktery knihovnu vidi poprve), uloha, data, prostredi,
pravidla: cist jen `~/coldbox/docs/` a `help()` balicku, nevystupovat
z `~/coldbox/`, nic neinstalovat. Na konec debrief podle sablony.

### Osnova zpravy (soucast zadani)

Zprava o pouzitelnosti knihovny a dokumentace, ne o vlastnich uvahach
agenta: osnova ve tvaru "kde jsi tapal, co jsi cekal" spoustela filtr
Sonnetu (`[reasoning_extraction]`, kolo 1). Platne zneni je v
`~/coldbox/runs/r1/*/zadani.md`:

```
0. Znal jsi knihovnu eis_analysis uz pred touto ulohou (napr. z jineho
   projektu)?
1. Vysledek: hodnota(y) a nejistota, a z ceho plynou.
2. Pouzita volani: funkce s parametry nebo CLI prikazy, v poradi.
3. Prectene soubory: uplny seznam (cesty).
4. Nejasna mista v dokumentaci: soubor a sekce, co v ni chybelo nebo
   slo vylozit vic zpusoby, a kterou variantu jsi pouzil.
5. Problemy pri pouziti: vyjimky, chybove hlasky a varovani (doslovne)
   a jak jsi je vyresil.
6. Nesoulad: tvrzeni v dokumentaci, ktera neodpovidala chovani knihovny.
7. Spolehlivost vysledku: co ji podporuje a co zpochybnuje.
```

Bod 0 je kontrola izolace: zminka o vyvoji projektu (M136, Wa u ZrO2,
stress test) znamena, ze kontext unikl. Vysledek ulohy se pak nepocita.

## Ulohy kola 1

Ctyri umele ulohy pokryvajici hlavni cesty (oxid, DRT, validace, CLI).
Konkretni parametry, sum (0.5-1 %, seed) a tolerance se urci v kroku 1
a zapisou do `~/coldbox-private/`; zde jen typ a smysl. Vsechny prvky
v okne viditelne, procesy oddelene aspon o 2 dekady.

| Id | Zadani (zkracene) | Data | Pravda | Co testuje |
|---|---|---|---|---|
| T1 | "Dielektricka vrstva, plocha A, eps_r dane (ne 22, A ne 1 cm^2). Urci tloustku a jeji nejistotu." | R_s - (R_f \| C_f) | d = eps0 eps_r A / C_f | fit, `analyze_oxide_layer` (predani eps_r a plochy, jednotky), nejistota d z chyby fitu C (vysledek ji nema) |
| T2 | "Kolik procesu je ve spektru, jake maji casove konstanty a odpory? Jaky je odpor elektrolytu?" | R_s - (R1 \| C1) - (R2 \| C2), R1 != R2 | 2 procesy, tau_i = R_i C_i, R_i, R_pol, R_inf = R_s | DRT (auto lambda, piky), cteni `DRTResult`, R_inf |
| T3 | "Dve mereni tehoz vzorku. Ktere je duveryhodne? Pokud nektere jen zcasti, odkud kam?" | A: spektrum jako T2; B: totez s Re(Z) pod f_x zkreslenym (porusuje KK) | A v poradku, B nevalidni pod f_x | Lin-KK/Z-HIT, cteni rezidui a verdiktu, lokalizace vadneho pasma |
| T4 | "Jen pres CLI: nafituj Randlesuv obvod s konecnou difuzi (transmisivni hranice), vysledky uloz." | R_s - (R_p \| Ws), parametry jine nez v `example/` | parametry generatoru | `eis --help`, syntaxe obvodu, Ws vs Wo vs W, `--save` (JSON/CSV) |

Pozdeji (kolo 2 nebo dal): T1 s Q misto C (pravda C_eff podle
Hsu-Mansfelda; vsimne si agent rozptylu proti Brugovi?), inverzni smer
(eps_r ze zname tloustky), R_inf se seriovou L, Warburgovo spektrum
s analytickou pravdou z driftletu.

## Kroky

### 1. Piskoviste a pravda

- Zalozit `~/coldbox/` a `~/coldbox-private/` podle sekce Piskoviste.
- Generator v `~/coldbox-private/`: pro kazdou ulohu obvod, parametry,
  frekvencni rozsah, sum se seedem, tolerance (napr. d v ramci 2x stderr
  fitu C a nejvyse 5 %, tau v ramci 0.15 dekady jako Epeak, R v ramci 5 %).
- **Overeni piskoviste:** ve venv `eis_analysis.__file__` lezi ve venv
  a verze odpovida commitu; `eis --help` bezi.
- **Overeni uloh:** kazdou ulohu vyresim sam v piskovisti jen z
  `~/coldbox/docs/` a `help()`, a dostanu pravdu v toleranci. U T3 navic:
  nase validace B opravdu oznaci a A propusti. Neda-li se to, uloha je
  spatne polozena (nebo je to uz nalez) a upravi se pred spustenim agentu.

### 2. Pilot

- Jeden agent na T4 (nejkratsi). Ucel: odladit prepinace `claude -p`,
  zadani a sablonu.
- **Overeni:** bod 0 debriefu bez znalosti vyvoje; seznam prectenych
  souboru jen z `~/coldbox/`; debrief vyplneny ve vsech bodech. Jinak
  upravit a pilot opakovat. Bez splneni se kolo 1 nespousti.

### 3. Kolo 1

- T1-T4, kazda jednim agentem, nejvyse 2 soucasne. Navic T4 na slabsim
  modelu (`sonnet`): co zvladne slabsi model jen z `--help`, je dobre
  zdokumentovane.
- Kazdy debrief ulozit doslovne do `~/coldbox-private/r1/`.

### 4. Triage

Kazdy nalez overit (princip 6) a zaradit:

| Kategorie | Priklad | Priorita |
|---|---|---|
| A: tichy spatny vysledek | knihovna vrati spatne cislo bez varovani pri rozumnem volani | vysoka (CLAUDE.md: chyby ovlivnujici vysledky) |
| B: nespravna dokumentace | priklad nefunguje, popis neodpovida chovani | vysoka (CLAUDE.md: nespravna dokumentace) |
| C: tren v API / chybejici informace | agent hadal, obchazel, varovani nepojmenovalo pricinu | stredni |
| D: chyba agenta | dokumentace je jasna, agent ji prehledl | bez akce; opakuje-li se u vice agentu, je to C |

Vystup: `doc/COLD_AGENT.md` (jako `doc/STRESS_TEST.md`): kolo, commit,
ulohy, vysledek proti pravde, nalezy s kategorii a dukazem (citace
debriefu, reprodukce), navrzena akce. Pravda ulohy se tam zapise az po
kole (ulohy kola 2 jsou nove).

### 5. Opravy

Mimo tento plan: ke kazdemu nalezu A-C navrh opravy, schvaleni,
samostatny commit (`fix(...)` nebo `docs(...)`). Nalezy B jsou vetsinou
jen dokumentace (bez bumpu verze).

### 6. Kolo 2 (overeni oprav)

Nove agenty (opet bez kontextu, nove piskoviste z noveho commitu) na
ulohy, ktere mely nalezy A-C, s novymi parametry, plus pripadne nove
ulohy. Nalez je uzavren, kdyz ho druhe kolo nezopakuje. Dale jen pri
vetsich zmenach API nebo dokumentace (napr. pred minor releasem).

## Rizika a omezeni

- **Mekka izolace:** agent muze vyjit z `~/coldbox/` a najit repozitar.
  Ochrana: nema duvod (vse potrebne je v piskovisti), bod 0 a seznam
  prectenych souboru v debriefu. Neda se vynutit.
- **Falesne nalezy:** agent nahlasi chybu, ktera neni. Ochrana: princip 6.
- **Umele ulohy jsou lehci nez realita:** nepokryji volbu modelu ani
  realna data (zamerne, viz Cil). Nalezy jsou dolni mez.
- **Pokryti podle uloh:** agent najde jen to, kam ho uloha zavede (Voigt
  chain, THD, export grafu nepokryto).
- **Siroka dokumentace:** 22 souboru proti jednomu pruvodci driftletu
  (`llms.txt` se spustitelnymi sablonami). Ztraci-li se agenti v ni, je to
  nalez C; samostatny pruvodce pro agenty jen pokud to kolo ukaze.
- **Cena:** v kole 1 (Sonnet, umele ulohy) 1-1.5 min a 0.2-0.4 USD na
  agenta, cele kolo asi 1.5 USD (`doc/COLD_AGENT.md`).
- **Agent jako uzivatel neni clovek:** LLM cte dokumentaci celou a doslova,
  clovek prehlizi. Co agent zvladne, nemusi zvladnout clovek; co agent
  nezvladne, clovek skoro jiste ne.
