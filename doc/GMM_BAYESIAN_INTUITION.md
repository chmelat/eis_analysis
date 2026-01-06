# Bayesovske rozhodovani v DRT analyze s GMM

Intuitivni vysvetleni principu automatickeho urceni poctu piku v DRT spektru.

## Zakladni problem

Mas DRT spektrum (gamma vs tau) a chces zjistit: **Kolik tam je skutecnych elektrochemickych procesu?**

To je tezke, protoze:
- Piky se mohou prekryvat
- Sum vytvari falesne "picky"
- Nevis predem, jestli mas 1, 2, 3 nebo 5 procesu

## GMM jako model

Gaussian Mixture Model rika: "Predpokladam, ze DRT spektrum vzniklo sectenim N Gaussovskych zvonu."

```
gamma(tau) = w1*G1(mu1,sigma1) + w2*G2(mu2,sigma2) + ... + wn*Gn(mun,sigman)
```

Kde kazdy Gaussian reprezentuje jeden elektrochemicky proces.

## Proc log(tau) transformace?

GMM analyza se provadi v logaritmickem prostoru `log10(tau)` misto v linearnim `tau`. Duvody:

### Rozsah casovych konstant

Elektrochemicke procesy maji casove konstanty pres mnoho radu:

```
Priklad realneho systemu:
  - Prenos naboje:          tau1 = 0.001 s  (1 ms)
  - Difuze v oxidove vrstve: tau2 = 0.1 s   (100 ms)
  - Pomala adsorpce:        tau3 = 10 s    (10000 ms)
```

V linearnim prostoru tau:
```
tau:  |*|                    |*|                                        |*|
      0                      0.1                                        10
      ^
      tau1=0.001 je prakticky na nule!
```

V log10(tau) prostoru:
```
log(tau): |    *    |    *    |    *    |
         -3       -1        1

Rovnomerne rozlozene!
```

### Tvar piku

**V linearnim tau:** Piky jsou **asymetricke** (log-normalni distribuce)
- Dlouhy ocas smerem k vyssim tau
- GMM predpoklada Gaussiany -> spatny fit

**V log(tau):** Piky jsou **priblizne symetricke** (Gaussovske)
- GMM predpoklad je splnen
- Lepsi kvalita fitu

```
Linearni tau:                    Log(tau):
     ___                          ___
    /   \___                     /   \
   /        \___                /     \
  /             \____          /       \
--+-------------------> tau  --+---------> log(tau)
  Asymetricky                  Symetricky
```

### Fyzikalni interpretace sirky

Sirka piku v **dekadach** ma primy fyzikalni vyznam:

- **sigma ~ 0.2 dekad**: Idealni Voigt clanek (R||C), jedna casova konstanta
- **sigma > 0.5 dekad**: Distribuovany proces (CPE), vice casovych konstant

V linearnim prostoru by sigma zaviselo na poloze piku:
- Pik pri tau = 0.001 s by mel sigma v milisekundach
- Pik pri tau = 10 s by mel sigma v sekundach
- Nelze je primo porovnat!

### Separace prekryvajicich se piku

Dva procesy s tau1 = 1 ms a tau2 = 5 ms:

- **Linearni:** Vzdalenost = 4 ms (zavisi na absolutni hodnote)
- **Logaritmicky:** Vzdalenost = log10(5) - log10(1) = 0.7 dekad (nezavisle na skale)

GMM v log-prostoru detekuje prekryv konzistentne bez ohledu na to, jestli jsou piky pri ms nebo pri sekundach.

### Konvence v EIS

DRT se standardne zobrazuje jako gamma vs log(tau) nebo gamma vs log(f), protoze:
- Bode plot pouziva log(f)
- Nyquist ma logaritmickou strukturu frekvenci
- Literatura pouziva log-skalu

GMM v log(tau) tedy pracuje ve stejnem prostoru, ve kterem se DRT interpretuje.

### Shrnuti transformace

| Vlastnost | Linearni tau | Log(tau) |
|-----------|--------------|----------|
| Rozsah | Stlaceny | Rovnomerny |
| Tvar piku | Asymetricky | Symetricky = Gaussovsky |
| Sirka sigma | Zavisi na poloze | Univerzalni (v dekadach) |
| Separace | Nekonzistentni | Konzistentni |
| Konvence | Nepouziva se | Standard v EIS |

Proto kod dela:
```python
log_tau = np.log10(tau)  # Transformace do log-prostoru
# ... GMM fit v log_tau ...
tau_center = 10**mu      # Zpet do linearniho prostoru pro vystup
```

## Bayesovsky problem: Kolik Gaussianu?

Tady prichazi **BIC (Bayesian Information Criterion)**. Intuitivne:

```
BIC = "Jak moc se model neshoduje s daty" + "Penalizace za slozitost"
```

Formalne:
```
BIC = -2 * ln(likelihood) + k * ln(n)
```

kde:
- `k` = pocet parametru modelu
- `n` = pocet dat

## Proc je to Bayesovske?

BIC vychazi z Bayesova teoremu. Zjednodusene:

```
P(Model | Data)  ~  P(Data | Model) x P(Model)
```

- **P(Data | Model)** = likelihood - jak dobre model vysvetluje data
- **P(Model)** = prior - a priori preference jednodussich modelu

BIC je aproximace log(P(Model | Data)). Nizsi BIC = pravdepodobnejsi model.

## Occamova britva v akci

Klicova intuice: **Pridani dalsiho Gaussianu vzdy zlepsi fit, ale ne vzdy zlepsi model.**

Priklad:

| Model | Pocet piku | Fit kvalita | Penalizace | BIC |
|-------|------------|-------------|------------|-----|
| M1 | 1 | 85% | nizka | 150 |
| M2 | 2 | 95% | stredni | 120 |
| M3 | 3 | 97% | vysoka | 125 |
| M4 | 4 | 98% | velmi vysoka | 140 |

M2 vyhrava, protoze pridani 3. piku nezlepsi fit dostatecne na to, aby to vyvazilo penalizaci.

## Early stopping v teto implementaci

Kod pouziva navic **threshold** (vychozi 10.0):

```python
if BIC_improvement > threshold:
    best_idx = curr_idx  # Prijmi slozitejsi model
else:
    break  # Zastavit - zlepseni je prilis male
```

Intuice: "Pridej dalsi pik jen pokud to BIC zlepsi o vic nez 10 bodu."

- **Nizsi threshold (2-5)** -> citlivejsi, najde vice piku
- **Vyssi threshold (15-20)** -> konzervativnejsi, jen jasne piky

CLI parametr: `--gmm-bic-threshold`

## Prakticky priklad

Predstav si EIS mereni s dvema RC clanky:

1. Fitujes 1 Gaussian -> BIC = 200 (spatny fit)
2. Fitujes 2 Gaussiany -> BIC = 120 (vyrazne zlepseni o 80!)
3. Fitujes 3 Gaussiany -> BIC = 115 (zlepseni jen o 5)

S threshold=10 se zastavi na 2 picich, protoze 80 > 10 ale 5 < 10.

## Shrnuti

```
         +-------------------------------------+
         |    Bayesovske rozhodovani = BIC     |
         +-------------------------------------+
                          |
         +----------------+----------------+
         |                                 |
         v                                 v
   Likelihood                         Penalizace
 (fit accuracy)                    (model complexity)
         |                                 |
         +----------------+----------------+
                          |
                          v
              "Optimalni pocet piku"
              = balance mezi
              presnosti a jednoduchosti
```

Chces model, ktery:

- Vysvetli data dobre
- Neni zbytecne slozity
- Je robustni vuci sumu

BIC tohle dela automaticky - proto je to objektivni a reprodukovatelne.

## Souvisejici dokumentace

- [GMM Peak Detection](GMM_PEAK_DETECTION.md) - technicke detaily implementace a pouziti
- [GCV Implementation](GCV_IMPLEMENTATION.md) - volba regularizacniho parametru lambda
