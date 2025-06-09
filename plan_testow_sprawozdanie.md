# PLAN TESTÓW DLA SPRAWOZDANIA

# Projekt: Sequencing by Hybridization (SBH)

## ✅ ZGODNOŚĆ Z WYMAGANIAMI

### Zaimplementowane komponenty:

- ✅ **Generator instancji** (DNAGenerator + SpectrumGenerator)
- ✅ **Generator rozwiązania początkowego** (ClassicSBH)
- ✅ **Algorytm metaheurystyczny** (ThreePhaseSBH - trzyfazowy adaptacyjny)
- ✅ **Miara jakości** (podobieństwo Levenshteina)
- ✅ **Testy jednostkowe** (19 testów - wszystkie przechodzą)

### Parametry zgodne z wymaganiami:

- **Długość DNA**: 300-1000 nt ✅
- **Rozmiar k-meru**: 7-10 ✅
- ✅ **Błędy**: 2-15% negatywne/pozytywne

---

## 📋 SYSTEMATYCZNY PLAN TESTÓW

### **FAZA 1: TESTY PARAMETRÓW ALGORYTMU** (Type 1)

_Cel: Uzasadnienie wybranych parametrów metaheurystyki_

#### Test 1.1: Wpływ parametru `candidate_size`

**Parametry stałe**: n=300, k=8, błędy=5%+5%, powtórzeń=5  
**Testowane wartości**: [3, 5, 8, 10, 15, 20, 25, 30]  
**Komenda**:

```bash
python test_candidate_size_custom.py --length 300 --k 8 --pos_error 0.05 --neg_error 0.05 --candidates "3,5,8,10,15,20,25,30" --repetitions 5
```

#### Test 1.2: Wpływ parametru `error_threshold`

**Parametry stałe**: n=300, k=8, błędy=10%+10%, candidate_size=15, powtórzeń=5  
**Testowane wartości**: [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]  
**Potrzebny nowy skrypt**: `test_error_threshold.py`

#### Test 1.3: Porównanie strategii algorytmów

**Parametry stałe**: n=400, k=8, błędy=8%+8%, powtórzeń=10  
**Testowane algorytmy**: ClassicSBH vs ThreePhaseSBH  
**Komenda**:

```bash
python test_two_phase.py --length 400 --k 8 --error 0.08 --trials 10
```

---

### **FAZA 2: TESTY JAKOŚCI SEKWENCJONOWANIA** (Type 2)

_Cel: Analiza wpływu parametrów problemulogienia na trudność rozwiązania_

#### Test 2.1: Wpływ długości sekwencji DNA

**Parametry stałe**: k=8, błędy=5%+5%, candidate_size=15, powtórzeń=8  
**Testowane długości**: [200, 300, 400, 500, 600, 700, 800]  
**Komendy**:

```bash
python test_candidate_size_custom.py --length 200 --k 8 --candidates "15" --repetitions 8
python test_candidate_size_custom.py --length 300 --k 8 --candidates "15" --repetitions 8
# ... (dla każdej długości)
```

#### Test 2.2: Wpływ rozmiaru k-meru

**Parametry stałe**: n=400, błędy=6%+6%, candidate_size=15, powtórzeń=8  
**Testowane k**: [6, 7, 8, 9, 10]  
**Komendy**:

```bash
python test_candidate_size_custom.py --length 400 --k 6 --candidates "15" --repetitions 8
python test_candidate_size_custom.py --length 400 --k 7 --candidates "15" --repetitions 8
# ... (dla każdego k)
```

#### Test 2.3: Wpływ poziomu błędów negatywnych

**Parametry stałe**: n=350, k=8, błędy_poz=3%, candidate_size=15, powtórzeń=6  
**Testowane błędy negatywne**: [0%, 2%, 5%, 8%, 12%, 15%]  
**Komendy**:

```bash
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.03 --neg_error 0.0 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.03 --neg_error 0.02 --candidates "15" --repetitions 6
# ... (dla każdego poziomu błędów)
```

#### Test 2.4: Wpływ poziomu błędów pozytywnych

**Parametry stałe**: n=350, k=8, błędy_neg=3%, candidate_size=15, powtórzeń=6  
**Testowane błędy pozytywne**: [0%, 2%, 5%, 8%, 12%, 15%]

#### Test 2.5: Wpływ kombinacji błędów

**Parametry stałe**: n=300, k=8, candidate_size=15, powtórzeń=4  
**Testowane kombinacje**:

- (0%, 0%), (2%, 2%), (5%, 5%), (8%, 8%), (10%, 10%), (15%, 15%)
- (0%, 10%), (5%, 10%), (10%, 5%), (10%, 0%) - asymetryczne

---

### **FAZA 3: TESTY WYDAJNOŚCIOWE**

_Cel: Analiza czasu wykonania i skalowalności_

#### Test 3.1: Skalowalność względem długości DNA

**Parametry**: k=8, błędy=5%+5%, candidate_size=10, powtórzeń=3  
**Długości**: [100, 200, 400, 600, 800, 1000]  
**Miara**: czas wykonania + podobieństwo Levenshteina

#### Test 3.2: Skalowalność względem candidate_size

**Parametry**: n=500, k=8, błędy=8%+8%, powtórzeń=3  
**Candidate sizes**: [5, 10, 20, 40, 80, 160]

---

## 🎯 KONKRETNE KOMENDY DO WYKONANIA

### Dzień 1: Testy parametrów algorytmu

```bash
# Test 1.1 - candidate_size
python test_candidate_size_custom.py --length 300 --k 8 --pos_error 0.05 --neg_error 0.05 --candidates "3,5,8,10,15,20,25,30" --repetitions 5

# Test 1.3 - porównanie algorytmów
python test_two_phase.py --length 400 --k 8 --error 0.08 --trials 10
```

### Dzień 2: Wpływ długości i k-meru

```bash
# Test 2.1 - długość DNA (wybierz kilka reprezentatywnych)
python test_candidate_size_custom.py --length 200 --k 8 --candidates "15" --repetitions 8
python test_candidate_size_custom.py --length 400 --k 8 --candidates "15" --repetitions 8
python test_candidate_size_custom.py --length 600 --k 8 --candidates "15" --repetitions 8
python test_candidate_size_custom.py --length 800 --k 8 --candidates "15" --repetitions 8

# Test 2.2 - rozmiar k-meru
python test_candidate_size_custom.py --length 400 --k 6 --candidates "15" --repetitions 8
python test_candidate_size_custom.py --length 400 --k 7 --candidates "15" --repetitions 8
python test_candidate_size_custom.py --length 400 --k 9 --candidates "15" --repetitions 8
python test_candidate_size_custom.py --length 400 --k 10 --candidates "15" --repetitions 8
```

### Dzień 3: Wpływ błędów

```bash
# Test 2.3 - błędy negatywne
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.03 --neg_error 0.0 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.03 --neg_error 0.02 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.03 --neg_error 0.05 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.03 --neg_error 0.08 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.03 --neg_error 0.12 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.03 --neg_error 0.15 --candidates "15" --repetitions 6

# Test 2.4 - błędy pozytywne
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.0 --neg_error 0.03 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.02 --neg_error 0.03 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.05 --neg_error 0.03 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.08 --neg_error 0.03 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.12 --neg_error 0.03 --candidates "15" --repetitions 6
python test_candidate_size_custom.py --length 350 --k 8 --pos_error 0.15 --neg_error 0.03 --candidates "15" --repetitions 6
```

---

## 📊 ANALIZA WYNIKÓW

### Dla każdego testu zbierz:

1. **Średnie podobieństwo Levenshteina** (główna miara)
2. **Odchylenie standardowe** (stabilność)
3. **Czas wykonania** (wydajność)
4. **Najlepszy i najgorszy wynik** (zakres)

### Wykresy do utworzenia:

1. **Wykres liniowy**: candidate_size vs podobieństwo
2. **Wykres liniowy**: długość DNA vs podobieństwo
3. **Wykres liniowy**: poziom błędów vs podobieństwo
4. **Wykres słupkowy**: porównanie algorytmów
5. **Mapa ciepła**: kombinacje błędów neg/poz
6. **Wykres czasowy**: skalowalność

### Tabele podsumowujące:

- Najlepsze parametry dla każdego typu testu
- Porównanie algorytmów w różnych scenariuszach
- Rekomendacje użycia

---

## 🎯 KRYTERIA SUKCESU

### Minimalne wymagania:

- ✅ Przynajmniej 3 testy parametrów algorytmu
- ✅ Przynajmniej 4 testy jakości sekwencjonowania
- ✅ Użycie podobieństwa Levenshteina
- ✅ Testy na różnych długościach DNA (300-800)
- ✅ Testy na różnych poziomach błędów (2-15%)

### Docelowe wyniki:

- **Bez błędów**: >95% podobieństwa
- **Błędy 5%+5%**: >40% podobieństwa
- **Błędy 10%+10%**: >25% podobieństwa
- **Błędy 15%+15%**: >15% podobieństwa

---

_Plan utworzony na podstawie wymagań projektu - można modyfikować w zależności od dostępnego czasu i zasobów obliczeniowych._
