# Speedtests
Itt különféle módszerek sebességtesztjeit írom le a [2D4_sim](README.md), a [2D4](../README.md) főprogramjához kapcsolódóan.

### Mikor jó egy újítás
Sebességtesztek azért kellenek, hogy lehessen látni, hogy egyes változtatások valóban gyorsítanak-e a programon, vagy csak bonyolultabb és nehezebben átlátható kódot eredményeznek-e. Ha nincs jelentős sebességnövekedés és a kód sem egyszerűsödik, akkor a változtatásokat általában nem hajtom végre.

### Mennyire megbízható egy sebességteszt
A futási sebesség függ a számítógép lelki állapotától is, ami egy véletlen faktor, ezért többször kell mérni adott realizáció mellett is. Azonban különböző realizációk is szükségesek, illetve különböző rendszerméretek, illetve paraméterek. Ez olyan sok lehetőség, hogy legtöbbször csak megnézem 1 paraméter-együttesre, egyetlen realizációra és rendszerméretre. Ha itt van előnye az új kódnak, akkor jöhet, ha nem, akkor kuka.

## Mérések
Most pedig jöjjön a felsorolása a különféle méréseknek.

### Első ötletek <a name="első-ötletek"></a>

Már a kód legelső változata során eszembe jutott néhány fejlesztési lehetőség, amiket elsőkörben próbáltam ki. Ezek (amelyiket megtartottam, ott **✓** van):
   1. **✓ Release** kapcsoló: a fordítónak nincs konkrétan beadva a `Release` kapcsoló. cmake-ben ez nem is triviális.
   2. **✓ O3 march=native** kapcsoló: ez sincs beadva. Az `O3` csak gyorsít, a `march=native` viszont már platformfüggővé is teszi a kódot. Azon a gépen garantált csak a futása, amin fordítva volt.
   3. **✓ ffast-math** kapcsoló: float-okkal sok művelet nem asszociatív és disztributív úgy, mint ahogy a matekban megszoktuk. Pl. `1./a * b != b/a`, de ha be van kapcsolva a kapcsoló, akkor átrendezi a fordító, továbbá az `std::pow(double, int)`-et is szorzásokra bontja le. Ez a kapcsoló megváltoztatja a koordináták helyeit `10^{-8}` pontosság környékén.
   4. **✓ unordered - ordered**: a diszlokációkat rendezve tárolja, csökkenő Burger's vektor és csökkenő y érték szerint. Ez javítja a tippet a Burger's vektoron, és a memóriát is kicsit lokálisabbá teszi, mert nagyobb eséllyel hatnak kölcsön egymással y-ban is közeli részecskék. Segít továbbá abban, hogy később a diszlokáció Burgers' vektorát a sorszámából származtassa.
   5. az f_dx-ben okosan zárójelezett hatványok: nem kell külön elvégezni a `pow(x,6)`-ot, ha már a `pow(x,3)` el van végezve. Az `ffast-math` kapcsoló használata esetén úgy tűnik, a fordító okosabb nálam.
   6. **✓** Felesleges ellenőrizni, hogy a diszlokációk 1e-6-nál közelebb vannak-e, és az alapján számolni a kölcsönhatást, mert úgyis csak a derivált térben számít igazán. Így igaz, de a sebességteszt alapján egy `if` elágazás nem sokat lassít a CPU kódon, mert nem párhuzamos. Így viszont legalább korrekt a kód.

      Ugyanakkor a ANALYTIC_FIELD_N darab további tér számolásánál tényleg nem kell ellenőrizni. Nemcsak egy `if`-et spórol meg az ember, de két négyzetreemelést is, illetve a feltételt nem ellenőrző kód egyszerűbb (viszont ennyivel bővül a forráskód). Ha `USE_IEEE_HYPERBOLIC 0`, akkor nem ellenőrzi a távoli terekre a feltételt.

Az eredményeket táblázatba foglaltam, az eredmények másodpercben értendők, 64-es rendszerméret mellett, t=100 értékig futtatva.


| futtatás  | eredeti, nem release   | + Release (1.)         | + O3 native (2.)       | + fast-math (3.)       | + ordered (4.)         | (+ hatv átcsop) (5.)   | (+ if(true) x) (6.)    |
|:-----------:|:----------------------:|:----------------------:|:----------------------:|:----------------------:|:----------------------:|:----------------------:|:----------------------:|
| 1.        | 2.0118000030517578e+01 | 2.0917999982833862e+01 | 1.5968999862670898e+01 | 1.1220000028610229e+01 | 1.0274999856948853e+01 | 1.1256000041961670e+01 | 1.0776000022888184e+01 |
| 2.        | 2.0361999988555908e+01 | 2.0230999946594238e+01 | 1.3758000135421753e+01 | 1.0907999992370605e+01 | 1.0813000202178955e+01 | 1.0667000055313110e+01 | 1.0984999895095825e+01 |
| 3.        | 2.0273999929428101e+01 | 1.9748000144958496e+01 | 1.4757000207901001e+01 | 1.0744999885559082e+01 | 1.0612999916076660e+01 | 1.0647000074386597e+01 | 1.0644999980926514e+01 |

### Normalize <a name="normalize"></a>
Vajon egy árva while - if felcserélése segíthet? Történt ugyanis, hogy az `utility.h` fileban a `normalize(double&)` függvényben egy `while` van a szükséges `if` helyett. Ez lassít a programon, de mennyire? Semennyire. Alant a mérések.

|    index   |        1 |        2 |        3 |        4 |        5 |        6 |        7 |        8 |        9 |       10 |       11 |       12 |       13 |
|:----------:|:--------:|:--------:|:--------:|:--------:|:--------:|:---------:|:-------:|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
|  loop-pal  | 11.97199 | 11.14700 | 11.09599 | 9.999000 | 11.09400 | 11.22399 | 11.04999 | 11.02100 | 11.08200 | 11.11500 | 11.05200 | 11.15100 | 11.16000 |
| 1db if-fel | 11.35199 | 10.95499 | 11.05699 | 10.97300 | 10.93000 | 11.10100 | 11.09800 | 11.05700 | 11.01399 | 11.02399 | 10.98600 | 11.24300 | 11.23699 |
		

## További ötletek

### Nem IEEE kompatibilis trigonometrikus függvények <a name="ieee_hyperbolic"></a>
Ebben az esetben a sinh(x+i) és cosh(x+i) függvények lineárkombinációinak i-re való felösszegzésében felhasználtam a matematikai azonosságukat, így nem kell `ANALYTIC_FIELD_N` * 2 alkalommal kiszámolni őket. Ugyanis a cosh(x) és sinh(x) értékeiből konstansok szorzásával és összeadásával elérhető a kívánt eredmény. Eltérés ugyanakkor van a pontosságban, ami egy numerikus zajként jelentkezik. Ha ezzel játszanál, tekintsd meg kérlek a [játszótéren](../sandbox/ieee_hyperbolic/readme.md)!

| rendszerméret (mód) |        1. teszt        |        2. teszt        |        3. teszt        |        4. teszt        |        5. teszt        |
|:-------------------:|:----------------------:|:----------------------:|:----------------------:|:----------------------:|:----------------------:|
|    64 (azonosság)   | 5.2159998416900635e+00 | 5.3910000324249268e+00 | 5.1890001296997070e+00 | 5.2009999752044678e+00 | 5.3299999237060547e+00 |
|   64 (IEEE komp.)   | 1.0894000053405762e+01 | 1.0714999914169312e+01 | 1.0882000207901001e+01 | 1.1835000038146973e+01 | 1.0898999929428101e+01 |
|   1024 (azonosság)  | 1.5149408999919891e+04 | 1.3689345000028610e+04 | 1.5099381999969482e+04 | 1.5056825999975204e+04 | 1.4902827000141144e+04 |
|  1024 (IEEE komp.)  | 2.4513041000127792e+04 | 2.4512913999795914e+04 | 2.4554437000036240e+04 | 2.4696273000001907e+04 | 2.4696167000055313e+04 |
|  16384 (azonosság)  | 2.1527989997863770e+03 | 2.1368240001201630e+03 | 2.0477929999828339e+03 |                        |                        |
|  16384 (IEEE komp.) | 3.7827619998455048e+03 | 3.7910530002117157e+03 | 3.4341679999828339e+03 |                        |                        |

### Diszlokációk tárolása Burgers' vektor nélkül
Ha a diszlokáció típusát a sorszámából tudja származtatni, akkor nem kell eltárolni a Burgers' vektorát. Ez kevesebb memóriát igényel, viszont egy számolást (int összehasonlítás egy konstans számmal) el kell végezni. Nagy rendszerméretekre ez mindenképp előny kell(ene), hogy legyen, mert ott a sebességet a memóriaelérés limitálhatja. GPU-n még valószínűbb, hogy ez egy hasznos újítás.


A mérési **eredmények** az alábbi táblázatban láthatóak. A Burgers' vektorral sorban a program eltárolja a Burgers' vektort, az egyszerűsítettben kiszámolja. Mindkét esetben a cosh és a sinh különböző argumentumokra azonossággal van számolva. A program csak `-O3 -native=march -ffast-math` sebességoptimalizálással rendelkezik. Az eredmények szerint a várttal ellentétben lényegében nem gyorsított a programon, mi több, inkább lassított. Talán hiányzott csak az `-flto` kapcsoló? Talán nem inline-oltam? [Utóbbi biztos nem](https://stackoverflow.com/questions/1759300/when-should-i-write-the-keyword-inline-for-a-function-method).

|     rendszerméret (mód)    |        1. teszt        |        2. teszt        |        3. teszt        |        4. teszt        |        5. teszt        |
|:--------------------------:|:----------------------:|:----------------------:|:----------------------:|:----------------------:|:----------------------:|
|   64 (Burgers' vektorral)  | 5.3080000877380371e+00 | 5.3410000801086426e+00 | 5.2910001277923584e+00 | 5.2750000953674316e+00 | 5.2620000839233398e+00 |
|     64 (egyszerűsített)    | 5.3159999847412109e+00 | 5.3970000743865967e+00 | 5.3099999427795410e+00 | 5.3910000324249268e+00 | 5.3429999351501465e+00 |
|  1024 (Burgers' vektorral) | 1.5731082999944687e+04 | 1.5780119999885559e+04 | 1.7698270999908447e+04 | 1.7703957999944687e+04 | 1.7725728999853134e+04 |
|    1024 (egyszerűsített)   | 1.5880330999851227e+04 | 1.6714176000118256e+04 | 1.7654266999959946e+04 | 1.7829569000005722e+04 | 1.7845421000003815e+04 |
| 16384 (Burgers' vektorral) | 2.0485009999275208e+03 | 2.0983800001144409e+03 | 2.1508740000724792e+03 |                        |                        |
| 16384 (egyszerűsített)     | 2.0691749999523163e+03 | 2.0900899999141693e+03 | 2.1879570000171661e+03 |                        |                        |

### `-Ofast -flto`
A kapcsoló gyorsítja a programot? A gyorsulás a kicsi, közepes és nagy rendszerekre 3%, 0% és 9% adódik. Ezek nem hangsúlyos gyorsítások, de a kapcsolókat valószínűleg érdemes ebben az állapotban tartani.

|     rendszerméret (mód)     |        1. teszt        |        2. teszt        |        3. teszt        | 4. teszt               | 5. teszt               |
|:---------------------------:|:----------------------:|:----------------------:|:----------------------:|------------------------|------------------------|
|   64 (O3 fast-math native)  | 5.9229998588562012e+00 | 5.7869999408721924e+00 | 5.7030000686645508e+00 | 5.5079998970031738e+00 | 5.5789999961853027e+00 |
|    64 (Ofast native flto)   | 5.6670000553131104e+00 | 5.5780000686645508e+00 | 5.5190000534057617e+00 | 5.4210000038146973e+00 | 5.5039999485015869e+00 |
| 1024 (O3 fast-math native)  | 1.5133529999971390e+04 | 1.5076518000125885e+04 | 1.5133529999971390e+04 | 1.4964429000139236e+04 | 1.4958436999797821e+04 |
| 1024 (Ofast native flto)    | 1.5152983999967575e+04 | 1.4931483999967575e+04 | 1.5093271000146866e+04 | 1.5103620000123978e+04 | 1.4904933000087738e+04 |
| 16384 (O3 fast-math native) | 2.0485009999275208e+03 | 2.0983800001144409e+03 | 2.1508740000724792e+03 |                        |                        |
| 16384 (Ofast native flto)   | 1.8966349999904633e+03 | 1.9895490000247955e+03 | 1.9001370000839233e+03 |                        |                        |