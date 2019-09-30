# 2D4::2D_DDD_simulation
This is the main component of 2D4, that evolvs discrete dislocation system according to the interaction between the single particles.

At this moment, all the information on how to use this component can be found in the [main README file](https://github.com/danieltuzes/2D4) for the complete repo.

## Roadmap
Fejleszteni lehet a kód hatékonyságát tudományos eszközök segítségével, úgy mint cache hit és matematikai azonosságok használatával. Ezen túl kell fordító szintjén is a fejlesztés, valamint új funkció implementálása, amihez meg kell érteni a kód egy részét, és egyéni szájíz szerint át kell írni. A tervezett lépések:

1. Gábor sdddst kódját github másolni és futtatni 64-es és 1024-es méretekre, 1-1-et
2. A saját szíjíz szerint átírtat (de std::pow-val), új funkciót nem tartalmazót is lefuttatni ugyanazokra (hiba javítva)
3. Sebességtesztet csinálni, és megmérni, mennyit számít

   1. a Release kapcsoló?
   1. az O3 march=native kapcsoló?
   2. az ffast-math
   4. unordered - ordered diszlokáció
   3. az f_dx-ben okosan zárójelezett hatványok

   Az eredmények azt mutatják, hogy 1-3. gyorsít, 4. és 5. nem gyorsít a programon 64-es rendszerméret mellett.

    | #           | eredeti, nem release   | + Release              | + O3 native            | + fast-math            | + ordered              | (+ hatv átcsop)        | (+ if(true) x)         |
    |-------------|------------------------|------------------------|------------------------|------------------------|------------------------|------------------------|------------------------|
    | 1. futtatás | 2.0118000030517578e+01 | 2.0917999982833862e+01 | 1.5968999862670898e+01 | 1.1220000028610229e+01 | 1.0274999856948853e+01 | 1.1256000041961670e+01 | 1.0776000022888184e+01 |
    | 2. futtatás | 2.0361999988555908e+01 | 2.0230999946594238e+01 | 1.3758000135421753e+01 | 1.0907999992370605e+01 | 1.0813000202178955e+01 | 1.0667000055313110e+01 | 1.0984999895095825e+01 |
    | 3. futtatás | 2.0273999929428101e+01 | 1.9748000144958496e+01 | 1.4757000207901001e+01 | 1.0744999885559082e+01 | 1.0612999916076660e+01 | 1.0647000074386597e+01 | 1.0644999980926514e+01 |

4. Forráskódot egyszerűsíteni és szépíteni
   1. a *Simulation* és *SimulationData* class public függvényeinek `const` minősítéseit helyesen kiírni
   2. megjegyzéseket írni, hogy melyik változó és függvény micsoda

4. Új funkciókat implementálni

   1. dconf kiírásnál a Burgers vector értéke legyen csak 1 v -1, nem kell fixed scientific
   2. **✓** logfile-hoz headert írni és kerüljön bele az eltelt számítógépes idő 
   3. feszültségérték kezdeti értéke lehessen konstans
   4. a diszlokációkat beolvasásnál rendezze
   5. a diszlokációk burgersvectorát az ID-ből származtassa, sebességteszt 64-es, 1024-es és 16384-as rendszerméretre
   6. blokkosított módon iterálni végig a diszlokációkon, sebességteszt 64-es, 1024-es és 16384-as rendszerméretre
   7. sinh és cosh egyszerűsített számolása, sebességteszt 64-es, 1024-es és 16384-as rendszerméretre