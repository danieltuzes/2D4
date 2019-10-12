# 2D4::2D4_sim
This is the main component of 2D4, that evolvs discrete dislocation system according to the interaction between the single particles.

At this moment, all the information on how to use this component can be found in the [main README file](https://github.com/danieltuzes/2D4) for the complete repo.

## Roadmap
Fejleszteni lehet a kód hatékonyságát tudományos eszközök segítségével, úgy mint cache hit és matematikai azonosságok használatával. Ezen túl kell fordító szintjén is a fejlesztés, valamint új funkció implementálása, amihez meg kell érteni a kód egy részét, és egyéni szájíz szerint át kell írni. A tervezett lépések:

1. Gábor sdddst kódját github másolni és futtatni 64-es és 1024-es méretekre, 1-1-et. A saját szíjíz szerint átírtat (de std::pow-val), új funkciót nem tartalmazót is lefuttatni ugyanazokra, és megnézni, hogy ugyanazok-e. Két hiba volt, a precision handler 98. sorában és az analytic fieldben rossz zárójelezés volt a 32. (?) sor környékén. A `0203ff356cc3358aaebe232baaae672aaa42e307` már csak az előbbi bugot, az utána következő már azt sem tartalmazza. Az eredmények eltérnek az eredetitől az `ffast-math` kapcsoló miatt és az átzárójelezett műveletek miatt.

3. Sebességtesztet csinálni, és megmérni, mennyit számít

   1. a **Release** kapcsoló?
   1. az **O3 march=native** kapcsoló?
   2. az **ffast-math**
   4. **unordered - ordered** diszlokáció
   3. az f_dx-ben okosan zárójelezett hatványok

   Az [eredmények](speedtests.md#első-ötletek) azt mutatják, hogy 1-3. gyorsít, 4. és 5. nem gyorsít a programon 64-es rendszerméret mellett. Az első 3-at teljesítményi okokból, a 4-et kompatibilitási okokból tartom meg.

   A sebesség- és programkódoptimalizálások, különösképp az fast-math és a műveletek átzárójelezése (1/a * b helyett b/a pl.), megváltoztatja a program műküdését. Az ebből adódó eltérések **numerikus hibák**, amelyeket lehetne minimalizálni, ha szükséges volna. Az ebből adódó eltéréseket érdemes megmérni, hogy lehessen tudni, mennyire megbízhatóak az eredmények. Ennek a pontosabb méréséhez szükség volna arra, hogy garantálni lehessen, hogy két szimulációt ugyanabban az időpillanatban írjon ki a szimuláció. Enélkül csak a legközelebbieket tudjuk csak kiírni, és az időkülönbség nem lesz 0 a kicsit eltérő számolás miatt. Az időhibával együtt mért hiba nagysága mindenesetre **2e-13** nagyságrendbe esik.

4. Forráskódot egyszerűsíteni és szépíteni
   1. **✓** a `Simulation`, `SimulationData` és `StressProtocol` classok public függvényeinek `const` minősítéseit helyesen kiírni
   2. megjegyzéseket írni, hogy melyik változó és függvény micsoda

4. Új funkciókat implementálni
   1. **✓** dconf kiírásnál a Burgers vector értéke legyen csak 1 v -1, nem kell fixed scientific
   2. **✓** logfile-hoz headert írni és kerüljön bele az eltelt számítógépes időt 
   3. feszültségérték kezdeti értéke lehessen konstans, és lehessen ciklikusan terhelni
   4. Lehessen állapotokat kiíratni adott szimulációs időgyakorisággal
   5. **✓** a diszlokációkat beolvasásnál rendezze, majd kiírásnál rendezze vissza az eredeti rendbe
   6. a diszlokációk burgersvectorát az ID-ből származtassa, sebességteszt 64-es, 1024-es és 16384-as rendszerméretre
   7. sinh és cosh egyszerűsített számolása, sebességteszt 64-es, 1024-es és 16384-as rendszerméretre
   8. blokkosított módon iterálni végig a diszlokációkon, sebességteszt 64-es, 1024-es és 16384-as rendszerméretre
   9. **✓** a `normalize` nem kell loopot tartalmazzon, elég csak 1x ellenőrizni

### Új funkciók implementálása

   1. **✓** dconf kiírásnál a Burgers vector értéke legyen csak 1 v -1, nem kell fixed scientific
   
		A beolvasásnál ellenőrzi, hogy a Burger's vector kb egész-e, és ha igen, egészekre kerekíti. Ellenőrzi, hogy 0-e az össz Burger's vector. Ha bármelyikre a válasz nem, akkor kilép. Továbbá a pozíció kiírási típusa defaulton van pontosságát is csökkentettem, az utolsó jegyek úgysem érvényesek, így viszont könnyebben áttekinthető az érték.

   2. **✓** logfile-hoz headert írni és kerüljön bele az eltelt számítógépes időt

		Belekerült a header #-vel kezdve (gnuplot kihagyja a sort), és a cellahatárolók tabulátorok, sokkal könnyebb áttekinteni.

   3. feszültségérték kezdeti értéke lehessen konstans, és lehessen ciklikusan terhelni

		A kezdeti érték lehet konstans, és lineárisan növekvő is. A konstans 0 és a monoton növők ugyanazokat a diszlokációelrendeződéseket adják 64-es rendszerre, 1000 időre és 1e-4 rátájú feszültségemelésre. Lehet ciklikusan terhelni.

   4. Lehessen állapotokat kiíratni adott szimulációs időgyakorisággal
   5. **✓** a diszlokációkat beolvasásnál rendezze, majd kiírásnál rendezze vissza az eredeti rendbe

		Rendezi és elmenti a visszarendezés sorrendjét `disl_order` néven, kiírásnál pedig ennek segítségével íratja ki. Ugyanez érvényes a subconfigokra is.

   6. a diszlokációk Burger's vectorát az ID-ből származtassa, sebességteszt 64-es, 1024-es és 16384-as rendszerméretre
   7. sinh és cosh egyszerűsített számolása, sebességteszt 64-es, 1024-es és 16384-as rendszerméretre
   8. blokkosított módon iterálni végig a diszlokációkon, sebességteszt 64-es, 1024-es és 16384-as rendszerméretre
   9. **✓** A `normalize` nem kell loopot tartalmazzon, elég csak 1x ellenőrizni. A [sebességteszt](speedtests.md#normalize) szerint kb. 0.1% az előny. Így marad a `while`.

		
