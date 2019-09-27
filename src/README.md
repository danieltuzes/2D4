# 2D4::2D_DDD_simulation
This is the main component of 2D4, that evolvs discrete dislocation system according to the interaction between the single particles.

At this moment, all the information on how to use this component can be found in the [main README file](https://github.com/danieltuzes/2D4) for the complete repo.

## Roadmap
Fejleszteni lehet a kód hatékonyságát tudományos eszközök segítségével, úgy mint cache hit és matematikai azonosságok használatával. Ezen túl kell fordító szintjén is a fejlesztés, valamint új funkció implementálása, amihez meg kell érteni a kód egy részét, és egyéni szájíz szerint át kell írni. A tervezett lépések:

1. Gábor sdddst kódját github másolni és futtatni 64-es és 1024-es méretekre, 1-1-et
2. A saját szíjíz szerint átírtat, új funkciót nem tartalmazót is lefuttatni ugyanazokra
3. A compiler optimalizált, nativ march-olt kódot előállítani (assert-tel tesztelni lehet, hogy Release-s)
4. Ha nem ugyanaz az 1-es és a 3-as, akkor kideríteni az okát, egyébként sebességkülönbséget mérni 
   
   1. unordered esetre
   2. mennyit gyorsít, ha a diszlokációk rendezettek?

5. Új funkciókat implementálni

   1. dconf kiírásnál a Burgers vector értéke legyen cask 1 v -1, nem kell fixed scientific
   2. logfile-hoz nfo file-t írni
   3. feszültségérték kezdeti értéke lehessen konstans
   3. 